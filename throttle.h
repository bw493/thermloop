// throttle.h - thermloop, portable predictive thermal throttle
//
// Targets: ATmega328P (Arduino UNO Rev3), ESP32-S3, Linux/ARM (PiEEG).
// Freestanding by design: no heap, no exceptions, no STL, no iostream.
// The only header pulled in is <math.h>, for exp() at construction.
//
// Why this file does not reuse predict.h:
//   predict.h stores the convolution kernel h[m], which is 5*tau_max/dt
//   samples long. For the wearable BCI plant (tau_max = 60 s, dt = 0.1 s)
//   that is 3000 floats = 12 KB, against 2 KB of SRAM on an UNO. The
//   recursive state-space realization below is the SAME LTI system --
//   exactly, not approximately -- but costs one scalar per RC stage.
#pragma once

#include <math.h>

// ---------------------------------------------------------------------
// Platform profile. Only two things vary across the three targets: the
// scalar width and the static ceiling on stage count. Define TL_REAL and
// TL_MAX_STAGES before including this header to override.
// ---------------------------------------------------------------------
#if !defined(TL_REAL)
  #if defined(ARDUINO_ARCH_AVR)
    #define TL_REAL       float    // no FPU, 2 KB SRAM total
    #define TL_MAX_STAGES 4
  #elif defined(ARDUINO_ARCH_ESP32) || defined(ESP32)
    #define TL_REAL       float    // Xtensa LX7 FPU is single-precision only;
    #define TL_MAX_STAGES 8        // double would drop to software emulation
  #else
    #define TL_REAL       double   // Cortex-A: double costs nothing
    #define TL_MAX_STAGES 16
  #endif
#endif

typedef TL_REAL tl_real;

// expf()/powf() avoid a needless promotion on the two float targets.
#if defined(ARDUINO_ARCH_AVR) || defined(ARDUINO_ARCH_ESP32) || defined(ESP32)
  #define TL_EXP(x)    expf(x)
  #define TL_POW(x, y) powf(x, y)
#else
  #define TL_EXP(x)    exp(x)
  #define TL_POW(x, y) pow(x, y)
#endif

// Plan generator. Returns the power requested at horizon step k.
// The ctx pointer carries caller state, so no globals are required --
// which matters on an UNO, where a global would cost permanent SRAM.
typedef tl_real (*tl_plan_fn)(int k, const void* ctx);


// =====================================================================
// ThermalThrottle
//
// State: x_[i] is stage i's contribution to the rise above ambient, so
//        T[n] = T_amb + sum_i x_[i].
//
// Update (zero-order hold, exact for piecewise-constant power):
//        x_i[n+1] = a_i * x_i[n] + b_i * P[n]
//        a_i = exp(-dt / tau_i),  b_i = R_i * (1 - a_i)
//
// Unrolling that recursion reproduces the convolution kernel term by
// term, since h[m] = sum_i R_i * (1 - a_i) * a_i^m. Same system, two
// realizations; only the memory bill differs.
// =====================================================================
class ThermalThrottle {
public:
    // R in C/W, C in J/C, T_amb in C, dt in seconds.
    // Call once at setup(). This is the only place exp() is evaluated,
    // which matters on AVR where each call costs hundreds of cycles.
    void begin(const tl_real* R, const tl_real* C, int n,
               tl_real T_amb, tl_real dt) {
        n_ = (n > TL_MAX_STAGES) ? TL_MAX_STAGES : n;
        T_amb_ = T_amb;
        dt_ = dt;
        alphaMin_ = 0.1;
        stride_ = 1;
        syncGain_ = 0.0;          // sync disabled until setSyncTau()
        bias_ = 0.0;
        for (int i = 0; i < n_; ++i) {
            tl_real tau = R[i] * C[i];
            R_[i]  = R[i];
            a_[i]  = TL_EXP(-dt / tau);
            b_[i]  = R[i] * (1.0 - a_[i]);
            ac_[i] = a_[i];          // coarse poles, stride 1 by default
            bc_[i] = b_[i];
            x_[i]  = 0.0;
        }
    }

    void setAlphaMin(tl_real a) { alphaMin_ = a; }

    // -----------------------------------------------------------------
    // Planning stride: advance s ticks per horizon sample during solve().
    // This is the loading knob. Cost falls by exactly s, and the plant
    // model does NOT degrade -- advancing s ticks at constant power is
    // still exact, since x[n+s] = a^s * x[n] + R * (1 - a^s) * P. What
    // coarsens is the resolution of the PLAN, not the physics.
    //
    // The peak is still caught. Within a constant-power bin every stage
    // moves monotonically toward its own equilibrium R_i * P, so the sum
    // is monotone and its extremum sits on a bin boundary -- exactly
    // where the coarse solve samples. Verified in throttle_test.cpp.
    //
    // Call once after begin(); pow() is evaluated here, never in loop().
    // -----------------------------------------------------------------
    void setPlanStride(int s) {
        if (s < 1) s = 1;
        stride_ = s;
        for (int i = 0; i < n_; ++i) {
            ac_[i] = TL_POW(a_[i], (tl_real)s);
            bc_[i] = R_[i] * (1.0 - ac_[i]);
        }
    }

    int planStride() const { return stride_; }

    // -----------------------------------------------------------------
    // Measurement sync: an integral disturbance estimator.
    //
    // The state x_ is open loop -- push() accumulates whatever it is told,
    // so errors in R, C or T_amb persist indefinitely. The fan loop does
    // not suffer this because it closes on the sensor; the throttle's
    // PREDICTION does, and a prediction drifted above the truth throttles
    // a workload that was never in danger.
    //
    // The correction is a SEPARATE, non-decaying offset rather than a
    // nudge to the stage states. That distinction is the whole design:
    // a correction injected into x_ decays at the plant's own poles, so
    // its equilibrium loop gain is only g / (1 - a_i) -- far below unity
    // for a slow stage, leaving most of the bias uncorrected. A standing
    // offset integrates to whatever value drives the mean residual to
    // zero, and leaves the physical dynamics untouched.
    //
    // What it CANNOT do is fix a gain error. That residual is correlated
    // with the load and oscillates with it; only re-identifying R and C
    // removes it. See identify.cpp.
    //
    // Call once per control tick, AFTER push().
    // -----------------------------------------------------------------
    void syncTemperature(tl_real T_meas) {
        if (syncGain_ <= 0.0) return;
        bias_ += syncGain_ * (T_meas - temperature());
    }

    // Sync time constant, in seconds: the offset converges with exactly
    // this time constant, since g = 1 - exp(-dt / tauSync). Specified in
    // physical units rather than as a bare gain. Choose it SEVERAL TIMES
    // the slowest plant tau -- fast enough to track ambient drift, slow
    // enough not to chase sensor noise into the horizon. Zero disables.
    void setSyncTau(tl_real tauSync) {
        if (tauSync <= 0.0) {
            syncGain_ = 0.0;
            return;
        }
        syncGain_ = 1.0 - TL_EXP(-dt_ / tauSync);
    }

    tl_real syncGain() const { return syncGain_; }

    void reset() {
        for (int i = 0; i < n_; ++i) x_[i] = 0.0;
        bias_ = 0.0;
    }

    // Record one executed power sample. Call once per control tick.
    // Always advances a single dt, independent of the planning stride.
    void push(tl_real P) {
        for (int i = 0; i < n_; ++i)
            x_[i] = a_[i] * x_[i] + b_[i] * P;
    }

    tl_real temperature() const {
        tl_real T = T_amb_ + bias_;
        for (int i = 0; i < n_; ++i) T += x_[i];
        return T;
    }

    tl_real bias() const { return bias_; }

    // -----------------------------------------------------------------
    // Largest uniform scale alpha in [alphaMin, 1] such that the planned
    // trace never drives the plant past T_limit.
    //
    // Superposition removes the search entirely. Split the horizon into
    //   past[k] : free decay of the CURRENT state, planned power off
    //   fut[k]  : response to the planned trace from a ZERO state
    // Only fut scales with alpha, so
    //   T[k](alpha) = past[k] + alpha * fut[k]
    // is affine in alpha, and the binding constraint is a plain minimum:
    //   alpha = min over k of (T_limit - past[k]) / fut[k]
    //
    // One pass, no bisection, no tolerance parameter, and the result is
    // exact rather than converged. Scratch space is 2 * n_ scalars; the
    // horizon is never stored, so K may be arbitrarily long even on AVR.
    //
    // Generator form: the plan is produced on demand, so an UNO can
    // synthesize its schedule from the live focus level instead of
    // holding an array it has no room for.
    // -----------------------------------------------------------------
    tl_real solve(tl_plan_fn plan, const void* ctx, int K,
                  tl_real T_limit) const {
        tl_real xf[TL_MAX_STAGES];   // free decay from current state
        tl_real xp[TL_MAX_STAGES];   // forced response from zero state
        for (int i = 0; i < n_; ++i) {
            xf[i] = x_[i];
            xp[i] = 0.0;
        }

        tl_real alpha = 1.0;
        for (int k = 0; k < K; ++k) {
            tl_real P    = plan(k, ctx);
            tl_real past = T_amb_ + bias_;   // the offset persists ahead
            tl_real fut  = 0.0;
            for (int i = 0; i < n_; ++i) {
                xf[i] = ac_[i] * xf[i];
                xp[i] = ac_[i] * xp[i] + bc_[i] * P;
                past += xf[i];
                fut  += xp[i];
            }
            // fut > 0 whenever any power has been planned; a zero-power
            // prefix imposes no constraint and is skipped.
            if (fut > 0.0) {
                tl_real cap = (T_limit - past) / fut;
                if (cap < alpha) alpha = cap;
            }
        }

        // A negative cap means the stored state ALREADY exceeds the limit
        // on its own. Clamping to alphaMin is the correct response: throttle
        // as hard as the floor allows and let the plant decay.
        if (alpha > 1.0) alpha = 1.0;
        if (alpha < alphaMin_) alpha = alphaMin_;
        return alpha;
    }

    // Array form, for the Pi and the ESP32-S3 where the trace fits in RAM.
    // Delegates to the generator so there is exactly one implementation.
    tl_real solve(const tl_real* planned, int K, tl_real T_limit) const {
        return solve(&arrayPlan_, (const void*)planned, K, T_limit);
    }

    // -----------------------------------------------------------------
    // Discrete variants for the BCI nodes. A classifier pipeline drops
    // whole DSP blocks rather than scaling a clock, so alpha must land on
    // a k/nBlocks grid. Rounds DOWN: T_limit is a ceiling, not a target.
    // -----------------------------------------------------------------
    int solveBlocks(tl_plan_fn plan, const void* ctx, int K,
                    tl_real T_limit, int nBlocks) const {
        return quantize_(solve(plan, ctx, K, T_limit), nBlocks);
    }

    int solveBlocks(const tl_real* planned, int K,
                    tl_real T_limit, int nBlocks) const {
        return quantize_(solve(planned, K, T_limit), nBlocks);
    }

    int stages() const { return n_; }
    tl_real ambient() const { return T_amb_; }
    tl_real interval() const { return dt_; }

private:
    static tl_real arrayPlan_(int k, const void* ctx) {
        return ((const tl_real*)ctx)[k];
    }

    static int quantize_(tl_real alpha, int nBlocks) {
        int k = (int)(alpha * nBlocks);
        if (k > nBlocks) k = nBlocks;
        if (k < 0) k = 0;
        return k;
    }

    tl_real R_[TL_MAX_STAGES];   // stage resistances, kept for restride
    tl_real a_[TL_MAX_STAGES];   // poles exp(-dt/tau_i)
    tl_real b_[TL_MAX_STAGES];   // input weights R_i * (1 - a_i)
    tl_real ac_[TL_MAX_STAGES];  // coarse poles a_i^stride
    tl_real bc_[TL_MAX_STAGES];  // coarse input weights
    tl_real x_[TL_MAX_STAGES];   // per-stage rise above ambient
    tl_real T_amb_;
    tl_real dt_;
    tl_real alphaMin_;
    tl_real syncGain_;
    tl_real bias_;               // learned offset, does NOT decay
    int n_;
    int stride_;
};
