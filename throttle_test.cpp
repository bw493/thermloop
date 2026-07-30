// throttle_test.cpp - thermloop, validation driver for throttle.h
//
// Four checks, each empirical rather than symbolic:
//   1. DC steady state          -> T_amb + P * sum(R)
//   2. Recursion vs convolution -> the two realizations agree
//   3. Closed-form alpha vs bisection -> the affine argument holds
//   4. Wearable BCI scenario    -> 1 C tissue rise budget is respected
//
// Build (desktop): g++ -std=c++17 -O2 throttle_test.cpp -o throttle_test
#include "throttle.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace std;

// ---------------------------------------------------------------------
// Ground truth 1: Fable's convolution kernel, evaluated directly.
//   h[m] = sum_i R_i * (exp(-m*dt/tau_i) - exp(-(m+1)*dt/tau_i))
// ---------------------------------------------------------------------
vector<double> buildKernel(const vector<double>& R, const vector<double>& C,
                           double dt, int N) {
    vector<double> h(N, 0.0);
    for (int m = 0; m < N; ++m) {
        double v = 0.0;
        for (size_t i = 0; i < R.size(); ++i) {
            double tau = R[i] * C[i];
            v += R[i] * (exp(-m * dt / tau) - exp(-(m + 1) * dt / tau));
        }
        h[m] = v;
    }
    return h;
}

// Convolution temperature after the power sequence P has been applied,
// starting from a cold plant. T[n] = T_amb + sum_m h[m] * P[n-m].
double convolveAt(const vector<double>& h, const vector<double>& P,
                  int n, double T_amb) {
    double T = T_amb;
    for (size_t m = 0; m < h.size(); ++m) {
        int idx = n - (int)m;
        if (idx < 0) break;
        T += h[m] * P[idx];
    }
    return T;
}

// ---------------------------------------------------------------------
// Ground truth 2: Fable's bisection on alpha. Rolls the plant forward for
// a candidate alpha and reports the peak, then binary searches.
// ---------------------------------------------------------------------
double peakForAlpha(const ThermalThrottle& base, const vector<double>& planned,
                    double alpha) {
    ThermalThrottle sim = base;          // copy carries the current state
    double peak = -1e30;
    for (size_t k = 0; k < planned.size(); ++k) {
        sim.push(alpha * planned[k]);
        double T = sim.temperature();
        if (T > peak) peak = T;
    }
    return peak;
}

double bisectAlpha(const ThermalThrottle& base, const vector<double>& planned,
                   double T_limit, double alphaMin, double tol) {
    if (peakForAlpha(base, planned, 1.0) <= T_limit) return 1.0;
    double lo = alphaMin;
    double hi = 1.0;
    while (hi - lo > tol) {
        double mid = 0.5 * (lo + hi);
        if (peakForAlpha(base, planned, mid) > T_limit) hi = mid;
        else lo = mid;
    }
    return lo;
}

// ---------------------------------------------------------------------
int main() {
    cout << fixed << setprecision(6);

    // The MCU plant from plant.h, sampled well below the smallest tau.
    vector<double> R = {2.0, 8.0, 20.0};
    vector<double> C = {0.25, 1.0, 6.0};
    double T_amb = 25.0;
    double dt = 0.05;

    ThermalThrottle plant;
    plant.begin(R.data(), C.data(), (int)R.size(), T_amb, dt);

    double sumR = 0.0;
    for (size_t i = 0; i < R.size(); ++i) sumR += R[i];

    // --- Check 1: DC steady state -------------------------------------
    cout << "[1] DC steady state" << endl;
    double P_dc = 2.0;
    plant.reset();
    for (int k = 0; k < 200000; ++k) plant.push(P_dc);
    cout << "    settled T   = " << plant.temperature() << endl;
    cout << "    expected    = " << T_amb + P_dc * sumR << endl;
    cout << "    error       = "
         << fabs(plant.temperature() - (T_amb + P_dc * sumR)) << endl << endl;

    // --- Check 2: recursion vs convolution ----------------------------
    cout << "[2] Recursion vs convolution kernel" << endl;
    double tauMax = 0.0;
    for (size_t i = 0; i < R.size(); ++i) {
        double t = R[i] * C[i];
        if (t > tauMax) tauMax = t;
    }
    int N = (int)ceil(5.0 * tauMax / dt);
    vector<double> h = buildKernel(R, C, dt, N);

    // A bursty trace: 4 W bursts against a 0.2 W idle floor.
    int steps = 4000;
    vector<double> P(steps);
    for (int k = 0; k < steps; ++k)
        P[k] = ((k / 200) % 3 == 0) ? 4.0 : 0.2;

    plant.reset();
    double maxDiff = 0.0;
    for (int k = 0; k < steps; ++k) {
        plant.push(P[k]);
        double d = fabs(plant.temperature() - convolveAt(h, P, k, T_amb));
        if (d > maxDiff) maxDiff = d;
    }
    cout << "    kernel taps = " << N << "  (" << N * sizeof(float)
         << " bytes as float)" << endl;
    cout << "    recursion   = " << R.size() << " stages  ("
         << 3 * R.size() * sizeof(float) << " bytes as float)" << endl;
    cout << "    max |diff|  = " << maxDiff << " C" << endl << endl;

    // --- Check 3: closed-form alpha vs bisection ----------------------
    cout << "[3] Closed-form alpha vs bisection" << endl;
    double T_limit = 70.0;
    int horizon = 1200;
    vector<double> planned(horizon);
    for (int k = 0; k < horizon; ++k)
        planned[k] = ((k / 100) % 2 == 0) ? 4.0 : 0.5;

    // Warm the plant so the "past" term is genuinely nonzero.
    plant.reset();
    for (int k = 0; k < 600; ++k) plant.push(1.5);
    cout << "    T at plan time = " << plant.temperature() << " C" << endl;

    double aClosed = plant.solve(planned.data(), horizon, T_limit);
    double aBisect = bisectAlpha(plant, planned, T_limit, 0.1, 1e-9);
    cout << "    closed form    = " << aClosed << endl;
    cout << "    bisection      = " << aBisect << endl;
    cout << "    difference     = " << fabs(aClosed - aBisect) << endl;
    cout << "    peak at alpha  = " << peakForAlpha(plant, planned, aClosed)
         << " C  (limit " << T_limit << ")" << endl << endl;

    // --- Check 4: wearable BCI, 1 C tissue rise budget ----------------
    cout << "[4] Wearable BCI node, 1 C rise budget" << endl;
    vector<double> Rb = {8.0, 12.0};
    vector<double> Cb = {0.5, 5.0};
    double skin = 33.0;
    double dtb = 0.1;

    ThermalThrottle pod;
    pod.begin(Rb.data(), Cb.data(), (int)Rb.size(), skin, dtb);

    // 5 mW idle telemetry, 0.6 s bursts of 120 mW classifier inference.
    int nb = 600;
    vector<double> want(nb);
    for (int k = 0; k < nb; ++k)
        want[k] = (k % 10 < 6) ? 0.120 : 0.005;

    double aBci = pod.solve(want.data(), nb, skin + 1.0);
    int blocks = pod.solveBlocks(want.data(), nb, skin + 1.0, 8);

    double meanWant = 0.0;
    for (int k = 0; k < nb; ++k) meanWant += want[k];
    meanWant /= nb;

    cout << "    requested mean = " << meanWant * 1e3 << " mW" << endl;
    cout << "    executed mean  = " << aBci * meanWant * 1e3 << " mW" << endl;
    cout << "    alpha          = " << aBci << endl;
    cout << "    DSP blocks     = " << blocks << " of 8" << endl;
    cout << "    peak rise      = "
         << peakForAlpha(pod, want, aBci) - skin << " C  (budget 1.000)"
         << endl << endl;

    // --- Check 5: generator overload vs array overload ----------------
    // The UNO cannot hold the trace, so it synthesizes each sample on
    // demand. Both paths must agree bit for bit.
    cout << "[5] Generator overload vs array overload" << endl;

    struct BurstPlan {
        double idle;
        double burst;
        int period;
        int duty;
    };
    BurstPlan bp = {0.5, 4.0, 200, 100};

    // Same schedule as check 3, produced procedurally instead of stored.
    struct Gen {
        static double at(int k, const void* ctx) {
            const BurstPlan* p = (const BurstPlan*)ctx;
            return ((k % p->period) < p->duty) ? p->burst : p->idle;
        }
    };

    vector<double> mirror(horizon);
    for (int k = 0; k < horizon; ++k) mirror[k] = Gen::at(k, &bp);

    plant.reset();
    for (int k = 0; k < 600; ++k) plant.push(1.5);

    double aArray = plant.solve(mirror.data(), horizon, T_limit);
    double aGen   = plant.solve(&Gen::at, &bp, horizon, T_limit);
    cout << "    array form     = " << aArray << endl;
    cout << "    generator form = " << aGen << endl;
    cout << "    difference     = " << fabs(aArray - aGen) << endl;
    cout << "    stored samples = " << horizon << " vs 0 ("
         << sizeof(BurstPlan) << " bytes of context)" << endl << endl;

    // --- Check 6: planning stride -------------------------------------
    // Claim: when the plan is piecewise constant on stride boundaries,
    // the coarse solve is EXACT, because every stage moves monotonically
    // within a bin so the peak lands on a bin boundary.
    cout << "[6] Planning stride (loading reduction)" << endl;
    int stride = 10;

    struct CoarseGen {
        static double at(int k, const void* ctx) {
            return Gen::at(k * 10, ctx);   // sample the plan on bin edges
        }
    };

    plant.setPlanStride(1);
    double aFine = plant.solve(&Gen::at, &bp, horizon, T_limit);
    plant.setPlanStride(stride);
    double aCoarse = plant.solve(&CoarseGen::at, &bp, horizon / stride, T_limit);
    plant.setPlanStride(1);

    cout << "    stride 1       = " << aFine
         << "   (" << horizon << " steps)" << endl;
    cout << "    stride " << stride << "      = " << aCoarse
         << "   (" << horizon / stride << " steps)" << endl;
    cout << "    difference     = " << fabs(aFine - aCoarse) << endl;
    cout << "    work reduced   = " << stride << "x" << endl << endl;

    // --- Check 7: measurement sync on a mis-specified model -----------
    // The rig is the truth; the model is deliberately wrong in BOTH the
    // gain (R low by 20 percent) and the offset (ambient high by 2 C).
    //
    // A complementary filter removes BIAS. It cannot remove a gain error,
    // because that residual is correlated with the load and oscillates
    // with it. So the honest metrics are mean error (the bias the filter
    // is meant to kill) and RMS error (which retains the oscillation).
    cout << "[7] Measurement sync vs open-loop drift" << endl;

    vector<double> Rt = {4.0, 9.0},  Ct = {1.5, 14.0};
    vector<double> Rm = {3.2, 7.2},  Cm = {1.5, 14.0};
    double ambTrue = 24.0, ambModel = 26.0, dtr = 0.5;
    double tauSlow = Rt[1] * Ct[1];        // 126 s

    ThermalThrottle rig, openLoop, slow, fast;
    rig.begin(Rt.data(), Ct.data(), 2, ambTrue, dtr);
    openLoop.begin(Rm.data(), Cm.data(), 2, ambModel, dtr);
    slow.begin(Rm.data(), Cm.data(), 2, ambModel, dtr);
    fast.begin(Rm.data(), Cm.data(), 2, ambModel, dtr);
    slow.setSyncTau(5.0 * tauSlow);        // 630 s, well clear of the plant
    fast.setSyncTau(1.0 * tauSlow);        // 126 s, contending with it

    cout << "    slowest tau    = " << tauSlow << " s" << endl;

    int ticks = 8000;                      // 4000 s of operation
    int half = ticks / 2;                  // score the settled half only
    double sOpen = 0.0, sSlow = 0.0, sFast = 0.0;
    double qOpen = 0.0, qSlow = 0.0, qFast = 0.0;
    int scored = 0;

    for (int k = 0; k < ticks; ++k) {
        double P = ((k / 400) % 2 == 0) ? 1.2 : 0.2;
        rig.push(P);
        openLoop.push(P);
        slow.push(P);
        fast.push(P);
        slow.syncTemperature(rig.temperature());   // the thermistor
        fast.syncTemperature(rig.temperature());

        if (k >= half) {
            double T = rig.temperature();
            double eO = openLoop.temperature() - T;
            double eS = slow.temperature() - T;
            double eF = fast.temperature() - T;
            sOpen += eO;  sSlow += eS;  sFast += eF;
            qOpen += eO * eO;  qSlow += eS * eS;  qFast += eF * eF;
            ++scored;
        }
    }

    cout << "                      mean err      RMS err" << endl;
    cout << "    open loop      = " << sOpen / scored
         << "     " << sqrt(qOpen / scored) << endl;
    cout << "    sync 5x tau    = " << sSlow / scored
         << "     " << sqrt(qSlow / scored) << endl;
    cout << "    sync 1x tau    = " << sFast / scored
         << "     " << sqrt(qFast / scored) << endl;
    cout << "    bias removed   = "
         << fabs(sOpen / scored) / fabs(sSlow / scored) << "x (5x tau)" << endl;

    return 0;
}
