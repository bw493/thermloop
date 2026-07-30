#pragma once

// <math.h> rather than <cmath>: this header is now included by AVR
// firmware, where the C header is the portable choice. HUGE_VAL comes
// from here and supplies the unbounded default limits.
#include <math.h>

#if !defined(ARDUINO_ARCH_AVR)
#include <cmath>
using namespace std;
#endif

// avr-libc does not export HUGE_VAL under a strict -std= dialect.
#if !defined(HUGE_VAL)
#define HUGE_VAL (__builtin_huge_val())
#endif

class PIController {
public:
    // Output limits default to unbounded. Existing callers -- including the
    // Tustin validation in controller_test.cpp -- keep their exact previous
    // behaviour, since the clamps can never fire at +/- HUGE_VAL.
    PIController(double Kp, double Ki, double dt,
                 double uMin = -HUGE_VAL, double uMax = HUGE_VAL) {
        Kp_ = Kp;
        Ki_ = Ki;
        dt_ = dt;
        uMin_ = uMin;
        uMax_ = uMax;

        // Tustin (bilinear) discretization of a continuous PI controller:
        //   u(t) = Kp*e(t) + Ki*integral(e)
        // maps to the difference equation:
        //   u[k] = u[k-1] + a0*e[k] + a1*e[k-1]
        a0_ = Kp_ + Ki_ * dt_ / 2.0;
        a1_ = -Kp_ + Ki_ * dt_ / 2.0;

        prevError_ = 0.0;
        prevOutput_ = 0.0;
    }

    double update(double error) { // handling the error from trap sum
        double output = prevOutput_ + a0_ * error + a1_ * prevError_;

        // Conditional integration (anti-windup).
        //
        // Clamping the RETURNED value alone is not enough. The recursion
        // carries u[k-1], so if the unclamped value were stored the loop
        // would keep accumulating while the actuator sits pinned at its
        // limit, and the controller could not respond when the error
        // finally reverses. Storing the SATURATED output is what makes
        // the loop recoverable. With the default unbounded limits this
        // branch never fires and the recursion is untouched.
        if (output > uMax_) output = uMax_;
        if (output < uMin_) output = uMin_;

        prevError_ = error;
        prevOutput_ = output;
        return output;
    }

    void setLimits(double uMin, double uMax) {
        uMin_ = uMin;
        uMax_ = uMax;
    }

    void reset() {
        prevError_ = 0.0;
        prevOutput_ = 0.0;
    }

private:
    double Kp_, Ki_, dt_;
    double a0_, a1_;
    double uMin_, uMax_;
    double prevError_, prevOutput_;
};

class ContinuousPIReference {
public:
    ContinuousPIReference(double Kp, double Ki, double dt) {
        Kp_ = Kp;
        Ki_ = Ki;
        dt_ = dt;
        integral_ = 0.0;
    }

    double update(double error) {
        // Continuous reference integrator (fine-grained ground truth).
        // Deliberately left unclamped: this class is the validation
        // baseline, not a controller meant to drive an actuator.
        integral_ += error * dt_;
        return Kp_ * error + Ki_ * integral_;
    }

private:
    double Kp_, Ki_, dt_;
    double integral_;
};
