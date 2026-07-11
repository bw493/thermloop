#pragma once

#include <cmath>

using namespace std;

class PIController {
public:
    PIController(double Kp, double Ki, double dt) {
        Kp_ = Kp;
        Ki_ = Ki;
        dt_ = dt;

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
        prevError_ = error;
        prevOutput_ = output;
        return output;
    }

private:
    double Kp_, Ki_, dt_;
    double a0_, a1_;
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
        // Continuous reference integrator (fine-grained ground truth)
        integral_ += error * dt_;
        return Kp_ * error + Ki_ * integral_;
    }

private:
    double Kp_, Ki_, dt_;
    double integral_;
};