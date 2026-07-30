// loopshape.h - thermloop, loop shaping for the fan control loop
//
// Composes the controller with the fan-side plant and asks the classical
// frequency-domain questions:
//
//     L(s) = C(s) * H_fan(s)
//     C(s) = Kp + Ki / s              (the same PI controller.h discretizes)
//
//   * crossoverFrequency: where |L(j*omega)| falls through 1
//   * phaseMarginDegrees: 180 + arg(L) at that crossover
//   * suggestPI:          invert the problem -- pick a crossover and a
//                         phase margin, get the (Kp, Ki) that deliver them
//
// The design loop this enables: identify.cpp fits the plant, fanplant.h
// linearizes it, suggestPI() picks gains, controller_test.cpp's Tustin
// machinery discretizes them, and focus_throttle.ino runs them onboard.
#pragma once
#include "fanplant.h"
#include <complex>
#include <cmath>

class LoopShaper {
public:
    struct PIGains {
        double Kp;
        double Ki;
    };

    // Analysis is small-signal, so the operating point is part of the
    // shaper's identity: same plant, different (P0, u0), different loop.
    LoopShaper(const FanPlant& plant, double P0, double u0)
        : plant_(&plant), P0_(P0), u0_(u0) {}

    // C(j*omega) = Kp - j * Ki / omega
    std::complex<double> controllerResponse(double Kp, double Ki,
                                            double omega) const {
        return std::complex<double>(Kp, -Ki / omega);
    }

    // L(j*omega) = C(j*omega) * H_fan(j*omega)
    std::complex<double> loopResponse(double Kp, double Ki,
                                      double omega) const {
        return controllerResponse(Kp, Ki, omega)
             * plant_->response(omega, P0_, u0_);
    }

    // Gain crossover: |L(j*omega)| = 1. For PI * first-order-lag the
    // magnitude is strictly decreasing in omega (integrator falls, lag
    // falls, nothing rises), so a bisection between a low and a high
    // frequency bracket is safe.
    double crossoverFrequency(double Kp, double Ki) const {
        double lo = 1.0e-6;
        double hi = 1.0e4;
        for (int i = 0; i < 200; ++i) {
            double mid = std::sqrt(lo * hi);   // bisect in log frequency
            if (std::abs(loopResponse(Kp, Ki, mid)) > 1.0)
                lo = mid;
            else
                hi = mid;
        }
        return std::sqrt(lo * hi);
    }

    // Phase margin: how far arg(L) sits above -180 degrees at crossover.
    double phaseMarginDegrees(double Kp, double Ki) const {
        const double pi = 3.14159265358979323846;
        double wc = crossoverFrequency(Kp, Ki);
        double phase = std::arg(loopResponse(Kp, Ki, wc)) * 180.0 / pi;
        return 180.0 + phase;
    }

    // Inverse design. Given a wanted crossover omega_c and phase margin
    // PM (degrees), solve for (Kp, Ki):
    //
    //   1. The controller must supply whatever phase the plant does not:
    //        phi_C = (PM - 180) - arg(H_fan(j*omega_c))     [degrees]
    //      A PI can only LAG (phi_C in (-90, 0]); asking it to lead means
    //      the crossover is too fast for this plant and margin. In that
    //      case the ratio below would come out negative -- the caller can
    //      check Ki > 0 as the feasibility test.
    //
    //   2. Phase of PI:  arg(C) = -atan(Ki / (Kp * omega_c))
    //        =>  Ki / Kp = -omega_c * tan(phi_C)
    //
    //   3. Magnitude closes the design:  |C| * |H_fan| = 1 at omega_c
    //        =>  Kp = 1 / ( |H_fan| * sqrt(1 + (Ki/Kp)^2 / omega_c^2) )
    PIGains suggestPI(double omegaC, double pmDegrees) const {
        const double pi = 3.14159265358979323846;

        std::complex<double> H = plant_->response(omegaC, P0_, u0_);
        double plantPhase = std::arg(H) * 180.0 / pi;
        double phiC = (pmDegrees - 180.0) - plantPhase;      // degrees
        double phiCrad = phiC * pi / 180.0;

        double ratio = -omegaC * std::tan(phiCrad);          // Ki / Kp
        double scale = std::sqrt(1.0 + (ratio / omegaC) * (ratio / omegaC));
        double Kp = 1.0 / (std::abs(H) * scale);
        double Ki = ratio * Kp;

        PIGains g;
        g.Kp = Kp;
        g.Ki = Ki;
        return g;
    }

private:
    const FanPlant* plant_;
    double P0_, u0_;
};
