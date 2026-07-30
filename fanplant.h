// fanplant.h - thermloop, fan actuator model
//
// Gives the LAST Foster stage -- the convection stage, the one that
// touches ambient -- a controllable gain. Airflow does not change the
// silicon or the sink; it changes how easily heat leaves the surface,
// which is exactly the convection resistance:
//
//     R_conv(u) = R_conv0 / (1 + kFan * u),   u in [0, 1]
//
// u = 0 is still air (the plant identify.cpp fitted, fan OFF), and
// kFan is the fan's authority: kFan = 3 means full fan cuts the
// convection resistance to a quarter of its still-air value.
//
// Composition is by pointer, same as ThermalImpedance: FosterNetwork
// has no default constructor, so this class holds a const pointer to
// the still-air base plant and derives everything from it. C_conv is
// recovered as tau_N / R_N, so plant.h needs no new accessor.
#pragma once
#include "plant.h"
#include <complex>
#include <vector>
#include <cstddef>

class FanPlant {
public:
    FanPlant(const FosterNetwork& base, double kFan)
        : base_(&base), kFan_(kFan) {}

    // Convection resistance at fan command u.
    double convResistance(double u) const {
        const std::vector<double>& R = base_->resistances();
        double R0 = R[R.size() - 1];
        return R0 / (1.0 + kFan_ * u);
    }

    // Total DC resistance with the fan at u: inner stages unchanged,
    // convection stage squeezed by the airflow.
    double totalResistance(double u) const {
        const std::vector<double>& R = base_->resistances();
        double sum = 0.0;
        for (std::size_t i = 0; i + 1 < R.size(); ++i)
            sum += R[i];
        return sum + convResistance(u);
    }

    // Where the plant settles under constant power P and fan command u.
    double steadyTemperature(double P, double u) const {
        return base_->ambient() + P * totalResistance(u);
    }

    // ---- small-signal linearization about an operating point (P0, u0) ----
    //
    // With steady power P0 flowing out through the convection stage, the
    // temperature drop across that stage is P0 * R_conv(u). Differentiate
    // in u:
    //
    //     dT/du = P0 * dR_conv/du = -P0 * kFan * R0 / (1 + kFan*u0)^2
    //
    // More fan -> cooler, so the raw slope is negative. gain() returns the
    // MAGNITUDE (degrees C of cooling authority per unit of fan command);
    // the loop error convention T_meas - T_target, as used in
    // focus_throttle.ino, folds the sign back in: hotter than target ->
    // positive error -> fan up -> temperature down. Negative feedback.
    double gain(double P0, double u0) const {
        const std::vector<double>& R = base_->resistances();
        double R0 = R[R.size() - 1];
        double denom = 1.0 + kFan_ * u0;
        return P0 * kFan_ * R0 / (denom * denom);
    }

    // The disturbance travels through the convection stage's own RC pair.
    // C_conv is fixed hardware; R_conv shrinks with airflow, so the stage
    // gets FASTER as the fan spins up:  tau(u0) = R_conv(u0) * C_conv.
    double timeConstant(double u0) const {
        const std::vector<double>& R   = base_->resistances();
        const std::vector<double>& tau = base_->taus();
        std::size_t last = R.size() - 1;
        double C_conv = tau[last] / R[last];
        return convResistance(u0) * C_conv;
    }

    // H_fan(s): first-order transfer function from fan command to
    // temperature deviation, linearized at (P0, u0). This is the plant
    // the loop-shaping analysis composes with C(s).
    std::complex<double> response(double omega, double P0, double u0) const {
        std::complex<double> s(0.0, omega);
        double K   = gain(P0, u0);
        double tau = timeConstant(u0);
        return K / (1.0 + s * tau);
    }

    const FosterNetwork& base() const { return *base_; }
    double authority() const { return kFan_; }

private:
    const FosterNetwork* base_;
    double kFan_;
};
