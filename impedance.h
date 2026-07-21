// impedance.h
#pragma once
#include "plant.h"
#include <complex>

class ThermalImpedance {
public:
    explicit ThermalImpedance(const FosterNetwork& plant)
        : plant_(&plant) {}

        std::complex<double> evaluate(std::complex<double> s) const {
        std::complex<double> Z(0.0, 0.0);
        const std::vector<double>& R = plant_->resistances();
        const std::vector<double>& tau = plant_->taus();

        for (std::size_t i = 0; i < R.size(); ++i) {
            // accumulate this stage's contribution to Z.
            // Each term is R[i] / (1.0 + s * tau[i])
            Z += R[i] / (1.0 + s * tau[i]);
        }

        return Z;
        }

private:
    const FosterNetwork* plant_;
};