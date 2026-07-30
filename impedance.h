// impedance.h - thermloop, driving-point thermal impedance Z(s)
#pragma once
#include "plant.h"
#include <complex>
#include <vector>
#include <cstddef>
#include <cmath>

class ThermalImpedance {
public:
    explicit ThermalImpedance(const FosterNetwork& plant)
        : plant_(&plant) {}

    // Z(s) = sum_i R_i / (1 + s * tau_i)
    std::complex<double> evaluate(std::complex<double> s) const {
        std::complex<double> Z(0.0, 0.0);
        const std::vector<double>& R = plant_->resistances();
        const std::vector<double>& tau = plant_->taus();

        for (std::size_t i = 0; i < R.size(); ++i) {
            Z += R[i] / (1.0 + s * tau[i]);
        }
        return Z;
    }

    // Frequency response Z(j*omega), with omega in rad/s.
    std::complex<double> response(double omega) const {
        return evaluate(std::complex<double>(0.0, omega));
    }

    double magnitude(double omega) const {
        return std::abs(response(omega));
    }

    double magnitudeDB(double omega) const {
        return 20.0 * std::log10(magnitude(omega));
    }

    // Phase in degrees. Negative here, since every pole is real and in the LHP.
    double phaseDegrees(double omega) const {
        const double pi = 3.14159265358979323846;
        return std::arg(response(omega)) * 180.0 / pi;
    }

    // DC value: Z(0) should equal the sum of the R_i.
    double dcResistance() const {
        return evaluate(std::complex<double>(0.0, 0.0)).real();
    }

private:
    const FosterNetwork* plant_;
};