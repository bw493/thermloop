// plant.cpp - thermloop, C++ port, first module
//
// Foster thermal network. The step response is the closed-form integral of the
// impulse response, so nothing needs to be integrated numerically. For each
// first-order stage the pair
//
//     h_i(t) = (1/C_i) e^{-t/tau_i}      (impulse response)
//     s_i(t) = R_i (1 - e^{-t/tau_i})    (step response)
//
// satisfies s_i = integral of h_i from 0 to t, because tau_i / C_i = R_i. The
// network response is the superposition of these independent stages.
#include <vector>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <utility>

class FosterNetwork {
public:
    FosterNetwork(std::vector<double> R, std::vector<double> C, double T_amb = 25.0)
        : R_(std::move(R)), C_(std::move(C)), T_amb_(T_amb) {
        tau_.resize(R_.size());
        for (std::size_t i = 0; i < R_.size(); ++i)
            tau_[i] = R_[i] * C_[i];
    }
// resistance unit: degrees Celsius per watt (temperature divided by power)
// capacitance unit: joules per degree Celsius (energy divided by temperature)

    double impulse(double t) const {
        if (t < 0.0) return 0.0;
        double h = 0.0;
        for (std::size_t i = 0; i < R_.size(); ++i)
            h += (1.0 / C_[i]) * std::exp(-t / tau_[i]);
        return h;
    }

    double step(double t, double p = 1.0) const {         // p is the power scaling factor
        // stage is a physical component of the model, namely one R-C pair
        //  step is the temperature response of the system to a power input impulse
        if (t < 0.0) return T_amb_;
        double rise = 0.0;
        for (std::size_t i = 0; i < R_.size(); ++i)
            rise += R_[i] * (1.0 - std::exp(-t / tau_[i]));
        return T_amb_ + p * rise;
    }

    // Total DC thermal resistance. Because each s_i is the integral of h_i, the
    // area under the impulse response equals this sum exactly, with no
    // integration required.
    double total_resistance() const {
        double sum = 0.0;
        for (double r : R_) sum += r;
        return sum;
    }

    double ambient() const { return T_amb_; }

private:
    std::vector<double> R_, C_, tau_;
    double T_amb_;
};

FosterNetwork example_mcu_plant() {
    return FosterNetwork({2.0, 8.0, 20.0}, {0.25, 1.0, 6.0}, 25.0);
}

int main() {
    const FosterNetwork plant = example_mcu_plant();
    const double Rsum = plant.total_resistance();

    // The area under h(t) is Rsum by construction
    // s(∞) = ∫₀^∞ h(σ) dσ = Σ R_i.
    std::printf("sum of R (= area under h) = %.4f\n", Rsum);
    std::printf("steady step               = %.4f (expect %.4f)\n",
                plant.step(1.0e6, 2.0), plant.ambient() + 2.0 * Rsum);
    return 0;
}