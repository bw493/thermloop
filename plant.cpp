// plant.cpp - thermloop, C++ port, first module
#include <vector>
#include <cmath>
#include <cstddef>
#include <utility>
#include <cstdio>

class FosterNetwork {
public:
    FosterNetwork(std::vector<double> R, std::vector<double> C, double T_amb = 25.0) {
        R_ = R;
        C_ = C;
        T_amb_ = T_amb;
        tau_.resize(R_.size());
        for (std::size_t i = 0; i < R_.size(); ++i)
            tau_[i] = R_[i] * C_[i];
    }

    // Time constants tau_i = R_i * C_i, precomputed once at construction.
    const std::vector<double>& taus() const { return tau_; }

    double impulse(double t) const {
        if (t < 0.0) return 0.0;
        double h = 0.0;
        for (std::size_t i = 0; i < R_.size(); ++i)
            h += (R_[i] / tau_[i]) * std::exp(-t / tau_[i]);
        return h;
    }

    double step(double t, double p = 1.0) const { // p is the power scaling factor
        if (t < 0.0) return T_amb_;
        double rise = 0.0;
        for (std::size_t i = 0; i < R_.size(); ++i)
            rise += R_[i] * (1.0 - std::exp(-t / tau_[i]));
        return T_amb_ + p * rise;
    }

    const std::vector<double>& resistances() const { return R_; }
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

    double Rsum = 0.0;
    for (double r : plant.resistances())
        Rsum += r;

    // The area under h(t) equals the sum of the thermal resistances.
    const double dt = 0.05;
    double area = 0.0;
    double prev = plant.impulse(0.0);
    for (double t = dt; t <= 2000.0; t += dt) {
        double v = plant.impulse(t);
        area += 0.5 * (v + prev) * dt;
        prev = v;
    }

    std::printf("sum of R        = %.4f\n", Rsum);
    std::printf("area under h(t) = %.4f\n", area);
    std::printf("steady step     = %.4f (expect %.4f)\n",
                plant.step(1.0e6, 2.0), plant.ambient() + 2.0 * Rsum);
    return 0;
}