#include <cmath>
#include <cstdio>

class FosterStage { // a single R-C pair
    public:
    FosterStage(double R, double C, double T_amb = 25.0) {
        R_ = R;       // how hard it is for heat to escape
        C_ = C;       // heat needed per degree of temperature rise
        T_amb_ = T_amb;

        tau_ = R_ * C_;      // time constant
        invC_ = 1.0 / C_;    // temperature jump per unit of instantaneous heat
    }

    double resistance() const { return R_; }
    double capacitance() const { return C_; }
    double tau() const { return tau_; }
    double ambient() const { return T_amb_; }

    double impulse(double t) const {
        if (t < 0.0) return 0.0;
        return invC_ * std::exp(-t / tau_);
    }

    double step(double t, double p = 1.0) const {
        if (t < 0.0) return T_amb_;
        return T_amb_ + p * R_ * (1.0 - std::exp(-t / tau_));
    }

    private:
    double R_;
    double C_;
    double tau_;
    double invC_;
    double T_amb_;
};