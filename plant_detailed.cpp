#include <vector>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <functional>

class FosterNetwork {
    public:
    FosterNetwork(std::vector<double> R, std::vector<double> C, double T_amb = 25.0) {
        R_ = R;
        C_ = C;
        T_amb_ = T_amb;
    
        tau_.resize(R_.size());

        for (std::size_t i=0; i<R_.size(); i++) {
            tau_[i]= R_[i] * C_[i];
        }
    }

    const std::vector<double>& taus() const { return tau_; }
    const std::vector<double>& impulseWeights() const { return invC_; }
    const std::vector<double>& resistances() const { return R_; }
    const std::vector<double>& capacitances() const { return C_; }
    double ambient() const { return T_amb_; }

    double totalResistance() const {
        double sum = 0.0;
        for (std::size_t i = 0; i < R_.size(); i++) {
            sum += R_[i];
        }
        return sum;
    }

    // Phase three: per stage responses
    double stageImpulse(std::size_t i, double t) const {
        if (t < 0.0) return 0.0;
        return invC_[i] * std::exp(-t / tau_[i]);
    }

    double stageStep(std::size_t i, double t) const {
        if (t < 0.0) return 0.0;
        return R_[i] * (1.0 - std::exp(-t / tau_[i]));
    }

    // Phase three: summed responses
    double impulse(double t) const {
        if (t < 0.0) return 0.0;
        double h = 0.0;
        for (std::size_t i = 0; i < R_.size(); i++) {
            h += stageImpulse(i, t);
        }
        return h;
    }

    double step(double t, double p = 1.0) const {
        if (t < 0.0) return T_amb_;
        double rise = 0.0;
        for (std::size_t i = 0; i < R_.size(); i++) {
            rise += p * stageStep(i, t);
        }
        return T_amb_ + rise;
    }

    private:
    std::vector<double> R_;
    std::vector<double> C_;
    std::vector<double> tau_;
    std::vector<double> invC_;
    double T_amb_;
};