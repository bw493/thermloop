#pragma once

#include "plant.h"
#include <vector>
#include <cmath>
#include <cstddef>

using namespace std;

class predictor {
public:
    predictor(const FosterNetwork& plant, double dt)
        : plant_(plant), dt_(dt), t_(0.0) { // dt is sampling interval

        // track data type of resistance and tau
        const vector<double>& R = plant_.resistances();
        const vector<double>& tau = plant_.taus();

        // track time constant maximum to determine how many time steps to simulate
        double tau_max = 0.0;
        for (size_t i = 0; i < tau.size(); ++i) {
            if (tau[i] > tau_max) tau_max = tau[i];
        }

        // decide the window of computation:
        // beyond 5 time constants the step response has effectively stopped decaying
        double raw_N = 5.0 * tau_max / dt_;
        double N_rounded_up = ceil(raw_N);
        size_t N = (size_t)N_rounded_up;
        h_.resize(N); // kernel size is N, the number of time steps
        // the kernel is a fixed-size array of weights describing
        // how the system responds over time to a single brief input

        for (size_t m = 0; m < N; ++m) {
            double val = 0.0;
            for (size_t i = 0; i < R.size(); ++i) {
                val += R[i] * (exp(-(double)m * dt_ / tau[i])
                             - exp(-(double)(m + 1) * dt_ / tau[i]));
            }
            h_[m] = val;
        }
    }

    double predict(double p) {
        double T = plant_.ambient();
        for (double h_i : h_) {
            T += p * h_i;
        }
        return T;
    }

private:
    FosterNetwork plant_;
    double dt_;
    double t_;
    vector<double> h_;
};
