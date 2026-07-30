// identify.cpp - thermloop, system identification from a step response
//
// Removes the last guesswork from the platform: PLANT_R and PLANT_C in
// focus_throttle.ino must come from the rig, not from a plausible-looking
// constant. Run the sketch's CAL mode, capture the CSV, feed it here.
//
// Method: separable (variable projection) least squares.
//
//   The model  y(t) = sum_i R_i * (1 - exp(-t / tau_i))
//   is NONLINEAR in tau but LINEAR in R. So for any fixed pair of taus the
//   best R is a closed-form 2x2 solve, and only the two taus need to be
//   searched. That collapses a four-parameter nonlinear fit into a
//   two-parameter grid search with an exact inner solve -- far more robust
//   than throwing all four at a gradient method, and simple enough to read.
//
// Build: g++ -std=c++17 -O2 identify.cpp -o identify
// Usage: ./identify step.csv <P_watts> [T_amb]
//        ./identify --selftest
//
// CSV format: two columns, "t_seconds,temperature_C", one sample per line.
// A leading header line is skipped. T_amb defaults to the first sample.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct Fit {
    double R1, R2, tau1, tau2, rms;
};

// ---------------------------------------------------------------------
// Inner solve: taus fixed, R by ordinary least squares on the 2x2 normal
// equations. Returns false if the two basis functions are collinear, which
// happens when the taus are pushed too close together.
// ---------------------------------------------------------------------
bool solveR(const vector<double>& t, const vector<double>& y,
            double tau1, double tau2, double& R1, double& R2, double& sse) {
    double S11 = 0.0, S12 = 0.0, S22 = 0.0, b1 = 0.0, b2 = 0.0;
    for (size_t j = 0; j < t.size(); ++j) {
        double p1 = 1.0 - exp(-t[j] / tau1);
        double p2 = 1.0 - exp(-t[j] / tau2);
        S11 += p1 * p1;
        S12 += p1 * p2;
        S22 += p2 * p2;
        b1  += p1 * y[j];
        b2  += p2 * y[j];
    }
    double det = S11 * S22 - S12 * S12;
    if (fabs(det) < 1e-12) return false;

    R1 = (b1 * S22 - b2 * S12) / det;
    R2 = (S11 * b2 - S12 * b1) / det;

    // A negative resistance is not physical; reject rather than report it.
    if (R1 < 0.0 || R2 < 0.0) return false;

    sse = 0.0;
    for (size_t j = 0; j < t.size(); ++j) {
        double m = R1 * (1.0 - exp(-t[j] / tau1))
                 + R2 * (1.0 - exp(-t[j] / tau2));
        double d = y[j] - m;
        sse += d * d;
    }
    return true;
}

// ---------------------------------------------------------------------
// Outer search: log-spaced grid over the two taus, then two refinement
// passes around the winner. Log spacing matters -- thermal time constants
// span decades, so a linear grid wastes nearly all its points.
// ---------------------------------------------------------------------
Fit identify(const vector<double>& t, const vector<double>& y,
             double tauLo, double tauHi) {
    Fit best;
    best.R1 = best.R2 = best.tau1 = best.tau2 = 0.0;
    best.rms = 1e300;

    for (int pass = 0; pass < 3; ++pass) {
        int N = 80;
        double lo = log(tauLo), hi = log(tauHi);
        double bestSse = 1e300;
        double b1t = tauLo, b2t = tauHi, bR1 = 0.0, bR2 = 0.0;

        for (int i = 0; i < N; ++i) {
            double ta = exp(lo + (hi - lo) * i / (N - 1));
            for (int k = i + 1; k < N; ++k) {          // enforce tau1 < tau2
                double tb = exp(lo + (hi - lo) * k / (N - 1));
                if (tb < ta * 1.5) continue;           // keep them separable
                double R1, R2, sse;
                if (!solveR(t, y, ta, tb, R1, R2, sse)) continue;
                if (sse < bestSse) {
                    bestSse = sse;
                    b1t = ta; b2t = tb; bR1 = R1; bR2 = R2;
                }
            }
        }

        best.tau1 = b1t; best.tau2 = b2t;
        best.R1 = bR1;   best.R2 = bR2;
        best.rms = sqrt(bestSse / (double)t.size());

        // Narrow the window around the winner and search again.
        tauLo = b1t / 3.0;
        tauHi = b2t * 3.0;
        if (tauLo < 1e-3) tauLo = 1e-3;
    }
    return best;
}

void report(const Fit& f, double P) {
    double C1 = f.tau1 / f.R1;
    double C2 = f.tau2 / f.R2;
    cout << fixed << setprecision(4);
    cout << "  step power     = " << P << " W" << endl;
    cout << "  R              = {" << f.R1 << ", " << f.R2 << "} C/W" << endl;
    cout << "  C              = {" << C1 << ", " << C2 << "} J/C" << endl;
    cout << "  tau            = {" << f.tau1 << ", " << f.tau2 << "} s" << endl;
    cout << "  total R        = " << f.R1 + f.R2 << " C/W" << endl;
    cout << "  fit RMS        = " << f.rms << " C/W" << endl;
    cout << endl;
    cout << "  Paste into focus_throttle.ino:" << endl;
    cout << "    const float PLANT_R[2] = {" << f.R1 << ", " << f.R2 << "};"
         << endl;
    cout << "    const float PLANT_C[2] = {" << C1 << ", " << C2 << "};"
         << endl;
}

// ---------------------------------------------------------------------
// NTC beta from two calibration points. Beta is defined by
//   R(T) = R_ref * exp(beta * (1/T - 1/T_ref)),  T in kelvin
// so two measurements determine it directly. Use ice water and a warm
// bath if you have no datasheet, and prefer points far apart.
// ---------------------------------------------------------------------
double ntcBeta(double R1, double T1_C, double R2, double T2_C) {
    double T1 = T1_C + 273.15;
    double T2 = T2_C + 273.15;
    return log(R1 / R2) / (1.0 / T1 - 1.0 / T2);
}

// ---------------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc >= 2 && string(argv[1]) == "--selftest") {
        // Generate a step response from KNOWN parameters, add realistic
        // 10 bit quantisation noise, and check that the fit recovers them.
        double R1 = 4.0, R2 = 9.0, C1 = 1.5, C2 = 14.0;
        double tau1 = R1 * C1, tau2 = R2 * C2;
        double P = 1.5, T_amb = 24.0, dt = 0.5;

        vector<double> t, y;
        unsigned seed = 12345;
        for (int k = 0; k < 1600; ++k) {           // 800 s of data
            double tt = k * dt;
            double rise = P * (R1 * (1.0 - exp(-tt / tau1))
                             + R2 * (1.0 - exp(-tt / tau2)));
            seed = seed * 1103515245u + 12345u;    // deterministic noise
            double noise = ((double)((seed >> 16) & 0x7fff) / 32767.0 - 0.5)
                           * 0.25;                 // +/- 0.125 C
            t.push_back(tt);
            y.push_back((T_amb + rise + noise - T_amb) / P);
        }

        cout << "=== self-test: recover known parameters ===" << endl;
        cout << fixed << setprecision(4);
        cout << "  truth R        = {" << R1 << ", " << R2 << "} C/W" << endl;
        cout << "  truth C        = {" << C1 << ", " << C2 << "} J/C" << endl;
        cout << "  truth tau      = {" << tau1 << ", " << tau2 << "} s"
             << endl << endl;

        Fit f = identify(t, y, 0.5, 400.0);
        report(f, P);

        cout << endl << "=== NTC beta from two points ===" << endl;
        // 10k at 25 C, 3.6k at 50 C -- typical of a 3950 part.
        double beta = ntcBeta(10000.0, 25.0, 3603.0, 50.0);
        cout << "  beta           = " << beta << " K" << endl;
        return 0;
    }

    if (argc < 3) {
        cout << "usage: identify <step.csv> <P_watts> [T_amb_C]" << endl;
        cout << "       identify --selftest" << endl;
        return 1;
    }

    string path = argv[1];
    double P = atof(argv[2]);
    ifstream in(path.c_str());
    if (!in) {
        cout << "cannot open " << path << endl;
        return 1;
    }

    vector<double> t, T;
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string a, b;
        if (!getline(ss, a, ',')) continue;
        if (!getline(ss, b, ',')) continue;
        char* endA = 0;
        double ta = strtod(a.c_str(), &endA);
        if (endA == a.c_str()) continue;            // header line, skip
        t.push_back(ta);
        T.push_back(atof(b.c_str()));
    }
    if (t.size() < 20) {
        cout << "need at least 20 samples, got " << t.size() << endl;
        return 1;
    }

    double T_amb = (argc >= 4) ? atof(argv[3]) : T[0];
    double t0 = t[0];
    vector<double> y;
    for (size_t j = 0; j < t.size(); ++j) {
        t[j] -= t0;
        y.push_back((T[j] - T_amb) / P);
    }

    cout << "=== identified from " << path << " ("
         << t.size() << " samples) ===" << endl;
    cout << fixed << setprecision(4) << "  T_amb          = " << T_amb
         << " C" << endl;
    Fit f = identify(t, y, 0.5, t[t.size() - 1]);
    report(f, P);
    return 0;
}
