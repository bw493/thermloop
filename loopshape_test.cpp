// loopshape_test.cpp - thermloop, validation for fanplant.h + loopshape.h
//
// Three checks, in the usual style: analytic where a closed form exists,
// empirical round-trip where one does not.
//
//   1. FanPlant statics: R_conv(u), total resistance, and steady
//      temperature against hand-computed values.
//   2. suggestPI round-trip: ask for a crossover and a phase margin,
//      then hand the suggested gains back to crossoverFrequency() /
//      phaseMarginDegrees() and confirm they deliver what was ordered.
//   3. Closed loop in the time domain: the REAL PIController from
//      controller.h (clamped to [0,1], anti-windup active) drives the
//      NONLINEAR fan plant integrated stage-by-stage with Euler steps.
//      The linear analysis promised stability; the nonlinear simulation
//      has to keep that promise at the setpoint.
//
// Plant numbers are the identify.cpp self-test defaults, same as
// focus_throttle.ino, so every layer of the platform is talking about
// the same piece of hardware.
#include "plant.h"
#include "fanplant.h"
#include "loopshape.h"
#include "controller.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

int main() {
    // The identify.cpp self-test plant: 2 stages, last one is convection.
    FosterNetwork rig({4.1474, 8.8781}, {1.5880, 14.6102}, 24.0);
    double kFan = 3.0;                    // full fan quarters R_conv
    FanPlant fanned(rig, kFan);

    cout << "=== 1. FanPlant statics ===" << endl;
    double P0 = 1.0;                      // W, operating-point power
    double R_conv_still = fanned.convResistance(0.0);
    double R_conv_full  = fanned.convResistance(1.0);
    cout << "R_conv(u=0)   = " << R_conv_still
         << "  (expect 8.8781)" << endl;
    cout << "R_conv(u=1)   = " << R_conv_full
         << "  (expect " << 8.8781 / 4.0 << ")" << endl;
    cout << "T_ss(P=1,u=0) = " << fanned.steadyTemperature(P0, 0.0)
         << "  (expect " << 24.0 + 1.0 * (4.1474 + 8.8781) << ")" << endl;
    cout << "T_ss(P=1,u=1) = " << fanned.steadyTemperature(P0, 1.0)
         << "  (expect " << 24.0 + 1.0 * (4.1474 + 8.8781 / 4.0) << ")"
         << endl;

    cout << endl << "=== 2. suggestPI round-trip ===" << endl;
    // Linearize where the loop will actually sit: partial fan.
    double u0 = 0.3;
    LoopShaper shaper(fanned, P0, u0);

    double wantWc = 0.05;                 // rad/s, gentle for a slow rig
    double wantPM = 60.0;                 // degrees
    LoopShaper::PIGains g = shaper.suggestPI(wantWc, wantPM);
    cout << "suggested Kp  = " << g.Kp << endl;
    cout << "suggested Ki  = " << g.Ki
         << "  (feasible iff > 0)" << endl;

    double gotWc = shaper.crossoverFrequency(g.Kp, g.Ki);
    double gotPM = shaper.phaseMarginDegrees(g.Kp, g.Ki);
    cout << "crossover     = " << gotWc << "  (asked " << wantWc << ")"
         << endl;
    cout << "phase margin  = " << gotPM << "  (asked " << wantPM << ")"
         << endl;

    cout << endl << "=== 3. nonlinear closed loop ===" << endl;
    // Setpoint must live between the fan-off and fan-full steady states.
    double T_target = 33.0;
    double dt = 0.5;                      // same tick as focus_throttle.ino
    PIController fan(g.Kp, g.Ki, dt, 0.0, 1.0);

    // Per-stage Euler integration of the nonlinear plant. Stage i holds
    // x[i], its rise above ambient:  C_i dx_i/dt = P - x_i / R_i(u),
    // with only the last R_i depending on the fan command.
    vector<double> R = {4.1474, 8.8781};
    vector<double> C = {1.5880, 14.6102};
    vector<double> x = {0.0, 0.0};        // cold start at ambient
    double u = 0.0;

    double T = 24.0;
    double h = 0.05;                      // s, fine ODE step under dt
    int ticks = (int)(3000.0 / dt);       // 3000 s of simulated time
    for (int k = 0; k < ticks; ++k) {
        // controller runs once per tick, on the measured temperature
        u = fan.update(T - T_target);
        // plant integrates at the fine step between ticks
        for (double s = 0.0; s < dt; s += h) {
            double Rc = R[1] / (1.0 + kFan * u);
            x[0] += h * (P0 - x[0] / R[0]) / C[0];
            x[1] += h * (P0 - x[1] / Rc)  / C[1];
        }
        T = 24.0 + x[0] + x[1];
    }
    cout << "T after 3000 s = " << T << "  (target " << T_target << ")"
         << endl;
    cout << "fan command    = " << u << "  (must sit inside [0,1])"
         << endl;

    return 0;
}
