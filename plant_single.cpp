#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

class FosterStage { // a single R-C pair
    public:
    FosterStage(double R, double C, double T_amb = 25.0) {
        R_ = R;       // how hard it is for heat to escape
        C_ = C;       // heat needed per degree of temperature rise
        T_amb_ = T_amb; // room temp

        tau_ = R_ * C_;      // time constant
        invC_ = 1.0 / C_;    // temperature jump per unit of instantaneous heat
    }

    double resistance() const { return R_; }
    double capacitance() const { return C_; }
    double tau() const { return tau_; }
    double ambient() const { return T_amb_; }

    double impulse(double t) const { // the power input at time t = 0
        if (t < 0.0) return 0.0;
        return invC_ * exp(-t / tau_);
    }

    double step(double t, double p = 1.0) const { // temperature response to a sustained, constant power input turned on at t=0
        if (t < 0.0) return T_amb_;
        return T_amb_ + p * R_ * (1.0 - exp(-t / tau_));
    }

    private:
    double R_;
    double C_;
    double tau_;
    double invC_;
    double T_amb_;
};

int main() {
    const FosterStage stage(10.0, 2.0, 25.0);
    const double p = 3.0;

    cout << fixed << setprecision(4);

    // Preliminary Test 1: right after turning on the power, has the
    // temperature barely changed from room temperature?
    double temp_at_start = stage.step(0.0, p);
    cout << "Preliminary Test 1: Temperature right when power turns on" << endl;
    cout << "Actual temperature:   " << temp_at_start << endl;
    cout << "Room temperature:     " << stage.ambient() << endl;
    cout << "These should match, since no time has passed yet." << endl;
    cout << endl;

    // Preliminary Test 2: does the impulse response start at its
    // highest point, exactly 1/C, at t = 0?
    double impulse_at_start = stage.impulse(0.0);
    double expected_impulse_start = 1.0 / stage.capacitance();
    cout << "Preliminary Test 2: Impulse response at the instant of the heat burst" << endl;
    cout << "Actual value:      " << impulse_at_start << endl;
    cout << "Expected value:    " << expected_impulse_start << endl;
    cout << "These should match, since the burst deposits heat all at once." << endl;
    cout << endl;

    // Preliminary Test 3: does the impulse response shrink as time
    // goes on, since the heat leaks away and is not replenished?
    double impulse_early = stage.impulse(1.0);
    double impulse_later = stage.impulse(10.0);
    cout << "Preliminary Test 3: Impulse response decays over time" << endl;
    cout << "Value at t = 1:    " << impulse_early << endl;
    cout << "Value at t = 10:   " << impulse_later << endl;
    cout << "The second number should be smaller, since the heat is leaking away." << endl;
    cout << endl;

    // Preliminary Test 4: does the step response only ever rise, never
    // overshoot or dip, as time goes on?
    double temp_early = stage.step(1.0, p);
    double temp_later = stage.step(10.0, p);
    cout << "Preliminary Test 4: Temperature keeps rising over time" << endl;
    cout << "Temperature at t = 1:    " << temp_early << endl;
    cout << "Temperature at t = 10:   " << temp_later << endl;
    cout << "The second number should be larger, since heat keeps building up." << endl;
    cout << endl;

    // Final Test 1: the area under the impulse response over all time
    // should equal the resistance R, since integrating (1/C) * exp(-t/RC)
    // from 0 to infinity yields exactly R.
    const double dt = 0.01;
    double area = 0.0;
    double prev = stage.impulse(0.0);
    for (double t = dt; t <= 2000.0; t += dt) {
        double v = stage.impulse(t);
        area += 0.5 * (v + prev) * dt;
        prev = v;
    }
    cout << "Final Test 1: Does the impulse response add up to R?" << endl;
    cout << "Resistance R is set to:      " << stage.resistance() << endl;
    cout << "Sum of impulse response:     " << area << endl;
    cout << "These two numbers should match." << endl;
    cout << endl;

    // Final Test 2: as t grows very large, the step response should
    // converge to T_amb + p * R.
    double steady = stage.step(1.0e6, p);
    double expected_steady = stage.ambient() + p * stage.resistance();
    cout << "Final Test 2: Does the temperature stop changing eventually?" << endl;
    cout << "Temperature after a long time:   " << steady << endl;
    cout << "Expected final temperature:      " << expected_steady << endl;
    cout << "These two numbers should match." << endl;
    cout << endl;

    // Final Test 3: at exactly one time constant (t = tau), the step
    // response should have risen to (1 - 1/e) of its total rise, a
    // standard property of first-order RC systems.
    double at_tau = stage.step(stage.tau(), p);
    double expected_at_tau = stage.ambient() + p * stage.resistance() * (1.0 - exp(-1.0));
    cout << "Final Test 3: Is the temperature correct at one time constant?" << endl;
    cout << "Temperature at t = tau:          " << at_tau << endl;
    cout << "Expected temperature at t = tau:  " << expected_at_tau << endl;
    cout << "These two numbers should match." << endl;

    return 0;
}