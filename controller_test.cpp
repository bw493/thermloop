#include <iostream>
#include "controller.h"

using namespace std;

int main() {
    // Continuous PI reference: fine-grained integrator as "ground truth"
    double Kp = 2.0;
    double Ki = 0.5;
    double dt_fine = 0.0001;   // fine step for continuous approximation
    double dt_actual = 0.01;   // actual discrete controller step

    double e0 = 25.0;     // initial error magnitude (matches original 50 - 25)
    double tau = 0.001;    // decay time constant, seconds

    double t = 0.0;
    double t_end = 50.0;

    PIController discretePI(Kp, Ki, dt_actual);
    ContinuousPIReference continuousPI(Kp, Ki, dt_fine);

    double nextControlUpdate = 0.0;
    double u_discrete = 0.0;

    while (t <= t_end) {
        double error = e0 * exp(-t / tau);   // <-- exponentially decaying error

        double u_continuous = continuousPI.update(error);

        // Discrete Tustin update only fires at actual dt intervals
        if (t >= nextControlUpdate) {
            u_discrete = discretePI.update(error);
            nextControlUpdate += dt_actual;
        }

        cout << "t = " << t
             << "  u_continuous = " << u_continuous
             << "  u_discrete = " << u_discrete << endl;

        t += dt_fine;
    }

    return 0;
}