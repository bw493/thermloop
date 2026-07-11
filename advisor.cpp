// advisor.cpp - thermloop, coffee-advisor composition layer
// Stitches the mug (plant.h) and the stove dial (controller.h) into one
// customer-facing object. timeToDrinkable() answers the bar-side question:
// "how long until my coffee is ready to sip?"
#include "plant.h"
#include "controller.h"
#include <cmath>
#include <iostream>

using namespace std;

class CoffeeAdvisor {
public:
    CoffeeAdvisor(const FosterNetwork& cup, double T_serve, double T_drink,
                  const PIController& stove, double dt)
        : cup_(cup), T_serve_(T_serve), T_drink_(T_drink),
          stove_(stove), dt_(dt) {}

    // How long the cup takes to fall from T_serve to T_drink on its own.
    // The descent is pure natural cooling: the error stays negative, so the
    // stove is clamped off the whole way and contributes nothing. That makes
    // this a closed-form exponential decay - no loop, no controller.
    //   t = tau * ln( (T_serve - T_amb) / (T_drink - T_amb) )
    double timeToDrinkable() const {
        double tau   = cup_.taus()[0];
        double T_amb = cup_.ambient();
        return tau * log((T_serve_ - T_amb) / (T_drink_ - T_amb));
    }

private:
    FosterNetwork cup_;
    double T_serve_, T_drink_;
    PIController stove_;   // held for the hold phase; idle during the descent
    double dt_;
};

int main() {
    // A single-stage mug cooling toward a 25 C room.
    // R = 30 C/W, C = 10 J/C  ->  tau = 300 s.
    FosterNetwork mug({30.0}, {10.0}, 25.0);

    double dt = 1.0;                     // one-second control tick
    PIController stove(0.8, 0.05, dt);   // the dial the advisor holds with

    double T_serve = 85.0;               // fresh off the roaster
    double T_drink = 60.0;               // the customer's ideal first sip

    CoffeeAdvisor advisor(mug, T_serve, T_drink, stove, dt);

    double wait = advisor.timeToDrinkable();

    cout << "Fresh cup on the bar at " << T_serve << " C." << endl;
    cout << "You'd like it at " << T_drink << " C." << endl;
    cout << "Give it about " << wait << " s ("
         << wait / 60.0 << " min) and it's ready to sip." << endl;

    return 0;
}