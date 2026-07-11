#include "plant.h"
#include "controller.h"

class CoffeeAdvisor {
public:
    CoffeeAdvisor(const FosterNetwork& cup, double T_serve, double T_drink, const PIController& blow)
        : cup_(cup), T_serve_(T_serve), T_drink_(T_drink), blow_(blow) {}

private:
    FosterNetwork cup_;
    double T_serve_;
    double T_drink_;
    PIController blow_;
};

double timeToDrinkable() {
    // Implementation for calculating time to drinkable temperature
    T = T_serve_;
    while (T > T_drink_) {
        // Calculate the temperature change using the cup's step response
        T = cup_.step(time, power_input); // power_input can be determined based on the controller's output
        time += dt; // Increment time by the time step
    }
}