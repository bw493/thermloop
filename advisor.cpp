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