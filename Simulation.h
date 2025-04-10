
#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include "Particle.h"


class Simulation {

public:
    float smoothingFunction(float r, float h);
    std::vector<Particle> findNeighbors(const Particle &particle, const std::vector<Particle> &particles, float radius);
    void updateParticles(float dt);

private:
    std::vector<Particle> particles;

};



#endif //SIMULATION_H
