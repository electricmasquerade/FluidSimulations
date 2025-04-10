
#include "Simulation.h"

float Simulation::smoothingFunction(const float r, const float h) { //r is the radius, h is the smoothing length
    // Implement the smoothing function here
    // using poly6 kernel for now
    if (r < h) {
        float weight = 315.0f / (64.0f * M_PI * std::pow(h, 9)) * std::pow(h * h - r * r, 3);
        return weight * pow((h*h - r*r), 3);
    }
    return 0.0f;
}

std::vector<Particle> Simulation::findNeighbors(const Particle &particle, const std::vector<Particle> &particles,
    const float radius) {
    std::vector<Particle> neighbors;
    for (const auto &other : particles) {
        if (&particle != &other && particle.getNeighborDistance(other) < radius) {
            neighbors.push_back(other);
        }
    }
    return neighbors;
}

void Simulation::updateParticles(float dt) {
}
