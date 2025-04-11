
#include "Simulation.h"

float Simulation::smoothingFunction(const float r, const float h) { //r is the radius, h is the smoothing length
    // Implement the smoothing function here
    // using poly6 kernel for now
    if (r < h) {
        const float weight = 315.0f / (64.0f * M_PI * std::pow(h, 9)) * std::pow(h * h - r * r, 3);
        return weight * pow((h*h - r*r), 3);
    }
    return 0.0f;
}

void Simulation::buildSpatialMap(const std::vector<Particle> &particles) {
    for (size_t i = 0; i < particles.size(); i++) {
        Vec3 pos = particles[i].getPosition();
        int cellX = static_cast<int>(pos[0] / cellSize);
        int cellY = static_cast<int>(pos[1] / cellSize);
        CellKey key{cellX, cellY};
        spatialMap[key].push_back(particles[i]);
    }

}

std::vector<Particle> Simulation::findNeighbors(const Particle &particle, const std::vector<Particle> &particles,
                                                const float radius) {
    std::vector<Particle> neighbors;
    //TODO: use spatial map to find neighbors
    return neighbors;
}

void Simulation::updateParticles(float dt) {
}
