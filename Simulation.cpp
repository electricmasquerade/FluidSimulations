
#include "Simulation.h"

float Simulation::smoothingFunction(const float r, const float h) { //r is the radius, h is the smoothing length
    // Implement the smoothing function here
    // using poly6 kernel for now, only defined where r < h
    if (r < h) {
        const float weight = 315.0f / (64.0f * M_PI * std::pow(h, 9)) * std::pow(h * h - r * r, 3);
        return weight;
    }
    return 0.0f;
}

void Simulation::buildSpatialMap(const std::vector<std::shared_ptr<Particle>> &particles) {
    spatialMap.clear();
    for (const auto & particle : particles) {
        Vec3 pos = particle->getPosition();
        const int cellX = static_cast<int>(std::floor(pos[0] / cellSize));
        const int cellY = static_cast<int>(std::floor(pos[1] / cellSize));
        CellKey key{cellX, cellY};
        spatialMap[key].push_back(particle);
    }

}

std::vector<std::shared_ptr<Particle>> Simulation::findNeighbors(const Particle &particle, const std::vector<std::shared_ptr<Particle>> &particles,
                                                const float radius) {
    std::vector<std::shared_ptr<Particle>> neighbors;
    Vec3 pos = particle.getPosition();
    int cellX = static_cast<int>(pos[0] / cellSize);
    int cellY = static_cast<int>(pos[1] / cellSize);
    // Check home cell and all neighboring cells, only one row deep for now
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            CellKey key{cellX + dx, cellY + dy};
            if (spatialMap.find(key) != spatialMap.end()) {
                for (const auto & neighbor : spatialMap[key]) {
                    if (&particle == neighbor.get()) {
                        continue; // Skip self
                    }
                    float distance = particle.getNeighborDistance(*neighbor);
                    if (distance < radius) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }
    }
    return neighbors;
}

void Simulation::updateParticles(const float dt) {
    // Build the spatial map for the current particles
    buildSpatialMap(particles);

    //Update all particles based on their neighbors and the kernel function
    for (auto &particle : particles) {
        particle->reset();
        // Find neighbors
        std::vector<std::shared_ptr<Particle>> neighbors = findNeighbors(*particle, particles, particle->getSmoothingLength());
        // Calculate forces based on neighbors
        for (const auto &neighbor : neighbors) {
            float distance = particle->getNeighborDistance(*neighbor);
            if (distance < particle->getSmoothingLength()) {
                float weight = smoothingFunction(distance, particle->getSmoothingLength());
                //Calculate the density of the particle
                particle->setDensity(particle->getDensity() + neighbor->getDensity() * weight);
                // Calculate the force based on the weight and the distance
                Vec3 force = weight * (neighbor->getPosition() - particle->getPosition());
                particle->addForce(force);
            }
        }
        // Update the particle
        particle->update(dt);
    }
}
