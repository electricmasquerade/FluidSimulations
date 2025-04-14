
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

float Simulation::calculateDensity(const std::shared_ptr<Particle> &particle, const std::vector<std::shared_ptr<Particle>> &neighbors) {
    float density = 0.0f;
    const float smoothingLength = particle->getSmoothingLength();
    for (const auto &neighbor : neighbors) {
        float distance = particle->getNeighborDistance(*neighbor);
        if (distance < smoothingLength) {
            const float weight = smoothingFunction(distance, smoothingLength);
            density += neighbor->getMass() * weight;
        }
    }
    return density;

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

void Simulation::initParticles(std::vector<std::shared_ptr<Particle>> &particles) {
    if (particles.empty()) return;
    this->particles = particles;
    const int grid_size = static_cast<int>(std::sqrt(particles.size()));
    const float spacing = domainSize / grid_size;
    //Uniformly distribute particles across domain based on spacing for initial position
    for (int i = 0; i < particles.size(); ++i) {
        const float x = (i % static_cast<int>(grid_size)) * spacing;
        const float y = (i / static_cast<int>(grid_size)) * spacing;
        particles[i]->setPosition(Vec3(x, y, 0.0f));
        particles[i]->setMass(1.0f); // Set mass to 1.0 for all particles
        particles[i]->setSmoothingLength(cellSize); // Set smoothing length based on spacing
        particles[i]->setDensity(restDensity); // Set rest density for all particles
    }


}

void Simulation::updateParticles(const float dt) {
    // Build the spatial map for the current particles
    buildSpatialMap(particles);

    //Update all particles based on their neighbors and the kernel function
    for (auto &particle : particles) {
        float density = 0.0f;
        particle->reset();
        //only get variables like smoothing length once if possible
        const float smoothingLength = particle->getSmoothingLength();
        // Find neighbors
        std::vector<std::shared_ptr<Particle>> neighbors = findNeighbors(*particle, particles, particle->getSmoothingLength());
        // Calculate forces based on neighbors
        for (const auto &neighbor : neighbors) {
            float distance = particle->getNeighborDistance(*neighbor);
            if (distance < smoothingLength) {
                const float weight = smoothingFunction(distance, smoothingLength);
                //Calculate and update the density of the particle
                particle->setDensity(calculateDensity(particle, neighbors));
                //Calculate and update pressure at particle
                float pressure = stiffness*( particle->getDensity() - restDensity);

                //Now calculate forces based on pressure and density, navier stokes
                //TODO: add viscosity and surface tension - how on earth
                Vec3 pressureForce = (pressure / (particle->getDensity() * particle->getDensity())) * weight * (particle->getPosition() - neighbor->getPosition());

                particle->addForce(pressureForce);


            }
        }
        // Update the particle
        particle->update(dt);
    }
}
