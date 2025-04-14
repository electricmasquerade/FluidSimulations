
#include "Simulation.h"

float Simulation::smoothingFunction(const float r, const float h) { //r is the radius, h is the smoothing length
    // Implement the smoothing function here
    // using poly6 kernel for now, only defined where r < h
    if (r < h) {
        const float weight = 315.0f / (64.0f * M_PI * std::pow(h, 9)) * std::pow(h * h - r * r, 3);
        return weight * 2000.0f; // Scale the weight to a reasonable value
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
    density = std::max(density, 0.01f);  // or even 1e-3f

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
    //std::cout << "Neighbors: " << neighbors.size() << std::endl;

    return neighbors;
}

void Simulation::initParticles(std::vector<std::shared_ptr<Particle>> &particles) {
    if (particles.empty()) return;
    this->particles = particles;
    const int grid_size = static_cast<int>(std::sqrt(particles.size()));
    const float spacing = domainSize / static_cast<float>(grid_size);
    //Uniformly distribute particles across domain based on spacing for initial position
    for (int i = 0; i < particles.size(); ++i) {
        const float x = (i % static_cast<int>(grid_size)) * spacing;
        const float y = (i / static_cast<int>(grid_size)) * spacing;
        particles[i]->setPosition(Vec3(x, y, 0.0f));
        particles[i]->setMass(1.0f); // Set mass to 1.0 for all particles
        particles[i]->setSmoothingLength(spacing * 1.5); // Set smoothing length based on spacing
        particles[i]->setDensity(restDensity); // Set rest density for all particles
    }


}

void Simulation::updateParticles(const float dt) {
    float viscosity = 0.1f;
    // Build the spatial map for the current particles
    buildSpatialMap(particles);

    //Update all particles with their density before doing anything else.
    for (auto &particle : particles) {
        particle->reset();
        // Find neighbors
        std::vector<std::shared_ptr<Particle>> neighbors = findNeighbors(*particle, particles, particle->getSmoothingLength());
        // Calculate density for the particle based on neighbors
        particle->setDensity(calculateDensity(particle, neighbors));
    }
    //Now update pressure for all particles.
    for (auto &particle : particles) {
        float pressure = stiffness * (particle->getDensity() - restDensity);
        pressure = std::clamp(pressure, 0.0f, 100.0f);  // You can tune this max value
        particle->setPressure(pressure);

    }
    //Now update forces for all particles.
    for (auto &particle : particles) {
        // Find neighbors
        std::vector<std::shared_ptr<Particle>> neighbors = findNeighbors(*particle, particles, particle->getSmoothingLength());
        // Calculate forces based on neighbors
        for (const auto &neighbor : neighbors) {
            if (particle == neighbor) continue; // Skip self
            float distance = particle->getNeighborDistance(*neighbor);
            if (distance < particle->getSmoothingLength()) {
                // Calculate force based on pressure and density differences
                float forceMagnitude = stiffness * (particle->getPressure() / (particle->getDensity() * particle->getDensity()) +
                                                    neighbor->getPressure() / (neighbor->getDensity() * neighbor->getDensity()));
                Vec3 forceDirection = unitVector(particle->getPosition() - neighbor->getPosition());
                Vec3 viscosityForce = viscosity * (neighbor->getVelocity() - particle->getVelocity());

                Vec3 force = forceMagnitude * forceDirection + viscosityForce;
                //apply gravity
                force += Vec3(0, 0.05f*gravity * particle->getDensity(), 0);
                const Vec3 drag = -0.1f * particle->getVelocity();
                force += drag;

                float maxForce = 1000.0f;
                if (forceMagnitude > maxForce) {
                    force = maxForce * forceDirection;
                }
                //std::cout << "forceMagnitude: " << forceMagnitude << std::endl;
                particle->addForce(force);
            }
        }
    }
    //Now update all particles with their forces
    for (auto &particle : particles) {
        //std::cout << "ρ: " << particle->getDensity() << "  P: " << particle->getPressure() << std::endl;
        particle->update(dt);
    }
}
