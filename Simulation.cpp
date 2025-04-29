
#include "Simulation.h"

float Simulation::smoothingPoly6(const float r, const float h) { //r is the radius, h is the smoothing length
    // Implement the smoothing function here
    // using 2D poly6 kernel for now, only defined where r < h
    if (r < h) {
        const float weight = 4.0f / (M_PI * std::pow(h, 8)) * std::pow(h * h - r * r, 3);
        return weight; //* 2000.0f; // Scale the weight to a reasonable value
    }
    return 0.0f;
}

Vec3 Simulation::spikyGradient(const Vec3 &r_vec, float h) {
    float r = r_vec.length();
    if (r > 0 && r < h) {
        const float coeff = -30.0f / (M_PI * std::pow(h,5));

        const float magnitude = coeff * (h-r) * (h-r);

        return -magnitude * (r_vec/r);
    }
    return {0, 0, 0};
}

Vec3 Simulation::viscosityLaplacian(const Vec3 &r_vec, float h) {
    //laplacian for viscosity in 2D
    float r = r_vec.length();
    if (r > 0 && r < h) {
        const float coeff = 40.0f / (M_PI * std::pow(h, 5));
        const float magnitude = coeff * (h - r);
        return magnitude * (r_vec / r);
    }
    return {0, 0, 0};
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
    for (const auto & particle : ghostParticles) {
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
            const float weight = smoothingPoly6(distance, smoothingLength);
            density += neighbor->getMass() * weight;
            density += particle->getMass() * smoothingPoly6(0, smoothingLength); // Add the mass of the current particle
        }
    }
    density = std::max(density, 0.01f);  // or even 1e-3f

    return density;

}

std::vector<std::shared_ptr<Particle>> Simulation::findNeighbors(const Particle &particle, const std::vector<std::shared_ptr<Particle>> &particles,
                                                                 const float radius) {
    std::vector<std::shared_ptr<Particle>> neighbors;
    Vec3 pos = particle.getPosition();
    const int cellX = static_cast<int>(pos[0] / cellSize);
    const int cellY = static_cast<int>(pos[1] / cellSize);
    constexpr int depth = 1;
    // Check home cell and all neighboring cells, only one row deep for now
    for (int dx = -depth; dx <= depth; ++dx) {
        for (int dy = -depth; dy <= depth; ++dy) {
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
    neighborList.clear();
    neighborList.resize(particles.size());
    ghostParticles.clear();
    if (particles.empty()) return;
    this->particles = particles;
    restDensity = particles[0]->getMass();
    const int grid_size = static_cast<int>(std::sqrt(particles.size()));
    const float spacing = domainSize / static_cast<float>(grid_size);
    //Uniformly distribute particles across domain based on spacing for the initial position
    for (int i = 0; i < particles.size(); ++i) {
        const float x = (i % static_cast<int>(grid_size)) * spacing;
        const float y = (i / static_cast<int>(grid_size)) * spacing;
        particles[i]->setPosition(Vec3(x, y, 0.0f));
        particles[i]->setMass(1.0f); // Set mass to 1.0 for all particles
        particles[i]->setSmoothingLength(cellSize); // Set a smoothing length based on spacing
        auto neighbors = findNeighbors(*particles[i], particles, particles[i]->getSmoothingLength());
        particles[i]->setDensity(calculateDensity(particles[i], neighbors));
        //restDensity = calculateDensity(particles[i], neighbors);
        float pressure = calculatePressure(particles[i]);

        particles[i]->setPressure(pressure);

    }

    //const float h = spacing * 1.5f;
    const float ghostSpacing = spacing * 0.25f;
    const int numGhostLayers = ceil(cellSize/ghostSpacing);
    const float offset = ghostSpacing;
    for (int i = 0; i < grid_size; ++i) {
        const float x = i * spacing;
        for (int layer = 1; layer <= numGhostLayers; ++layer) {
            constexpr float ghostMass = 1.0f;
            //bottom
            auto ghostBottom = std::make_shared<Particle>();
            ghostBottom->setPosition(Vec3(x, -layer * offset, 0.0f));
            ghostBottom->setMass(ghostMass); // Set mass to 1.0 for all particles
            ghostBottom->setSmoothingLength(cellSize); // Set smoothing length based on spacing
            auto ghostNeighbors = findNeighbors(*ghostBottom, particles, ghostBottom->getSmoothingLength());
            ghostBottom->setDensity(calculateDensity(ghostBottom, ghostNeighbors)); // Set rest density for all particles
            ghostBottom->setPressure(0.0f); // Set pressure to 0 for ghost particles
            ghostParticles.push_back(ghostBottom);

            //top
            auto ghostTop = std::make_shared<Particle>();
            ghostTop->setPosition(Vec3(x, domainSize + layer * offset, 0.0f));
            ghostTop->setMass(ghostMass); // Set mass to 1.0 for all particles
            ghostTop->setSmoothingLength(cellSize); // Set smoothing length based on spacing
            ghostNeighbors = findNeighbors(*ghostTop, particles, ghostTop->getSmoothingLength());
            ghostTop->setDensity(calculateDensity(ghostTop, ghostNeighbors)); // Set rest density for all particles
            ghostTop->setPressure(0.0f); // Set pressure to 0 for ghost particles
            ghostParticles.push_back(ghostTop);

            //left
            auto ghostLeft = std::make_shared<Particle>();
            ghostLeft->setPosition(Vec3(-layer * offset, x, 0.0f));
            ghostLeft->setMass(ghostMass); // Set mass to 1.0 for all particles
            ghostLeft->setSmoothingLength(cellSize); // Set smoothing length based on spacing
            ghostNeighbors = findNeighbors(*ghostLeft, particles, ghostLeft->getSmoothingLength());
            ghostLeft->setDensity(calculateDensity(ghostLeft, ghostNeighbors));
            ghostLeft->setPressure(0.0f); // Set pressure to 0 for ghost particles
            ghostParticles.push_back(ghostLeft);

            //right
            auto ghostRight = std::make_shared<Particle>();
            ghostRight->setPosition(Vec3(domainSize + layer * offset, x, 0.0f));
            ghostRight->setMass(ghostMass); // Set mass to 1.0 for all particles
            ghostRight->setSmoothingLength(cellSize); // Set smoothing length based on spacing
            ghostNeighbors = findNeighbors(*ghostRight, particles, ghostRight->getSmoothingLength());
            ghostRight->setDensity(calculateDensity(ghostRight, ghostNeighbors));
            ghostRight->setPressure(0.0f); // Set pressure to 0 for ghost particles
            ghostParticles.push_back(ghostRight);



        }
        std::cout << "Initialized " << ghostParticles.size() << " bottom ghost particles\n";
    }
    //add ghost particles outside the boundary evenly spaced
    // for (int i = 0; i < grid_size; ++i) {
    //     float h = spacing * 1.5f;
    //     float inside = spacing * 0.5f;
    //     float x = i * spacing;
    //     float y = -spacing * 1.5f;
    //     auto ghostParticle = std::make_shared<Particle>();
    //     ghostParticle->setPosition(Vec3(x, -inside, 0.0f));
    //     ghostParticle->setMass(1.0f); // Set mass to 1.0 for all particles
    //     ghostParticle->setSmoothingLength(spacing*1.5f); // Set smoothing length based on spacing
    //     ghostParticle->setDensity(restDensity); // Set rest density for all particles
    //     ghostParticle->setPressure(0.0f); // Set pressure to 0 for ghost particles
    //     ghostParticles.push_back(ghostParticle);
    //     std::cout << "Ghost Particle Position: " << ghostParticle->getPosition()[0] << ", " << ghostParticle->getPosition()[1] << std::endl;
    // }


}

void Simulation::updateParticles(const float dt) {
    // Build the spatial map for the current particles
    buildSpatialMap(particles);

    //reset neighbor list
    for (auto &neighbors : neighborList) {
        neighbors.clear();
    }

    size_t index = 0;
    //Update all particles with their density before doing anything else.
    for (auto &particle : particles) {
        particle->reset();
        // Find neighbors
        std::vector<std::shared_ptr<Particle>> neighbors = findNeighbors(*particle, particles, particle->getSmoothingLength());
        // Calculate density for the particle based on neighbors
        const auto density = calculateDensity(particle, neighbors);
        particle->setDensity(density);
        if (density > maxDensity) {
            maxDensity = density;
        }
        //Cache neighbors for later use to avoid recalculating
        //std::cout << "Density: " << particle->getDensity() << std::endl;
        neighborList[index] = neighbors;

        ++index;
    }

    //Now update pressure for all particles.
    for (auto &particle : particles) {
        //update the maximum density of each particle for accurate coloring
        particle->setMaxDensity(maxDensity);

        float pressure = calculatePressure(particle);

        //float pressure = stiffness * (particle->getDensity() - restDensity);
        //pressure = std::clamp(pressure, -1000.0f, 1000.0f);  // You can tune this max value
        particle->setPressure(pressure);

    }
    // Update ghost boundary particles’ densities and pressures
    for (auto &ghost : ghostParticles) {
        auto neighborsGhost = findNeighbors(*ghost, particles, ghost->getSmoothingLength());
        float ghostDensity = calculateDensity(ghost, neighborsGhost);
        ghostDensity = std::max(ghostDensity, restDensity);  // or even 1e-3f
        ghost->setDensity(ghostDensity);
        float pressure = calculatePressure(ghost);
        //clamp to avoid negative pressure
        //pressure = std::max(0.0f, pressure);  // You can tune this max value
        ghost->setPressure(pressure);
    }
    //Now update forces for all particles.
    index = 0;
    for (auto &particle : particles) {
        // Find neighbors
        const auto &neighbors = neighborList[index];
        // Calculate forces based on neighbors
        for (const auto &neighbor : neighbors) {
            if (particle == neighbor) continue; // Skip self
            float distance = particle->getNeighborDistance(*neighbor);
            if (distance < particle->getSmoothingLength()) {
                // Calculate force based on pressure and density differences

                //Test spiky kernel here
                Vec3 r_vec = particle->getPosition() - neighbor->getPosition();
                float h = particle->getSmoothingLength();

                //Compute spiky gradient for pressure force
                Vec3 spiky = spikyGradient(r_vec, h);

                //Create pressure term
                const float pressure = (particle->getPressure() + neighbor->getPressure()) / (2.0f * neighbor->getDensity());

                //Calculate pressure force
                Vec3 pressureForce = neighbor->getMass() * pressure * spiky;

                //Add viscosity force using laplacian
                Vec3 laplacian = viscosityLaplacian(r_vec, h);
                //velocity difference
                Vec3 velocityDiff = neighbor->getVelocity() - particle->getVelocity();
                //Calculate viscosity force

                Vec3 viscosityForce = viscosity * neighbor->getMass()/neighbor->getDensity() * laplacian * velocityDiff;
                Vec3 force = pressureForce + viscosityForce;

                particle->addForce(force);
            }

        }
        //Apply gravity
        particle->addForce(Vec3(0, particle->getMass() * gravity, 0));

        // Apply drag
        const Vec3 drag = -0.1f * particle->getVelocity();
        particle->addForce(drag);
        index++;
    }
    //Now update all particles with their forces
    for (auto &particle : particles) {
        //std::cout << "ρ: " << particle->getDensity() << "  P: " << particle->getPressure() << std::endl;
        particle->update(dt);
    }
}

float Simulation::calculatePressure(const std::shared_ptr<Particle>& particle) const
 {
    constexpr float y = 7.0f;
    constexpr float c0 = 150;
    const float B = (c0 * c0 * restDensity)/ y;
    float pressure = B* (std::pow(particle->getDensity() / restDensity, y)-1.0f);
    return pressure;
}

void Simulation::calibrateRestDensity() {
    float sum = 0;
    for (size_t i = 0; i < particles.size(); ++i) {
        auto &p = particles[i];
        auto nbrs = findNeighbors(*p, particles, p->getSmoothingLength());
        sum += calculateDensity(p, nbrs);
    }
    restDensity = sum / particles.size();
    std::cout << "Calibrated restDensity = " << restDensity << "\n";
}
