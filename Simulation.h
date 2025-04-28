
#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include "Particle.h"
#include <cmath>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <iostream>

//Create hash functions and key structs for spatial partitioning
struct CellKey {
    int x, y;
    bool operator==(const CellKey &other) const {
        return x == other.x && y == other.y;
    }
};

//create hash function for cell key
struct CellKeyHash {
    std::size_t operator()(const CellKey &key) const {
        return std::hash<int>()(key.x) ^ std::hash<int>()(key.y);
    }
};

class Simulation {

public:
    Simulation() = default;
    explicit Simulation(std::vector<std::shared_ptr<Particle>> &particles) : particles(particles) {
        // Initialize the simulation with the provided particles
        initParticles(particles);
    }
    Simulation(std::vector<std::shared_ptr<Particle>> &particles, const float domainSize, const float cellSize, const float stiffness,
        const float restDensity) : domainSize(domainSize), cellSize(cellSize), stiffness(stiffness), restDensity(restDensity), particles(particles) {
        // Initialize the simulation with the provided particles and custom parameters
        initParticles(particles);
    }
    ~Simulation() = default;
    float smoothingPoly6(float r, float h);
    Vec3 spikyGradient(const Vec3 &r_vec, float h);
    void buildSpatialMap(const std::vector<std::shared_ptr<Particle>> &particles);
    float calculateDensity(const std::shared_ptr<Particle> &particle, const std::vector<std::shared_ptr<Particle>> &neighbors);
    std::vector<std::shared_ptr<Particle>> findNeighbors(const Particle &particle, const std::vector<std::shared_ptr<Particle>> &particles, float radius);
    void initParticles(std::vector<std::shared_ptr<Particle>> &particles);
    void updateParticles(float dt);

private:
    float gravity{9.81f}; // Gravitational acceleration, always constant
    float domainSize{500.0f}; // Size of the simulation domain
    float cellSize{20.0f}; // Size of each cell in the grid, replace with smoothing length later
    float stiffness{1.0f}; // Stiffness of the fluid, tune for stability
    float restDensity{1.0f}; // Rest density of the fluid, calculate during init particles.
    std::vector<std::shared_ptr<Particle>> particles;
    std::vector<std::shared_ptr<Particle>> ghostParticles; //for boundary handling
    std::unordered_map<CellKey, std::vector<std::shared_ptr<Particle>>, CellKeyHash> spatialMap; // Hash map for spatial partitioning
    //Create a struct to store the keys for the particle grid during spatial partitioning


};




#endif //SIMULATION_H
