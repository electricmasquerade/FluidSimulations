
#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include "Particle.h"
#include <cmath>
#include <unordered_map>

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
    float smoothingFunction(float r, float h);
    void buildSpatialMap(const std::vector<std::shared_ptr<Particle>> &particles);
    std::vector<std::shared_ptr<Particle>> findNeighbors(const Particle &particle, const std::vector<std::shared_ptr<Particle>> &particles, float radius);
    void updateParticles(float dt);

private:
    std::vector<std::shared_ptr<Particle>> particles;
    float cellSize{20.0f}; // Size of each cell in the grid, replace with smoothing length later
    std::unordered_map<CellKey, std::vector<std::shared_ptr<Particle>>, CellKeyHash> spatialMap; // Hash map for spatial partitioning
    //Create a struct to store the keys for the particle grid during spatial partitioning


};




#endif //SIMULATION_H
