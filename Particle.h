
#ifndef PARTICLE_H
#define PARTICLE_H

#include "Vec3.h"
#include <vector>

/**
 * \class Particle
 * \brief Handles the individual particle data for the simulation.
 *
 * The Particle class represents a simulation particle with physical properties
 * such as mass, position, velocity, force, acceleration, density, pressure,
 * and smoothing length. It provides methods to update and reset its state, as
 * well as applying forces that affect its motion and other simulation properties.
 */

class Particle {
public:
    Particle();
    Particle(const float mass, const Vec3 &position, const Vec3 &velocity, const float smoothingLength = 1.0f)
        : mass(mass), position(position), velocity(velocity), smoothingLength(smoothingLength) {
        force = Vec3(0, 0, 0);
        acceleration = Vec3(0, 0, 0);
        density = 0;
        pressure = 0;
        color = {255,255,255}; // Default color white
    }

    void reset();
    void update(float dt);
    void addForce(const Vec3 &force);
    [[nodiscard]] float getNeighborDistance(const Particle &other) const {
        return(position - other.position).length();
    }

    //Getters
    [[nodiscard]] Vec3 getPosition() const { return position; }
    [[nodiscard]] Vec3 getVelocity() const { return velocity; }
    [[nodiscard]] Vec3 getForce() const { return force; }
    [[nodiscard]] Vec3 getAcceleration() const { return acceleration; }
    [[nodiscard]] float getDensity() const { return density; }
    [[nodiscard]] float getPressure() const { return pressure; }
    [[nodiscard]] float getMass() const { return mass; }
    [[nodiscard]] float getSmoothingLength() const { return smoothingLength; }

    [[nodiscard]] std::vector<int> getColor() const {
        return color;
    }


    //Normalized getters
    [[nodiscard]] Vec3 getNormalizedVelocity() const {
        return unitVector(velocity);
    }
    [[nodiscard]] Vec3 getNormalizedForce() const {
        return unitVector(force);
    }
    [[nodiscard]] Vec3 getNormalizedAcceleration() const {
        return unitVector(acceleration);
    }

    //Setters
    void setMass(const float mass) { this->mass = mass; }
    void setPosition(const Vec3& position){ this->position = position; }
    void setDensity(const float density){ this->density = density; }
    void setVelocity(const Vec3& velocity){ this->velocity = velocity; }
    void setForce(const Vec3& force){ this->force = force; }
    void setAcceleration(const Vec3& acceleration){ this->acceleration = acceleration; }
    void setPressure(const float pressure){ this->pressure = pressure; }
    void setSmoothingLength(const float smoothingLength){ this->smoothingLength = smoothingLength; }
    void setColor(const std::vector<int>& color) { this->color = color; }


private:
    float mass{1};
    Vec3 position{};
    Vec3 velocity{};
    Vec3 force{};
    Vec3 acceleration{};
    std::vector<int> color{};

    //Specific SPH simulation variables
    float density{};
    float pressure{};
    float smoothingLength{};

    //Maximum speed for particles
    float maxParticleSpeed = 50.0f; // Maximum speed for particles

};



#endif //PARTICLE_H
