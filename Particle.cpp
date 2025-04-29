
#include "Particle.h"
#include "HelperFunctions.h"

Particle::Particle() {
    //Set almost all to zero by default
    this->position = Vec3(0, 0, 0);
    this->velocity = Vec3(0, 0, 0);
    this->force = Vec3(0, 0, 0);
    this->acceleration = Vec3(0, 0, 0);
    this->density = 0;
    this->pressure = 0;
    this->mass = 1; // Default mass
    this->smoothingLength = 1; // Default smoothing length

    this->color = {0, 0, 0};
}

void Particle::reset() {
    this->density = 0;
    this->pressure = 0;
    this->force = Vec3(0, 0, 0);

}

void Particle::update(const float dt) {
    constexpr float domainMin = 0.0f;
    constexpr float domainMax = 800.0f;
    constexpr float bounceDamping = 0.3f;
    constexpr float buffer = 1.0f; // Keep particle slightly away from edge
    acceleration = force / mass;
    velocity += dt * acceleration;
    position += dt * velocity;
    if (position[0] < domainMin) {
        position[0] = domainMin + buffer;
        velocity[0] *= -bounceDamping;
    }
    if (position[0] > domainMax) {
        position[0] = domainMax - buffer;
        velocity[0] *= -bounceDamping;
    }
    if (position[1] < domainMin) {
        position[1] = domainMin + buffer;
        velocity[1] *= -bounceDamping;
    }
    if (position[1] > domainMax) {
        position[1] = domainMax - buffer;
        velocity[1] *= -bounceDamping;
    }

    //print if out of range for debugging
    if (position.x() > domainMax || position.x() < domainMin) {
        std::cout << position.x() << std::endl;
    }
    if (position.y() > domainMax || position.y() < domainMin) {
        std::cout << position.y() << std::endl;
    }





    //update color based on density
    // Normalize density to a range of 0 to 1
    const float densityRange = maxDensity - minDensity;
    const float normalizedDensity = (density - minDensity) / densityRange;


    const float hue = 240 * (1 - normalizedDensity); // Blue to red
    color = HelperFunctions::HSVtoRGB(static_cast<int>(hue), 1, 1);
    // Set the color based on the speed

}

void Particle::addForce(const Vec3 &force) {
    this->force += force;
}
