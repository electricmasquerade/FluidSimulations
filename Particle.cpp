
#include "Particle.h"

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
    acceleration = force / mass;
    velocity += dt * acceleration;
    position += dt * velocity;
}

void Particle::addForce(const Vec3 &force) {
    this->force += force;
}
