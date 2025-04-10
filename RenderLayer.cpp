
#include "RenderLayer.h"

/**
 * @brief Initializes the particle shapes based on the provided particles.
 * @param particles The vector of particles to initialize.
 */
void RenderLayer::initParticles(const std::vector<Particle> &particles) {
    for (const auto &particle : particles) {
        //
        sf::CircleShape shape(particle.getMass());
        // Set the shape's radius based on the particle's mass
        shape.setRadius(particle.getMass());
        shape.setOrigin({shape.getRadius(), shape.getRadius()});
        shape.setPosition({particle.getPosition()[0], particle.getPosition()[1]});

        //set a variable for a color based on the particle's color
        const sf::Color particleColor(
            particle.getColor()[0],
            particle.getColor()[1],
            particle.getColor()[2]
        );
        shape.setFillColor(particleColor);
        //Add the shapes to the list of particle objects to update and render
        particleShapes.push_back(shape);

    }
}


/**
 * @brief Updates the shapes of the particles based on their current positions.
 * @param particles The vector of particles to update.
 */
void RenderLayer::updateShapes(const std::vector<Particle> &particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto &particle = particles[i];
        particleShapes[i].setPosition({particle.getPosition()[0], particle.getPosition()[1]});
        // Set the shape's color based on the particle's color
        const sf::Color particleColor(
            particle.getColor()[0],
            particle.getColor()[1],
            particle.getColor()[2]
        );
        particleShapes[i].setFillColor(particleColor);
    }
}

/**
 * @brief Draws the particle shapes to the provided window.
 * @param window The SFML window to draw the particles on.
 */
void RenderLayer::drawParticles(sf::RenderWindow &window) const {
    for (const auto &shape : particleShapes) {
        window.draw(shape);
    }
}
