#include "RenderLayer.h"

/**
 * @brief Initializes the particle shapes based on the provided particles.
 * @param particles The vector of particles to initialize.
 */
void RenderLayer::initParticles(const std::vector<std::shared_ptr<Particle>> &particles) {
    particle_vertices.clear();
    particle_vertices.setPrimitiveType(sf::PrimitiveType::Points);
    particle_vertices.resize(particles.size());
}


/**
 * @brief Updates the shapes of the particles based on their current positions.
 * @param particles The vector of particles to update.
 */
void RenderLayer::updateShapes(const std::vector<std::shared_ptr<Particle>> &particles) {
    if (particles.empty()) {
        return;
    }
    particle_vertices.resize(particles.size());

    for (std::size_t i = 0; i < particles.size(); ++i) {
        Vec3 pos = particles[i]->getPosition();
        particle_vertices[i].position = sf::Vector2f(pos[0], pos[1]);
        std::vector<int> color = particles[i]->getColor();
        particle_vertices[i].color = sf::Color(color[0], color[1], color[2]);

    }
}

/**
 * @brief Draws the particle shapes to the provided window.
 * @param window The SFML window to draw the particles on.
 */
void RenderLayer::drawParticles(sf::RenderWindow &window) const {
    window.draw(particle_vertices);
}
