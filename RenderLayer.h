
#ifndef RENDERLAYER_H
#define RENDERLAYER_H

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>
#include <vector>

#include "Particle.h"


//this is an interface layer designed to handle all the rendering steps to draw particles.
class RenderLayer {
public:
    explicit RenderLayer(const std::vector<std::shared_ptr<Particle>>& particles) {
        this->initParticles(particles);
    }
    ~RenderLayer() = default;

    void initParticles(const std::vector<std::shared_ptr<Particle>> &particles);
    void updateShapes(const std::vector<std::shared_ptr<Particle>> &particles);
    void drawParticles(sf::RenderWindow &window) const;
    // [[nodiscard]] const std::vector<sf::CircleShape> &getParticleShapes() const {
    //     return particleShapes;
    // }
private:
    //store a list of particles
    // std::vector<Particle> particles;
    //store a list of shapes that correspond to particles
    //std::vector<sf::CircleShape> particleShapes;
    sf::VertexArray particle_vertices;


};



#endif //RENDERLAYER_H
