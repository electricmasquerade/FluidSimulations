#include "imgui.h"
#include "imgui-SFML.h"
#include "Particle.h"

#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>

#include "RenderLayer.h"

int main() {
    sf::RenderWindow window(sf::VideoMode({640, 480}), "ImGui + SFML = <3");
    window.setFramerateLimit(120);
    ImGui::SFML::Init(window);

    RenderLayer renderLayer;

    //Create a single particle for testing purposes
    Particle particle(1.0f, Vec3(50, 50, 0), Vec3(200, 175, 0), 1.0f);

    renderLayer.initParticles(std::vector<Particle>{particle});


    sf::Clock deltaClock;
    while (window.isOpen()) {
        while (const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);

            if (event->is<sf::Event::Closed>()) {
                window.close();

            }
            if (const auto* resized = event->getIf<sf::Event::Resized>())
            {
                // update the view to the new size of the window
                sf::FloatRect visibleArea({0.f, 0.f}, sf::Vector2f(resized->size));
                window.setView(sf::View(visibleArea));
            }

        }
        sf::Time deltaTime = deltaClock.restart();
        ImGui::SFML::Update(window, deltaTime);



        //Window for sliders and controls
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowBgAlpha(0.5f);
        ImGui::SetNextWindowSize(ImVec2(250, 100));
        ImGui::Begin("Controls");

        static bool isRunning = true;
        ImGui::Checkbox("Running", &isRunning);
        sf::CircleShape shape = renderLayer.getParticleShapes()[0];
        static float posX = shape.getPosition().x;
        static float posY = shape.getPosition().y;
        ImGui::Text("Particle Position: (%.1f, %.1f)", posX, posY);
        if (!isRunning) {
            ImGui::SliderFloat("X Position", &posX, 0.f, window.getSize().x);
            ImGui::SliderFloat("Y Position", &posY, 0.f, window.getSize().y);
            particle.setPosition(Vec3(posX, posY, 0.0f));
            shape.setPosition(sf::Vector2f(particle.getPosition()[0], particle.getPosition()[1]));
        }
        else {
            particle.update(deltaTime.asSeconds());
            shape.setPosition(sf::Vector2f(particle.getPosition()[0], particle.getPosition()[1]));
            posX = particle.getPosition()[0];
            posY = particle.getPosition()[1];

            // 2. Update the collision conditions in the "else" branch when isRunning is true.
            // Find the block where you currently handle collisions; replace it with the code below:

            float radius = shape.getRadius();
            if (posX - radius < 0) {
                particle.setPosition(Vec3(radius, posY, 0));
                particle.setVelocity(Vec3(-particle.getVelocity()[0], particle.getVelocity()[1], 0));
            }
            else if (posX + radius > window.getSize().x) {
                particle.setPosition(Vec3(window.getSize().x - radius, posY, 0));
                particle.setVelocity(Vec3(-particle.getVelocity()[0], particle.getVelocity()[1], 0));
            }
            if (posY - radius < 0) {
                particle.setPosition(Vec3(posX, radius, 0));
                particle.setVelocity(Vec3(particle.getVelocity()[0], -particle.getVelocity()[1], 0));
            }
            else if (posY + radius > window.getSize().y) {
                particle.setPosition(Vec3(posX, window.getSize().y - radius, 0));
                particle.setVelocity(Vec3(particle.getVelocity()[0], -particle.getVelocity()[1], 0));
            }
        }


        if (ImGui::Button("Change Color")) {
            shape.setFillColor(sf::Color(rand() % 256, rand() % 256, rand() % 256));
        }
        ImGui::End();

        //Window for displaying data about the simulation
        ImGui::SetNextWindowPos(ImVec2(250, 0));
        ImGui::SetNextWindowBgAlpha(0.5f);
        ImGui::SetNextWindowSize(ImVec2(200, 100));
        ImGui::Begin("Simulation Data");
        ImGui::Text("Frame Rate: %.1f FPS", 1.0f / ImGui::GetIO().DeltaTime);
        ImGui::End();



        window.clear();
        window.draw(shape);
        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
}
