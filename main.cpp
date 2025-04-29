#include "imgui.h"
#include "imgui-SFML.h"
#include "Particle.h"

#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/System/Clock.hpp>
#include <SFML/Window/Event.hpp>
#include "RenderLayer.h"

#include <vector>
#include <cstdlib>  // for rand() and srand()
#include <ctime>    // for time()

#include "Simulation.h"

int main() {
    // Create the window and set up ImGui-SFML
    constexpr unsigned windowWidth = 500;
    sf::RenderWindow window(sf::VideoMode({windowWidth, windowWidth}), "ImGui + SFML = <3");
    window.setFramerateLimit(120);
    ImGui::SFML::Init(window);

    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));


    // Create many particles for testing
    constexpr int numParticles = 2000; // You can adjust this number as needed.
    constexpr float domainSize = static_cast<int>(windowWidth);
    float spacing = domainSize / std::sqrt(static_cast<float>(numParticles));
    float smoothingLength = 1.5f * spacing;
    float cellSize = smoothingLength;
    constexpr float stiffness = 100;
    constexpr float restDensity = 1.0f;
    constexpr float viscosity = 20.0f;

    std::vector<std::shared_ptr<Particle>> particles;
    particles.reserve(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        //Create particles with default values for now to put into simulation
        Particle particle;
        //particle.setPosition(Vec3(rand() % window.getSize().x, rand() % window.getSize().y, 0.0f));
        particle.setVelocity(Vec3((rand() % 100), (rand() % 100), 0.0f));
        particle.setSmoothingLength(smoothingLength);
        particles.push_back(std::make_shared<Particle>(particle));

    }
    Simulation simulation(particles, domainSize, cellSize, stiffness, restDensity, viscosity);
    RenderLayer renderLayer(particles);



    sf::Clock deltaClock;
    static bool isRunning = true;
    // Variables to control the first particle when paused
    static float sliderX = particles[0]->getPosition()[0];
    static float sliderY = particles[0]->getPosition()[1];

    while (window.isOpen()) {
        while (const auto event = window.pollEvent()) {
            ImGui::SFML::ProcessEvent(window, *event);

            if (event->is<sf::Event::Closed>()) {
                window.close();
            }
            if (const auto* resized = event->getIf<sf::Event::Resized>()) {
                // Update the view to the new size of the window
                sf::FloatRect visibleArea({0.f, 0.f}, sf::Vector2f(resized->size));
                window.setView(sf::View(visibleArea));
            }
        }

        sf::Time deltaTime = deltaClock.restart();
        ImGui::SFML::Update(window, deltaTime);

        // ImGui window for simulation controls
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowBgAlpha(0.5f);
        ImGui::SetNextWindowSize(ImVec2(250, 120));
        ImGui::Begin("Controls");

        ImGui::Checkbox("Running", &isRunning);
        ImGui::Text("Particle 0 Position: (%.1f, %.1f)", particles[0]->getPosition()[0], particles[0]->getPosition()[1]);

        if (!isRunning) {
            // When paused, allow manual control of the first particle
            ImGui::SliderFloat("X Position", &sliderX, 0.f, static_cast<float>(window.getSize().x));
            ImGui::SliderFloat("Y Position", &sliderY, 0.f, static_cast<float>(window.getSize().y));
            particles[0]->setPosition(Vec3(sliderX, sliderY, 0.0f));
        } else {
            // Update all particles when simulation is running
            for (auto &p : particles) {
                p->update(deltaTime.asSeconds());
            }
        }

        if (ImGui::Button("Change Color of Particle 0")) {
            // Change the first particle's color to a random color
            particles[0]->setColor(std::vector<int>({
                rand() % 256, rand() % 256, rand() % 256
            }));
        }
        ImGui::End();

        // ImGui window for displaying simulation data (e.g., frame rate)
        ImGui::SetNextWindowPos(ImVec2(260, 0));
        ImGui::SetNextWindowBgAlpha(0.5f);
        ImGui::SetNextWindowSize(ImVec2(200, 100));
        ImGui::Begin("Simulation Data");
        ImGui::Text("Frame Rate: %.1f FPS", 1.0f / ImGui::GetIO().DeltaTime);
        ImGui::End();

        //Update simulation timestep
        if (isRunning) {
            simulation.updateParticles(deltaTime.asSeconds());
        }
        // Update particle shapes based on the current state of each particle
        renderLayer.updateShapes(particles);

        window.clear();
        renderLayer.drawParticles(window);
        ImGui::SFML::Render(window);
        window.display();
    }

    ImGui::SFML::Shutdown();
    return 0;
}
