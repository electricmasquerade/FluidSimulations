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
    constexpr unsigned windowWidth = 800;
    sf::RenderWindow window(sf::VideoMode({windowWidth, windowWidth}), "ImGui + SFML = <3");
    window.setFramerateLimit(120);
    ImGui::SFML::Init(window);
    static float accumulator = 0.0f;
    const float physicsDt = 1.0f / 60.0f;

    // Seed the random number generator
    std::srand(static_cast<unsigned int>(std::time(nullptr)));


    // Create many particles for testing
    constexpr int numParticles = 2000; // You can adjust this number as needed.
    constexpr float domainSize = static_cast<int>(windowWidth);
    float spacing = domainSize / std::sqrt(static_cast<float>(numParticles));
    float smoothingLength = 2 * spacing;
    float cellSize = smoothingLength;
    constexpr float stiffness = 10;
    constexpr float restDensity = 1.0f;
    constexpr float viscosity = 100.0f;

    std::vector<std::shared_ptr<Particle>> particles;
    particles.reserve(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        //Create particles with default values for now to put into simulation
        Particle particle;
        //particle.setPosition(Vec3(rand() % window.getSize().x, rand() % window.getSize().y, 0.0f));
        particle.setVelocity(Vec3((rand() % 50), (rand() % 50), 0.0f));
        particle.setSmoothingLength(smoothingLength);
        particle.setMaxDensity(numParticles/20);
        particles.push_back(std::make_shared<Particle>(particle));

    }
    Simulation simulation(particles, domainSize, cellSize, stiffness, restDensity, viscosity);
    RenderLayer renderLayer(particles);

    //arrange particles into a ball in the center
    float centerX = windowWidth / 2.0f;
    float centerY = windowWidth / 2.0f;
    float radius = windowWidth / 4.0f;  // Use 1/4 of window width as radius

    for (size_t i = 0; i < particles.size(); ++i) {
        // Calculate angle for uniform distribution
        float angle = (2.0f * M_PI * i) / particles.size();

        // Calculate random radius (smaller than max radius) for more natural distribution
        float randRadius = radius * std::sqrt(static_cast<float>(rand()) / RAND_MAX);

        // Calculate position using polar coordinates
        float x = centerX + randRadius * std::cos(angle);
        float y = centerY + randRadius * std::sin(angle);

        particles[i]->setPosition(Vec3(x, y, 0.0f));
        particles[i]->setVelocity(Vec3((rand() % 100), (rand() % 100), 0.0f));
        //particles[i]->setVelocity(Vec3(0.0f, 0.0f, 0.0f)); // Start with zero velocity
    }


    //simulation.calibrateRestDensity();
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
        //ImGui::SetNextWindowSize(ImVec2(250, 120));
        ImGui::Begin("Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

        ImGui::Checkbox("Running", &isRunning);

        if (!isRunning) {
            // When paused, allow control of various variables
            static float gravity = 10.0f;
            ImGui::SliderFloat("Gravity", &gravity, -10.0f, 10.0f);
            simulation.setGravity(gravity);
        }


        // if (ImGui::Button("Change Color of Particle 0")) {
        //     // Change the first particle's color to a random color
        //     particles[0]->setColor(std::vector<int>({
        //         rand() % 256, rand() % 256, rand() % 256
        //     }));
        // }
        //add slider for gravity, linked to simulation->set gravity
        ImVec2 controls_pos = ImGui::GetWindowSize();
        ImGui::End();

        // ImGui window for displaying simulation data (e.g., frame rate)
        //adjusts position based on size of a previous window, i.e. get size of controls



        ImGui::SetNextWindowPos(ImVec2(controls_pos.x, 0));
        ImGui::SetNextWindowBgAlpha(0.5f);
        //ImGui::SetNextWindowSize(ImVec2(200, 100));
        ImGui::Begin("Simulation Data", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
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