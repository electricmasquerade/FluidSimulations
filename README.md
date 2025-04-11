# FluidSimulation

FluidSimulation is a 2D fluid simulation project that demonstrates Smoothed Particle Hydrodynamics (SPH) techniques. The project is written in C++ and utilizes SFML for rendering along with ImGui for real-time controls and debugging.

## Overview

This project is designed to simulate fluid behavior using particle-based methods (SPH). It includes:
- **Particle Dynamics:** Simulation of fluid flow using particles.
- **Efficient Rendering:** Batch rendering of particles using an SFML Vertex Array for high performance.
- **Spatial Partitioning:** Implementation of a uniform grid (spatial map) for rapid neighbor searches, which is key for SPH calculations.
- **Real-Time Controls:** Integration with ImGui to adjust simulation parameters on the fly.

## Features

- **SPH Simulation:** Utilizes SPH techniques to calculate density, pressure, and forces between particles.
- **Batch Rendering:** Renders large numbers of particles efficiently.
- **Spatial Partitioning:** Uses a uniform grid to accelerate neighbor search operations.
- **Modular Design:** The simulation, rendering, and utility functions are organized into separate classes to facilitate future enhancements.

## Requirements

- **C++17** or later
- **SFML** (version 3.0 or higher)
- **ImGui** (version 1.91.9b or higher)
- **ImGui-SFML** for UI integration, version 3.0 or higher
- **CMake** for building the project

## Build Instructions

1. Clone the repository.
2. Create a build directory and run CMake to configure the project.
3. Build the project using your preferred build system (e.g., via CLion or the command line).

Example using the command line:

```bash
mkdir build
cd build
cmake ..
make
./FluidSimulation
