# asu-p3: Eulerian Fluid Simulation

This repository contains an Eulerian fluid simulation used to simulate smoke, implemented in C++ as a class project for the Animación y Simulación Avanzada 1 course at URJC. 

The core mathematical logic and fluid dynamics are entirely coded in `src/Fluid2Exercise.cpp`. The foundational engine code, window management, and base utilities were provided by the course professor.

## Overview
This project simulates smoke moving through a 2D grid using the Eulerian approach to fluid dynamics. The numerical methods follow the standard stable fluids algorithm, applying the Navier-Stokes equations for incompressible fluid flow over discrete time steps.

### Mathematical Implementation Details
The simulation steps located in `Fluid2Exercise.cpp` implement the following physics concepts:
- **Advection:** Computes the transport of velocity and density (smoke ink) using Semi-Lagrangian advection.
- **External Forces:** Introduces external forces like gravity and continuous fluid emission into the system.
- **Viscosity:** Applies an explicit finite difference approximation of the Laplacian to diffuse velocity.
- **Pressure Projection:** Ensures fluid incompressibility by computing the pressure field. This is accomplished by solving the resulting Poisson equation using a Preconditioned Conjugate Gradient (PCG) solver, followed by subtracting the pressure gradient from the intermediate velocity field.

## Dependencies
The project utilizes **CMake** for build configuration and **vcpkg** for dependency management. Primary external dependencies include:
- `OpenGL`
- `FreeGLUT`

## Build & Run Instructions

Since this project uses modern CMake and vcpkg, it is recommended to use an IDE with native support for both, such as **CLion** or **Visual Studio Code**, which will automatically configure the build environment.

To build manually via the command line:

```bash
# Clone the repository
git clone https://github.com/kiinami/asu-p3.git
cd asu-p3

# Configure the project with vcpkg integration
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Release

# Build the project
cmake --build build

# Run the executable
./build/asu-p3
```

## Results

<!-- Replace the placeholder below with an actual video or image capture of the simulation -->
![Smoke Simulation Capture](static/simulation.gif)
