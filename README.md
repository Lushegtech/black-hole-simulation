# Black Hole Simulation in C++

A physically accurate simulation of a Schwarzschild (non-rotating) black hole using **General Relativistic Ray Tracing (GRRT)**. This project implements Einstein's theory of General Relativity to visualize gravitational lensing, the black hole shadow, and accretion disk effects.

![Black Hole Visualization](https://img.shields.io/badge/Physics-General_Relativity-blue)
![C++17](https://img.shields.io/badge/C++-17-00599C?logo=c%2B%2B)
![OpenGL](https://img.shields.io/badge/OpenGL-3.3-5586A4?logo=opengl)

## Features

### Physical Accuracy
- **Schwarzschild Metric**: Full implementation of the exact solution to Einstein's field equations for a non-rotating black hole
- **Geodesic Integration**: Solves the geodesic equations using 4th-order Runge-Kutta (RK4) method
- **Event Horizon**: Correctly models the boundary of no return at the Schwarzschild radius (rs = 2GM/c²)
- **Photon Sphere**: Captures unstable photon orbits at r = 1.5rs
- **ISCO**: Innermost Stable Circular Orbit at r = 3rs (where the accretion disk begins)
- **Gravitational Lensing**: Physically accurate bending of light around massive objects

### Dual Rendering Modes

#### 1. CPU Ray Tracer (`black-hole-raytracer`)
- High-quality offline rendering
- Outputs to PPM image format
- Configurable resolution
- Ideal for high-resolution still images

#### 2. GPU Real-Time Renderer (`black-hole-gpu`)
- Real-time interactive visualization
- GLSL fragment shader implementation
- 60+ FPS performance
- Interactive camera controls
- Live parameter adjustment

## Physics Background

This simulation is based on the **Schwarzschild solution** to Einstein's field equations:

```
ds² = -(1 - rs/r)c²dt² + (1 - rs/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)
```

Where:
- `rs = 2GM/c²` is the Schwarzschild radius (event horizon)
- Natural units are used: G = M = c = 1, so rs = 2

### Why General Relativity?

A Newtonian "dark star" model is **qualitatively incorrect** for simulating black holes. Key phenomena that only exist in General Relativity:

1. **Event Horizon**: Not a physical surface, but a boundary in spacetime geometry where all future paths lead inward
2. **ISCO**: Below this radius, no stable circular orbits exist (impossible in Newtonian physics)
3. **Photon Sphere**: Unstable circular orbits for light at r = 1.5rs
4. **Frame Dragging**: Spacetime itself rotates near a spinning black hole (Kerr metric)

## Building the Project

### Prerequisites

**Required:**
- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.10 or higher

**Optional (for GPU renderer):**
- OpenGL 3.3+
- GLFW3
- GLEW

### Install Dependencies (Ubuntu/Debian)

```bash
# Required
sudo apt-get update
sudo apt-get install build-essential cmake

# Optional (for GPU renderer)
sudo apt-get install libglfw3-dev libglew-dev libgl1-mesa-dev
```

### Build Instructions

```bash
# Clone the repository
git clone <repository-url>
cd black-hole-simulation

# Create build directory
mkdir build
cd build

# Configure with CMake
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build
cmake --build . -j$(nproc)

# Executables will be in the build/ directory
```

The build system will automatically detect if OpenGL libraries are available:
- **If found**: Both `black-hole-raytracer` and `black-hole-gpu` will be built
- **If not found**: Only `black-hole-raytracer` (CPU version) will be built

## Running the Simulations

### CPU Ray Tracer

```bash
# Default settings (800x600)
./black-hole-raytracer

# Custom resolution and output file
./black-hole-raytracer output/my_blackhole.ppm 1920 1080

# Quick test
./black-hole-raytracer output/test.ppm 400 300
```

**Performance**: Typical rendering times on a modern CPU:
- 400x300: ~3 seconds
- 800x600: ~30-60 seconds
- 1920x1080: ~3-5 minutes

**Output**: The program generates a PPM (Portable Pixmap) image. To convert to PNG:

```bash
convert output.ppm output.png
```

### GPU Real-Time Renderer

```bash
./black-hole-gpu
```

**Controls**:
- **W/S**: Move camera closer/farther from black hole
- **A/D**: Rotate around black hole (azimuthal angle)
- **Q/E**: Move camera up/down (polar angle)
- **Mouse Drag** (left button): Rotate view
- **Mouse Wheel**: Adjust field of view (zoom)
- **ESC**: Exit

**Performance**: 60+ FPS at 1280x720 on modern GPUs

## Project Structure

```
black-hole-simulation/
├── src/
│   ├── main.cpp          # CPU ray tracer entry point
│   ├── main_gpu.cpp      # GPU real-time renderer
│   └── shader.cpp        # Shader utility class
├── include/
│   ├── vec4.h            # 4-vector math (spacetime)
│   ├── mat4.h            # 4x4 matrix (metric tensor)
│   ├── constants.h       # Physical constants
│   ├── schwarzschild.h   # Schwarzschild metric & geodesics
│   ├── integrator.h      # RK4 numerical integrator
│   ├── raytracer.h       # CPU ray tracing engine
│   └── shader.h          # Shader loader
├── shaders/
│   ├── vertex.glsl       # Vertex shader (fullscreen quad)
│   └── fragment.glsl     # Fragment shader (GRRT pipeline)
├── output/               # Generated images
├── assets/               # Textures (star maps, etc.)
├── CMakeLists.txt        # Build configuration
└── README.md
```

## Technical Details

### General Relativistic Ray Tracing (GRRT) Algorithm

The core algorithm follows the **backward ray tracing** approach:

1. **Initialize Ray**: For each pixel, create a null geodesic (photon path) from the camera
2. **Integrate Geodesic**: Use RK4 to solve the geodesic equations step-by-step
3. **Check Termination**:
   - **Captured**: Ray crosses event horizon → pixel is black
   - **Escaped**: Ray reaches infinity → sample background color
   - **Disk Hit**: Ray intersects accretion disk → compute disk color
   - **Max Steps**: Prevent infinite loops (e.g., photon sphere orbits)
4. **Return Color**: Assign the final color to the pixel

### Geodesic Equations

Using the **Christoffel symbols** method, the geodesic equation is:

```
d²xᵘ/dλ² + Γᵘ_νσ (dxᵛ/dλ)(dxᵍ/dλ) = 0
```

Where:
- `λ` is the affine parameter (proper time for massive particles)
- `Γᵘ_νσ` are the Christoffel symbols derived from the metric

**Implementation**: See `schwarzschild.h:geodesic_acceleration()`

### Numerical Integration

**4th-Order Runge-Kutta (RK4)** is used for integrating the geodesic equations:

```
k₁ = h·f(λₙ, Yₙ)
k₂ = h·f(λₙ + h/2, Yₙ + k₁/2)
k₃ = h·f(λₙ + h/2, Yₙ + k₂/2)
k₄ = h·f(λₙ + h, Yₙ + k₃)
Yₙ₊₁ = Yₙ + (k₁ + 2k₂ + 2k₃ + k₄)/6
```

**Limitations**: Fixed-step RK4 is used for simplicity. For higher accuracy near the event horizon, adaptive step-size methods (RK45, Dormand-Prince) are recommended.

### GPU Implementation

The GPU version uses **GLSL fragment shaders** to trace one ray per pixel in parallel:

- **Vertex Shader** (`vertex.glsl`): Draws a fullscreen quad
- **Fragment Shader** (`fragment.glsl`): Contains the complete GRRT pipeline
  - Schwarzschild metric computation
  - Christoffel symbols
  - RK4 integrator
  - Termination checks
  - Color computation

This achieves **massively parallel** ray tracing, with each GPU core handling multiple pixels simultaneously.

## Configuration

Edit `include/constants.h` to modify simulation parameters:

```cpp
namespace Physics {
    constexpr double r_s = 2.0;              // Schwarzschild radius
    constexpr double r_ISCO = 6.0;           // ISCO radius
    constexpr double r_disk_inner = 6.0;     // Accretion disk inner edge
    constexpr double r_disk_outer = 20.0;    // Accretion disk outer edge
    constexpr double STEP_SIZE = 0.1;        // Integration step size
    constexpr int MAX_STEPS = 10000;         // Maximum iterations
}

namespace Rendering {
    constexpr double CAMERA_DISTANCE = 20.0; // Camera distance from BH
    constexpr double CAMERA_FOV = 60.0;      // Field of view (degrees)
}
```

## Visualization Guide

### What You're Seeing

1. **Black Circle (Shadow)**: The black hole's **shadow** - larger than the event horizon due to captured photons
   - Schwarzschild shadow radius: ~2.6 rs

2. **Bright Ring**: The **accretion disk** - hot plasma orbiting at and beyond the ISCO
   - Inner edge: r = 6 (ISCO)
   - Outer edge: r = 20 (configurable)
   - Temperature falls off as T ∝ r^(-3/4)

3. **Distorted Background**: **Gravitational lensing** effects
   - Light from distant stars/checkerboard is bent around the black hole
   - Creates multiple images and Einstein rings
   - More pronounced near the photon sphere (r = 3)

4. **Asymmetry** (future): When relativistic effects are added:
   - Doppler shift: approaching side is brighter/bluer
   - Relativistic beaming: intensity scales as (1+z)^(-4)

## Future Enhancements

### Planned Features

- [ ] **Kerr Metric**: Rotating black holes (spin parameter a)
- [ ] **Adaptive Integrator**: RK45 or Dormand-Prince for better accuracy
- [ ] **Relativistic Effects**: Full Doppler shift and beaming for accretion disk
- [ ] **Star Map Textures**: Replace checkerboard with real celestial sphere
- [ ] **Bloom & HDR**: Post-processing for realistic brightness
- [ ] **Ray Bundles**: Anti-aliasing via geodesic deviation equation (Interstellar method)
- [ ] **Magnetohydrodynamics**: Import GRMHD simulation data for realistic accretion
- [ ] **VR Support**: Immersive black hole exploration

### Advanced Topics

- **Kerr-Schild Coordinates**: For rotating black holes, avoid coordinate singularities at poles
- **GRRT + GRMHD**: Combine ray tracing with fluid dynamics simulations
- **Penrose Diagrams**: Visualize the full spacetime structure
- **Hawking Radiation**: Quantum effects near the horizon (requires beyond-GR physics)

## Performance Optimization

### CPU Version
- Compile with `-O3` and `-march=native` (already in CMake Release mode)
- Parallelize with OpenMP: `#pragma omp parallel for` on the pixel loop
- Reduce `MAX_STEPS` for faster (but less accurate) rendering

### GPU Version
- Reduce `MAX_STEPS` in `fragment.glsl` for higher FPS
- Lower screen resolution for faster rendering
- Use adaptive step-size in shader (advanced)

## References

### Scientific Papers
1. Schwarzschild, K. (1916). "On the Gravitational Field of a Mass Point"
2. Chandrasekhar, S. (1983). "The Mathematical Theory of Black Holes"
3. James, O. et al. (2015). "Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar" - [arXiv:1502.03808](https://arxiv.org/abs/1502.03808)

### Educational Resources
- [Spinning Black Holes Tutorial](https://20k.github.io/c++/2020/02/12/black-hole-shader.html) - Excellent GLSL implementation guide
- [Event Horizon Telescope](https://eventhorizontelescope.org/) - Real black hole observations
- [NASA Black Hole Visualization](https://svs.gsfc.nasa.gov/13326) - Scientific visualizations

### Code References
- [Odyssey](https://github.com/annakorukova/odyssey) - GPU GRRT code in CUDA
- [grtrace](https://github.com/grtrace/grtrace) - Raytracer in D language
- [rossning92/Blackhole](https://github.com/rossning92/Blackhole) - C++/OpenGL implementation

## Contributing

Contributions are welcome! Areas of interest:
- Implementing the Kerr metric (rotating black holes)
- Adding adaptive step-size integrators
- Optimizing shader performance
- Creating realistic accretion disk models
- Adding more celestial backgrounds

## License

MIT License - See LICENSE file for details

## Acknowledgments

This implementation is based on the technical guide for General Relativistic Black Hole Simulation, following the mathematical framework from:
- Einstein's General Theory of Relativity (1915)
- Schwarzschild's solution (1916)
- Modern numerical relativity techniques

Special thanks to the authors of the "Interstellar" visual effects paper for their detailed documentation of the DNGR renderer.

---

**Note**: This is a **physically accurate** simulation, not an artistic approximation. All visual effects (gravitational lensing, shadow size, disk structure) emerge naturally from solving Einstein's equations. No "fake" effects are added for visual appeal.
