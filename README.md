# Black Hole Ray Tracer

A C++ program that visually simulates a non-rotating (Schwarzschild) black hole using ray tracing techniques. This project demonstrates gravitational lensing effects by tracing light rays through curved spacetime.

## Features

- **Physically Accurate Ray Tracing**: Implements geodesic equations for light paths in Schwarzschild metric
- **Gravitational Lensing**: Visualizes how gravity bends light around massive objects
- **Event Horizon Rendering**: Clearly displays the Schwarzschild radius
- **Background Distortion**: Shows lensing effects on a checkerboard or starfield pattern
- **Standard Library Only**: No external dependencies required

## Physics Background

This simulation uses the Schwarzschild metric to describe spacetime curvature around a non-rotating black hole:

```
ds² = -(1 - rs/r)dt² + (1 - rs/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)
```

Where `rs = 2GM/c²` is the Schwarzschild radius (event horizon).

## Building the Project

### Prerequisites

- C++17 compatible compiler (GCC 7+, Clang 5+, MSVC 2017+)
- CMake 3.10 or higher (optional, but recommended)

### Build with CMake

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

### Build with g++ directly

```bash
g++ -std=c++17 -O3 -o black-hole-raytracer src/main.cpp src/raytracer.cpp src/vector3d.cpp -Iinclude
```

## Running the Simulation

```bash
./black-hole-raytracer [output_file.ppm] [width] [height]
```

**Examples:**

```bash
# Default settings (800x600)
./black-hole-raytracer

# Custom resolution
./black-hole-raytracer output.ppm 1920 1080

# Quick test
./black-hole-raytracer test.ppm 400 300
```

## Output

The program generates a PPM (Portable Pixmap) image file that can be viewed with most image viewers or converted to other formats using tools like ImageMagick:

```bash
convert output.ppm output.png
```

## Project Structure

```
black-hole-raytracer/
├── src/                 # Source files
│   ├── main.cpp        # Entry point
│   ├── raytracer.cpp   # Ray tracing engine
│   └── vector3d.cpp    # 3D vector math
├── include/            # Header files
│   ├── raytracer.h
│   ├── vector3d.h
│   └── constants.h
├── tests/              # Unit tests
├── scripts/            # Build and utility scripts
├── output/             # Generated images
├── docs/               # Documentation
├── CMakeLists.txt      # CMake configuration
└── README.md           # This file
```

## Configuration

You can modify simulation parameters in `include/constants.h`:

- `BLACK_HOLE_MASS`: Mass of the black hole
- `CAMERA_DISTANCE`: Distance from black hole to camera
- `MAX_ITERATIONS`: Ray tracing accuracy
- `STEP_SIZE`: Integration step size

## Technical Details

### Ray Tracing Algorithm

1. For each pixel, cast a ray from the camera
2. Integrate the geodesic equations using RK4 method
3. Check for event horizon crossing (captured by black hole)
4. If ray escapes, sample the background at final direction
5. Render pixel color based on result

### Geodesic Integration

The code uses a 4th-order Runge-Kutta (RK4) integrator to solve the differential equations describing light paths in curved spacetime.

## Performance

Typical rendering times on a modern CPU:
- 800x600: ~30-60 seconds
- 1920x1080: ~3-5 minutes

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License - See LICENSE file for details

## References

- Schwarzschild, K. (1916). "On the Gravitational Field of a Mass Point"
- Chandrasekhar, S. (1983). "The Mathematical Theory of Black Holes"
- James et al. (2015). "Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar"

## Next Steps

After mastering this project, check out the companion **N-Body Simulation** project for orbital dynamics around black holes!
