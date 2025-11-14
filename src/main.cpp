/**
 * Black Hole Ray Tracer - Main Entry Point
 *
 * A physically accurate simulation of a Schwarzschild (non-rotating) black hole
 * using General Relativistic Ray Tracing (GRRT).
 *
 * This implementation follows the technical guide for simulating black holes
 * in C++ using the geodesic equations from Einstein's General Relativity.
 */

#include "raytracer.h"
#include "schwarzschild.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

/**
 * Write image to PPM file (Portable Pixmap format)
 */
void write_ppm(const std::string& filename,
               const std::vector<Color>& image,
               int width, int height) {
    std::ofstream file(filename, std::ios::binary);

    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // PPM header
    file << "P3\n" << width << " " << height << "\n255\n";

    // Write pixel data
    for (const auto& pixel : image) {
        int r = static_cast<int>(pixel.r * 255.0);
        int g = static_cast<int>(pixel.g * 255.0);
        int b = static_cast<int>(pixel.b * 255.0);

        // Clamp values
        r = std::max(0, std::min(255, r));
        g = std::max(0, std::min(255, g));
        b = std::max(0, std::min(255, b));

        file << r << " " << g << " " << b << "\n";
    }

    file.close();
    std::cout << "Image saved to " << filename << std::endl;
}

/**
 * Print simulation information
 */
void print_info() {
    std::cout << "========================================\n";
    std::cout << "  Black Hole Ray Tracer (Schwarzschild)\n";
    std::cout << "========================================\n\n";

    std::cout << "Physical Parameters:\n";
    std::cout << "  Schwarzschild radius (rs): " << Physics::r_s << "\n";
    std::cout << "  Event horizon: " << Physics::r_horizon << "\n";
    std::cout << "  Photon sphere: " << Physics::r_photon << "\n";
    std::cout << "  ISCO radius: " << Physics::r_ISCO << "\n";
    std::cout << "  Accretion disk: " << Physics::r_disk_inner
              << " to " << Physics::r_disk_outer << "\n\n";

    std::cout << "Numerical Parameters:\n";
    std::cout << "  Integration method: RK4\n";
    std::cout << "  Step size: " << Physics::STEP_SIZE << "\n";
    std::cout << "  Max steps: " << Physics::MAX_STEPS << "\n";
    std::cout << "  Escape radius: " << Physics::ESCAPE_RADIUS << "\n\n";
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    std::string output_file = "output.ppm";
    int width = Rendering::WIDTH;
    int height = Rendering::HEIGHT;

    if (argc > 1) output_file = argv[1];
    if (argc > 2) width = std::atoi(argv[2]);
    if (argc > 3) height = std::atoi(argv[3]);

    // Print information
    print_info();

    std::cout << "Render Settings:\n";
    std::cout << "  Output file: " << output_file << "\n";
    std::cout << "  Resolution: " << width << "x" << height << "\n";
    std::cout << "  Camera distance: " << Rendering::CAMERA_DISTANCE << "\n";
    std::cout << "  Field of view: " << Rendering::CAMERA_FOV << "°\n\n";

    // Create ray tracer
    Schwarzschild metric(Physics::r_s);
    Camera camera(Rendering::CAMERA_DISTANCE,
                  Rendering::CAMERA_FOV * M_PI / 180.0,
                  static_cast<double>(width) / height);
    RayTracer tracer(metric, camera);

    // Render
    std::cout << "Starting ray tracing...\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<Color> image = tracer.render(width, height);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(
        end_time - start_time);

    std::cout << "\nRendering completed in " << duration.count()
              << " seconds\n";

    // Save image
    write_ppm(output_file, image, width, height);

    std::cout << "\nVisualization tips:\n";
    std::cout << "  - Black region: Black hole shadow (event horizon)\n";
    std::cout << "  - Bright ring: Accretion disk with gravitational lensing\n";
    std::cout << "  - Distorted background: Gravitational lensing effects\n";
    std::cout << "\nTo convert to PNG: convert " << output_file
              << " output.png\n";

    return 0;
}
