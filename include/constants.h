#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

namespace Physics {
// Natural units where G = M = c = 1
// This simplifies the Schwarzschild radius to rs = 2
constexpr double G = 1.0; // Gravitational constant
constexpr double c = 1.0; // Speed of light
constexpr double M = 1.0; // Black hole mass

// Schwarzschild radius: rs = 2GM/c²
// In natural units: rs = 2
constexpr double r_s = 2.0 * G * M / (c * c);

// Event horizon radius (same as Schwarzschild radius for non-rotating BH)
constexpr double r_horizon = r_s;

// Photon sphere: where photons can orbit (unstable)
// For Schwarzschild: r_photon = 1.5 * r_s = 3.0
constexpr double r_photon = 1.5 * r_s;

// Innermost Stable Circular Orbit (ISCO)
// For Schwarzschild: r_ISCO = 3 * r_s = 6.0
constexpr double r_ISCO = 3.0 * r_s;

// Accretion disk parameters
constexpr double r_disk_inner = r_ISCO; // Inner edge at ISCO
constexpr double r_disk_outer = 20.0;   // Outer edge (arbitrary)

// Integration parameters
constexpr double STEP_SIZE = 0.1;       // Fixed step size for RK4
constexpr int MAX_STEPS = 10000;        // Maximum integration steps
constexpr double ESCAPE_RADIUS = 100.0; // Ray escapes if r > this

// Numerical tolerances
constexpr double EPSILON = 1e-10; // Small number for comparisons
} // namespace Physics

namespace Rendering {
// Default image resolution
constexpr int WIDTH = 800;
constexpr int HEIGHT = 600;

// Camera parameters
constexpr double CAMERA_DISTANCE = 20.0; // Distance from black hole
constexpr double CAMERA_FOV = 60.0;      // Field of view in degrees

// Background
constexpr double BACKGROUND_R = 0.0;
constexpr double BACKGROUND_G = 0.0;
constexpr double BACKGROUND_B = 0.0;

// Accretion disk colors (temperature-based)
constexpr double DISK_TEMP_BASE = 10000.0; // Base temperature in Kelvin
} // namespace Rendering

#endif // CONSTANTS_H
