#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "schwarzschild.h"
#include "integrator.h"
#include "vec4.h"
#include <cmath>
#include <vector>

/**
 * Color structure for RGB values
 */
struct Color {
    double r, g, b;

    Color(double red = 0.0, double green = 0.0, double blue = 0.0)
        : r(red), g(green), b(blue) {}

    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }

    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }

    // Clamp to [0, 1]
    void clamp() {
        r = std::max(0.0, std::min(1.0, r));
        g = std::max(0.0, std::min(1.0, g));
        b = std::max(0.0, std::min(1.0, b));
    }
};

/**
 * Camera for ray generation
 */
struct Camera {
    Vec4 position;      // Camera position in spacetime
    Vec4 direction;     // Looking direction
    Vec4 up;            // Up vector
    double fov;         // Field of view in radians
    double aspect;      // Aspect ratio (width/height)

    Camera(double distance = Rendering::CAMERA_DISTANCE,
           double field_of_view = Rendering::CAMERA_FOV * M_PI / 180.0,
           double aspect_ratio = static_cast<double>(Rendering::WIDTH) / Rendering::HEIGHT)
        : position(0.0, distance, M_PI/2.0, 0.0),  // (t, r, θ, φ)
          direction(0.0, -1.0, 0.0, 0.0),           // Looking toward BH
          up(0.0, 0.0, 1.0, 0.0),                   // Up in θ direction
          fov(field_of_view),
          aspect(aspect_ratio) {}
};

/**
 * General Relativistic Ray Tracer
 * Implements backward ray tracing through curved spacetime
 */
class RayTracer {
private:
    Schwarzschild metric;
    GeodesicIntegrator integrator;
    Camera camera;

    /**
     * Initialize a null geodesic (photon path) for a given pixel
     * @param pixel_x Pixel x coordinate (normalized to [-1, 1])
     * @param pixel_y Pixel y coordinate (normalized to [-1, 1])
     * @return Initial state vector for the geodesic
     */
    Schwarzschild::State initialize_ray(double pixel_x, double pixel_y) const {
        Schwarzschild::State state;

        // Camera position in spherical coordinates
        double r_cam = camera.position.x();      // radial distance
        double theta_cam = camera.position.y();  // polar angle
        double phi_cam = camera.position.z();    // azimuthal angle

        // Initial position
        state[0] = 0.0;         // t
        state[1] = r_cam;       // r
        state[2] = theta_cam;   // θ
        state[3] = phi_cam;     // φ

        // Calculate initial velocity (4-velocity for photon)
        // For a photon, we need a null vector: g_μν v^μ v^ν = 0

        // Simple approach: shoot ray toward black hole with slight offset
        // based on pixel coordinates

        double tan_fov = std::tan(camera.fov / 2.0);

        // Ray direction in camera space
        double offset_x = pixel_x * tan_fov * camera.aspect;
        double offset_y = pixel_y * tan_fov;

        // Initial 4-velocity components
        // For equatorial plane motion, start with:
        // - radial velocity pointing inward
        // - angular velocity based on pixel offset

        double v_r = -1.0;  // Moving toward black hole
        double v_theta = offset_y * 0.5;
        double v_phi = offset_x * 0.5;

        // Solve for v_t using null condition: g_μν v^μ v^ν = 0
        // For Schwarzschild in spherical coords:
        // -(1 - rs/r)(v_t)² + (1 - rs/r)⁻¹(v_r)² + r²(v_θ)² + r²sin²θ(v_φ)² = 0

        double f = 1.0 - metric.r_s / r_cam;
        double sin_theta = std::sin(theta_cam);

        double spatial_term = (v_r * v_r) / f
                            + r_cam * r_cam * (v_theta * v_theta
                            + sin_theta * sin_theta * v_phi * v_phi);

        double v_t = std::sqrt(spatial_term / f);

        state[4] = v_t;
        state[5] = v_r;
        state[6] = v_theta;
        state[7] = v_phi;

        return state;
    }

    /**
     * Get background color (simple checkerboard pattern for now)
     */
    Color get_background_color(double theta, double phi) const {
        // Simple checkerboard pattern
        int checker_x = static_cast<int>(std::floor(phi * 4.0 / M_PI));
        int checker_y = static_cast<int>(std::floor(theta * 4.0 / M_PI));

        bool is_white = (checker_x + checker_y) % 2 == 0;

        if (is_white) {
            return Color(0.8, 0.8, 0.8);
        } else {
            return Color(0.2, 0.2, 0.4);
        }
    }

    /**
     * Get accretion disk color at given radius
     * Uses temperature profile T(r) ∝ r^(-3/4)
     */
    Color get_disk_color(double r) const {
        if (r < Physics::r_disk_inner || r > Physics::r_disk_outer) {
            return Color(0.0, 0.0, 0.0);
        }

        // Temperature decreases with radius: T ∝ r^(-3/4)
        double temp_ratio = std::pow(Physics::r_disk_inner / r, 0.75);
        double temperature = Rendering::DISK_TEMP_BASE * temp_ratio;

        // Simple blackbody color approximation
        // Hot = blue-white, cooler = red-orange
        double t_norm = std::min(1.0, temperature / 15000.0);

        Color color;
        color.r = 1.0;
        color.g = 0.5 + 0.5 * t_norm;
        color.b = t_norm;

        // Brightness falls off with radius
        double brightness = std::pow(Physics::r_disk_inner / r, 2.0);
        brightness = std::min(1.0, brightness);

        return color * brightness;
    }

    /**
     * Check if ray intersects the accretion disk
     * Disk is in equatorial plane (θ = π/2)
     */
    bool check_disk_intersection(const Schwarzschild::State& state_old,
                                  const Schwarzschild::State& state_new,
                                  double& intersection_r) const {
        double theta_old = state_old[2];
        double theta_new = state_new[2];
        double r_old = state_old[1];
        double r_new = state_new[1];

        // Check if ray crossed the equatorial plane
        double pi_2 = M_PI / 2.0;
        bool crossed = (theta_old - pi_2) * (theta_new - pi_2) < 0;

        if (crossed) {
            // Linear interpolation to find intersection radius
            double t = (pi_2 - theta_old) / (theta_new - theta_old);
            intersection_r = r_old + t * (r_new - r_old);

            // Check if within disk bounds
            return intersection_r >= Physics::r_disk_inner &&
                   intersection_r <= Physics::r_disk_outer;
        }

        return false;
    }

public:
    RayTracer(const Schwarzschild& schw = Schwarzschild(),
              const Camera& cam = Camera())
        : metric(schw), integrator(schw), camera(cam) {}

    /**
     * Trace a single ray and return its color
     * Implements the GRRT pipeline from Section V
     */
    Color trace_ray(double pixel_x, double pixel_y) {
        // Initialize ray
        auto state = initialize_ray(pixel_x, pixel_y);
        auto state_old = state;

        double lambda = 0.0;

        // Integration loop
        for (int step = 0; step < Physics::MAX_STEPS; ++step) {
            // Take one integration step
            state_old = state;
            state = integrator.step(lambda, state);
            lambda += Physics::STEP_SIZE;

            double r = state[1];
            double theta = state[2];
            double phi = state[3];

            // Condition A: Ray captured by black hole
            if (metric.is_captured(r)) {
                return Color(0.0, 0.0, 0.0);  // Black
            }

            // Condition B: Ray escaped to infinity
            if (metric.has_escaped(r)) {
                // Sample background
                return get_background_color(theta, phi);
            }

            // Condition C: Ray intersects accretion disk
            double intersection_r;
            if (check_disk_intersection(state_old, state, intersection_r)) {
                return get_disk_color(intersection_r);
            }
        }

        // Maximum steps reached
        return Color(0.0, 0.0, 0.0);
    }

    /**
     * Render entire image
     */
    std::vector<Color> render(int width, int height) {
        std::vector<Color> image(width * height);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // Normalize pixel coordinates to [-1, 1]
                double px = 2.0 * (x + 0.5) / width - 1.0;
                double py = 1.0 - 2.0 * (y + 0.5) / height;  // Flip Y

                Color color = trace_ray(px, py);
                image[y * width + x] = color;
            }

            // Progress indicator
            if (y % 10 == 0) {
                double progress = 100.0 * y / height;
                std::cout << "\rRendering: " << static_cast<int>(progress)
                         << "% complete" << std::flush;
            }
        }

        std::cout << "\rRendering: 100% complete" << std::endl;

        return image;
    }

    void set_camera(const Camera& cam) { camera = cam; }
};

#endif // RAYTRACER_H
