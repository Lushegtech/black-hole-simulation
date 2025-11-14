#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "integrator.h"
#include "schwarzschild.h"
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
    Vec4 position;  // Camera position in spacetime
    Vec4 direction; // Looking direction
    Vec4 up;        // Up vector
    double fov;     // Field of view in radians
    double aspect;  // Aspect ratio (width/height)

    Camera(double distance = Rendering::CAMERA_DISTANCE,
           double field_of_view = Rendering::CAMERA_FOV * M_PI / 180.0,
           double aspect_ratio = static_cast<double>(Rendering::WIDTH) /
                                 Rendering::HEIGHT)
        : position(0.0, distance, M_PI / 2.0, 0.0), // (t, r, θ, φ)
          direction(0.0, -1.0, 0.0, 0.0),           // Looking toward BH
          up(0.0, 0.0, 1.0, 0.0),                   // Up in θ direction
          fov(field_of_view), aspect(aspect_ratio) {}
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
        double r_cam = camera.position.x();     // radial distance
        double theta_cam = camera.position.y(); // polar angle
        double phi_cam = camera.position.z();   // azimuthal angle

        // Initial position
        state[0] = 0.0;       // t
        state[1] = r_cam;     // r
        state[2] = theta_cam; // θ
        state[3] = phi_cam;   // φ

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

        double v_r = -1.0; // Moving toward black hole
        double v_theta = offset_y * 0.5;
        double v_phi = offset_x * 0.5;

        // Solve for v_t using null condition: g_μν v^μ v^ν = 0
        // For Schwarzschild in spherical coords:
        // -(1 - rs/r)(v_t)² + (1 - rs/r)⁻¹(v_r)² + r²(v_θ)² + r²sin²θ(v_φ)² = 0

        double f = 1.0 - metric.r_s / r_cam;
        double sin_theta = std::sin(theta_cam);

        double spatial_term =
            (v_r * v_r) / f +
            r_cam * r_cam * (v_theta * v_theta + sin_theta * sin_theta * v_phi * v_phi);

        double v_t = std::sqrt(spatial_term / f);

        state[4] = v_t;
        state[5] = v_r;
        state[6] = v_theta;
        state[7] = v_phi;

        return state;
    }

    /**
     * Hash function for procedural star generation
     */
    double hash(double x, double y) const {
        double val = std::sin(x * 127.1 + y * 311.7) * 43758.5453123;
        return val - std::floor(val);
    }

    /**
     * Get cinematic starfield background
     */
    Color get_background_color(double theta, double phi) const {
        // Normalize to [0,1] UV coordinates
        double u = phi / (2.0 * M_PI) + 0.5;
        double v = theta / M_PI;

        Color color(0.0, 0.0, 0.0);

        // Bright stars (100 stars)
        for (int i = 0; i < 100; ++i) {
            double star_u = hash(i, 0.0);
            double star_v = hash(i, 1.0);
            double dist = std::sqrt((u - star_u) * (u - star_u) + (v - star_v) * (v - star_v));

            if (dist < 0.002) {
                double brightness = 1.0 - (dist / 0.002);
                brightness = std::pow(brightness, 3.0);

                // Star color variation (blue to yellow-white)
                double star_hue = hash(i, 2.0);
                Color star_color = Color(0.8 + 0.2 * star_hue, 0.9 + 0.05 * star_hue,
                                        1.0 - 0.2 * (1.0 - star_hue));
                color = color + star_color * brightness;
            }
        }

        // Dim background stars (100 more)
        for (int i = 100; i < 200; ++i) {
            double star_u = hash(i, 0.0);
            double star_v = hash(i, 1.0);
            double dist = std::sqrt((u - star_u) * (u - star_u) + (v - star_v) * (v - star_v));

            if (dist < 0.001) {
                double brightness = 0.3 * (1.0 - (dist / 0.001));
                color = color + Color(brightness, brightness, brightness);
            }
        }

        // Milky Way-like glow
        double galaxy = std::pow(std::abs(std::sin(u * M_PI * 3.0)), 2.0) *
                       std::pow(std::abs(std::sin(v * M_PI)), 4.0);
        color = color + Color(0.15, 0.1, 0.2) * (galaxy * 0.3);

        // Base space color
        color = color + Color(0.01, 0.01, 0.015);

        color.clamp();
        return color;
    }

    /**
     * Blackbody color approximation
     */
    Color blackbody_color(double temperature) const {
        double t = std::min(2.0, std::max(0.1, temperature / 10000.0));

        Color color;
        if (t < 0.5) {
            // Red to orange (cooler)
            double mix_factor = t * 2.0;
            color.r = 1.0;
            color.g = 0.2 + 0.4 * mix_factor;
            color.b = 0.0 + 0.2 * mix_factor;
        } else {
            // Orange to yellow-white (hotter)
            double mix_factor = (t - 0.5) * 2.0;
            color.r = 1.0;
            color.g = 0.6 + 0.35 * mix_factor;
            color.b = 0.2 + 0.6 * mix_factor;
        }

        return color;
    }

    /**
     * Get cinematic accretion disk color with relativistic effects
     * Includes Doppler shift and relativistic beaming
     */
    Color get_disk_color(double r, double v_phi_photon = 0.0) const {
        if (r < Physics::r_disk_inner || r > Physics::r_disk_outer) {
            return Color(0.0, 0.0, 0.0);
        }

        // Temperature profile: T ∝ r^(-3/4)
        double temp_ratio = std::pow(Physics::r_disk_inner / r, 0.75);
        double base_temp = 8000.0; // Base temperature in Kelvin
        double temperature = base_temp * temp_ratio;

        // Keplerian orbital velocity for disk
        double v_phi_disk = std::sqrt(1.0 / r) * 0.5;

        // Relativistic Doppler shift
        double doppler_factor = 1.0 + (v_phi_disk - v_phi_photon) * 0.3;
        double observed_temp = temperature * doppler_factor;

        // Get blackbody color
        Color color = blackbody_color(observed_temp);

        // Relativistic beaming: brightness ∝ (doppler_factor)³
        double beaming = std::pow(std::max(0.1, doppler_factor), 3.0);

        // Brightness falloff with radius
        double brightness = std::pow(Physics::r_disk_inner / r, 2.5) * beaming;
        brightness = std::min(3.0, brightness);

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
    RayTracer(const Schwarzschild& schw = Schwarzschild(), const Camera& cam = Camera())
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
                return Color(0.0, 0.0, 0.0); // Black
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
                double py = 1.0 - 2.0 * (y + 0.5) / height; // Flip Y

                Color color = trace_ray(px, py);
                image[y * width + x] = color;
            }

            // Progress indicator
            if (y % 10 == 0) {
                double progress = 100.0 * y / height;
                std::cout << "\rRendering: " << static_cast<int>(progress) << "% complete"
                          << std::flush;
            }
        }

        std::cout << "\rRendering: 100% complete" << std::endl;

        return image;
    }

    void set_camera(const Camera& cam) {
        camera = cam;
    }
};

#endif // RAYTRACER_H
