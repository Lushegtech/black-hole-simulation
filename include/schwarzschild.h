#ifndef SCHWARZSCHILD_H
#define SCHWARZSCHILD_H

#include "vec4.h"
#include "mat4.h"
#include "constants.h"
#include <cmath>
#include <array>

/**
 * Schwarzschild metric and geodesic equations
 *
 * Metric: ds² = -(1 - rs/r)dt² + (1 - rs/r)⁻¹dr² + r²(dθ² + sin²θ dφ²)
 *
 * Using conserved quantities (Lagrangian method):
 * - E: Energy (conserved due to time-translation symmetry)
 * - L: Angular momentum (conserved due to rotational symmetry)
 *
 * For motion in equatorial plane (θ = π/2):
 * dt/dλ = E / (1 - 2/r)
 * dφ/dλ = L / r²
 * (dr/dλ)² = E² - (1 - 2/r)(δ + L²/r²)
 *
 * where δ = 0 for photons (null geodesics), δ = 1 for massive particles
 */
class Schwarzschild {
public:
    double r_s;  // Schwarzschild radius

    Schwarzschild(double schwarzschild_radius = Physics::r_s)
        : r_s(schwarzschild_radius) {}

    /**
     * State vector: [t, r, theta, phi, v_t, v_r, v_theta, v_phi]
     * Position: (t, r, θ, φ) in spherical coordinates
     * Velocity: (dt/dλ, dr/dλ, dθ/dλ, dφ/dλ)
     */
    using State = std::array<double, 8>;

    /**
     * Compute metric tensor g_μν at given position
     */
    Mat4 metric(double r, double theta) const {
        Mat4 g;

        double sin_theta = std::sin(theta);
        double f = 1.0 - r_s / r;  // (1 - rs/r)

        // Schwarzschild metric in spherical coordinates
        g(0, 0) = -f;                    // g_tt
        g(1, 1) = 1.0 / f;               // g_rr
        g(2, 2) = r * r;                 // g_θθ
        g(3, 3) = r * r * sin_theta * sin_theta;  // g_φφ

        // Off-diagonal elements are zero (diagonal metric)
        g(0, 1) = g(1, 0) = 0.0;
        g(0, 2) = g(2, 0) = 0.0;
        g(0, 3) = g(3, 0) = 0.0;
        g(1, 2) = g(2, 1) = 0.0;
        g(1, 3) = g(3, 1) = 0.0;
        g(2, 3) = g(3, 2) = 0.0;

        return g;
    }

    /**
     * Compute inverse metric g^μν
     */
    Mat4 inverse_metric(double r, double theta) const {
        Mat4 g_inv;

        double sin_theta = std::sin(theta);
        double f = 1.0 - r_s / r;

        g_inv(0, 0) = -1.0 / f;
        g_inv(1, 1) = f;
        g_inv(2, 2) = 1.0 / (r * r);
        g_inv(3, 3) = 1.0 / (r * r * sin_theta * sin_theta);

        // Off-diagonal elements are zero
        g_inv(0, 1) = g_inv(1, 0) = 0.0;
        g_inv(0, 2) = g_inv(2, 0) = 0.0;
        g_inv(0, 3) = g_inv(3, 0) = 0.0;
        g_inv(1, 2) = g_inv(2, 1) = 0.0;
        g_inv(1, 3) = g_inv(3, 1) = 0.0;
        g_inv(2, 3) = g_inv(3, 2) = 0.0;

        return g_inv;
    }

    /**
     * Compute Christoffel symbols Γ^μ_νσ at given position
     * Only non-zero components are computed
     */
    struct ChristoffelSymbols {
        double Gamma_r_tt, Gamma_t_tr, Gamma_r_rr;
        double Gamma_r_thetatheta, Gamma_r_phiphi;
        double Gamma_theta_rtheta, Gamma_theta_phiphi;
        double Gamma_phi_rphi, Gamma_phi_thetaphi;
    };

    ChristoffelSymbols christoffel(double r, double theta) const {
        ChristoffelSymbols Gamma;

        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double f = 1.0 - r_s / r;

        // Non-zero Christoffel symbols for Schwarzschild metric
        Gamma.Gamma_r_tt = (r_s * (r - r_s)) / (2.0 * r * r * r);
        Gamma.Gamma_t_tr = r_s / (2.0 * r * (r - r_s));
        Gamma.Gamma_r_rr = -r_s / (2.0 * r * (r - r_s));
        Gamma.Gamma_r_thetatheta = -(r - r_s);
        Gamma.Gamma_r_phiphi = -(r - r_s) * sin_theta * sin_theta;
        Gamma.Gamma_theta_rtheta = 1.0 / r;
        Gamma.Gamma_theta_phiphi = -sin_theta * cos_theta;
        Gamma.Gamma_phi_rphi = 1.0 / r;
        Gamma.Gamma_phi_thetaphi = cos_theta / sin_theta;

        return Gamma;
    }

    /**
     * Geodesic equations using Christoffel symbols
     * Returns d²x^μ/dλ² = -Γ^μ_νσ (dx^ν/dλ)(dx^σ/dλ)
     */
    std::array<double, 4> geodesic_acceleration(const State& y) const {
        double t = y[0], r = y[1], theta = y[2], phi = y[3];
        double v_t = y[4], v_r = y[5], v_theta = y[6], v_phi = y[7];

        auto Gamma = christoffel(r, theta);

        std::array<double, 4> accel;

        // d²t/dλ²
        accel[0] = -2.0 * Gamma.Gamma_t_tr * v_t * v_r;

        // d²r/dλ²
        accel[1] = -Gamma.Gamma_r_tt * v_t * v_t
                   - Gamma.Gamma_r_rr * v_r * v_r
                   - Gamma.Gamma_r_thetatheta * v_theta * v_theta
                   - Gamma.Gamma_r_phiphi * v_phi * v_phi;

        // d²θ/dλ²
        accel[2] = -2.0 * Gamma.Gamma_theta_rtheta * v_r * v_theta
                   - Gamma.Gamma_theta_phiphi * v_phi * v_phi;

        // d²φ/dλ²
        accel[3] = -2.0 * Gamma.Gamma_phi_rphi * v_r * v_phi
                   - 2.0 * Gamma.Gamma_phi_thetaphi * v_theta * v_phi;

        return accel;
    }

    /**
     * Compute conserved energy E for a geodesic
     * E = -g_tμ (dx^μ/dλ)
     */
    double conserved_energy(const State& y) const {
        double r = y[1], theta = y[2];
        double v_t = y[4];

        Mat4 g = metric(r, theta);
        return -g(0, 0) * v_t;  // E = -(1 - rs/r) * v_t
    }

    /**
     * Compute conserved angular momentum L for equatorial geodesic
     * L = g_φμ (dx^μ/dλ)
     */
    double conserved_angular_momentum(const State& y) const {
        double r = y[1], theta = y[2];
        double v_phi = y[7];

        Mat4 g = metric(r, theta);
        return g(3, 3) * v_phi;  // L = r² sin²θ * v_φ
    }

    /**
     * Check if ray has been captured by black hole
     */
    bool is_captured(double r) const {
        return r < r_s * 1.01;  // Small buffer for numerical stability
    }

    /**
     * Check if ray has escaped to infinity
     */
    bool has_escaped(double r) const {
        return r > Physics::ESCAPE_RADIUS;
    }
};

#endif // SCHWARZSCHILD_H
