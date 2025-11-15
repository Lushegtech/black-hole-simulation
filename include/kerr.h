#ifndef KERR_H
#define KERR_H

#include "constants.h"
#include "mat4.h"
#include "vec4.h"

#include <array>
#include <cmath>

/**
 * Kerr metric and geodesic equations for rotating black holes
 *
 * The Kerr metric describes a rotating black hole in Boyer-Lindquist coordinates.
 * This is the CORRECT metric for Sagittarius A*, which is a rotating supermassive
 * black hole (Schwarzschild is WRONG as it assumes no rotation).
 *
 * Metric parameters:
 * - M: Black hole mass (in geometric units where G=c=1)
 * - a: Spin parameter (0 ≤ a < M, dimensionless a* = a/M)
 *
 * Key functions of Boyer-Lindquist coordinates (r, θ):
 * - ρ²(r,θ) = r² + a²cos²θ
 * - Δ(r) = r² - 2Mr + a²
 * - Σ(r,θ) = (r² + a²)² - a²Δsin²θ
 *
 * Event horizon: r₊ = M + √(M² - a²)
 * Ergosphere: rₑ(θ) = M + √(M² - a²cos²θ)
 */
class Kerr {
  public:
    double M;  // Black hole mass
    double a;  // Spin parameter (a = a* × M, where a* ∈ [0,1])

    Kerr(double mass = Physics::M, double spin = 0.0) : M(mass), a(spin * mass) {}

    /**
     * State vector for Hamiltonian formulation: [t, r, θ, φ, p_t, p_r, p_θ, p_φ]
     * Position: (t, r, θ, φ) in Boyer-Lindquist coordinates
     * Momentum (covariant): (p_t, p_r, p_θ, p_φ)
     *
     * This 8-component formulation avoids sign ambiguities in the 4-component
     * formulation based on constants of motion (E, L, Q).
     */
    using State = std::array<double, 8>;

    /**
     * Compute ρ²(r, θ) = r² + a²cos²θ
     */
    double rho2(double r, double theta) const {
        double cos_theta = std::cos(theta);
        return r * r + a * a * cos_theta * cos_theta;
    }

    /**
     * Compute Δ(r) = r² - 2Mr + a²
     */
    double Delta(double r) const { return r * r - 2.0 * M * r + a * a; }

    /**
     * Compute Σ(r, θ) = (r² + a²)² - a²Δsin²θ
     */
    double Sigma(double r, double theta) const {
        double sin_theta = std::sin(theta);
        double r2_plus_a2 = r * r + a * a;
        return r2_plus_a2 * r2_plus_a2 - a * a * Delta(r) * sin_theta * sin_theta;
    }

    /**
     * Compute outer event horizon radius: r₊ = M + √(M² - a²)
     */
    double event_horizon() const { return M + std::sqrt(M * M - a * a); }

    /**
     * Compute ergosphere radius at angle θ: rₑ(θ) = M + √(M² - a²cos²θ)
     */
    double ergosphere(double theta) const {
        double cos_theta = std::cos(theta);
        return M + std::sqrt(M * M - a * a * cos_theta * cos_theta);
    }

    /**
     * Compute ISCO (Innermost Stable Circular Orbit) radius
     * This depends on black hole spin and is critical for accretion disk physics.
     *
     * For prograde orbits (most astrophysically relevant):
     * - Schwarzschild (a=0): r_isco = 6M
     * - Maximal Kerr (a=M): r_isco ≈ 1M
     *
     * Formula from Bardeen, Press, & Teukolsky (1972)
     */
    double isco_radius() const {
        double a_star = a / M; // Dimensionless spin

        // Z1 and Z2 functions
        double Z1 = 1.0 + std::pow(1.0 - a_star * a_star, 1.0 / 3.0) *
                              (std::pow(1.0 + a_star, 1.0 / 3.0) +
                               std::pow(1.0 - a_star, 1.0 / 3.0));
        double Z2 = std::sqrt(3.0 * a_star * a_star + Z1 * Z1);

        // ISCO radius (prograde orbit)
        double r_isco = M * (3.0 + Z2 - std::sqrt((3.0 - Z1) * (3.0 + Z1 + 2.0 * Z2)));

        return r_isco;
    }

    /**
     * Compute inverse metric g^μν in Boyer-Lindquist coordinates
     * Required for Hamiltonian formulation: dx^μ/dλ = g^μν p_ν
     */
    void inverse_metric_components(double r, double theta, double g_inv[4][4]) const {
        double rho_sq = rho2(r, theta);
        double Delta_val = Delta(r);
        double Sigma_val = Sigma(r, theta);
        double sin_theta = std::sin(theta);
        double sin2_theta = sin_theta * sin_theta;

        // Initialize to zero
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                g_inv[i][j] = 0.0;

        // Diagonal components
        g_inv[0][0] = -(r * r + a * a + 2.0 * M * r * a * a * sin2_theta / rho_sq) /
                      Delta_val; // g^tt
        g_inv[1][1] = Delta_val / rho_sq;                                 // g^rr
        g_inv[2][2] = 1.0 / rho_sq;                                       // g^θθ
        g_inv[3][3] = (Delta_val - a * a * sin2_theta) / (Delta_val * rho_sq * sin2_theta); // g^φφ

        // Off-diagonal (only t-φ coupling in Kerr)
        g_inv[0][3] = g_inv[3][0] = -2.0 * M * r * a / (Delta_val * rho_sq); // g^tφ
    }

    /**
     * Compute partial derivatives of inverse metric: ∂g^μν/∂x^σ
     * Required for Hamiltonian equations: dp_μ/dλ = -(1/2)(∂g^νσ/∂x^μ)p_ν p_σ
     *
     * Returns derivatives w.r.t. coordinates: [∂/∂t, ∂/∂r, ∂/∂θ, ∂/∂φ]
     * Note: ∂/∂t = ∂/∂φ = 0 (metric is stationary and axisymmetric)
     */
    void metric_derivatives(double r, double theta, double dg_dr[4][4],
                            double dg_dtheta[4][4]) const {
        double sin_theta = std::sin(theta);
        double cos_theta = std::cos(theta);
        double sin2_theta = sin_theta * sin_theta;
        double cos2_theta = cos_theta * cos_theta;

        double rho_sq = rho2(r, theta);
        double Delta_val = Delta(r);

        // Derivatives of basis functions
        double drho2_dr = 2.0 * r;
        double drho2_dtheta = -2.0 * a * a * sin_theta * cos_theta;

        double dDelta_dr = 2.0 * r - 2.0 * M;

        // Initialize to zero
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                dg_dr[i][j] = 0.0;
                dg_dtheta[i][j] = 0.0;
            }
        }

        // Compute derivatives (complex but standard formulas from Kerr metric)
        // These are algebraically intensive but necessary for accurate geodesics

        // ∂g^rr/∂r
        dg_dr[1][1] =
            (dDelta_dr * rho_sq - Delta_val * drho2_dr) / (rho_sq * rho_sq);

        // ∂g^rr/∂θ
        dg_dtheta[1][1] = -Delta_val * drho2_dtheta / (rho_sq * rho_sq);

        // ∂g^θθ/∂r
        dg_dr[2][2] = -drho2_dr / (rho_sq * rho_sq);

        // ∂g^θθ/∂θ
        dg_dtheta[2][2] = -drho2_dtheta / (rho_sq * rho_sq);

        // ∂g^tt/∂r (complex)
        double A = r * r + a * a + 2.0 * M * r * a * a * sin2_theta / rho_sq;
        double dA_dr = 2.0 * r + 2.0 * M * a * a * sin2_theta / rho_sq -
                       2.0 * M * r * a * a * sin2_theta * drho2_dr / (rho_sq * rho_sq);
        dg_dr[0][0] = -(dA_dr * Delta_val - A * dDelta_dr) / (Delta_val * Delta_val);

        // ∂g^tt/∂θ
        double dA_dtheta = 2.0 * M * r * a * a * 2.0 * sin_theta * cos_theta / rho_sq -
                           2.0 * M * r * a * a * sin2_theta * drho2_dtheta /
                               (rho_sq * rho_sq);
        dg_dtheta[0][0] = -dA_dtheta / Delta_val;

        // ∂g^φφ/∂r
        double B = Delta_val - a * a * sin2_theta;
        double dB_dr = dDelta_dr;
        dg_dr[3][3] = (dB_dr * Delta_val * rho_sq * sin2_theta -
                       B * (dDelta_dr * rho_sq * sin2_theta +
                            Delta_val * drho2_dr * sin2_theta)) /
                      (Delta_val * Delta_val * rho_sq * rho_sq * sin2_theta * sin2_theta);

        // ∂g^φφ/∂θ
        double dB_dtheta = -2.0 * a * a * sin_theta * cos_theta;
        dg_dtheta[3][3] =
            (dB_dtheta * Delta_val * rho_sq * sin2_theta +
             B * (-Delta_val * drho2_dtheta * sin2_theta -
                  Delta_val * rho_sq * 2.0 * sin_theta * cos_theta)) /
            (Delta_val * Delta_val * rho_sq * rho_sq * sin2_theta * sin2_theta);

        // ∂g^tφ/∂r (off-diagonal coupling)
        double C = -2.0 * M * r * a;
        double dC_dr = -2.0 * M * a;
        dg_dr[0][3] = dg_dr[3][0] =
            (dC_dr * Delta_val * rho_sq - C * (dDelta_dr * rho_sq + Delta_val * drho2_dr)) /
            (Delta_val * Delta_val * rho_sq * rho_sq);

        // ∂g^tφ/∂θ
        dg_dtheta[0][3] = dg_dtheta[3][0] =
            C * Delta_val * drho2_dtheta / (Delta_val * Delta_val * rho_sq * rho_sq);
    }

    /**
     * Hamiltonian formulation: compute dY/dλ where Y = [x^μ, p_μ]
     *
     * Hamilton's equations:
     * dx^μ/dλ = ∂H/∂p_μ = g^μν p_ν
     * dp_μ/dλ = -∂H/∂x^μ = -(1/2)(∂g^νσ/∂x^μ)p_ν p_σ
     *
     * This is the core physics engine for the simulation.
     */
    State hamiltonian_derivatives(const State& y) const {
        State dydt;

        // Extract state
        double r = y[1];
        double theta = y[2];
        double p_t = y[4], p_r = y[5], p_theta = y[6], p_phi = y[7];

        // Get inverse metric g^μν
        double g_inv[4][4];
        inverse_metric_components(r, theta, g_inv);

        // Position derivatives: dx^μ/dλ = g^μν p_ν
        dydt[0] = g_inv[0][0] * p_t + g_inv[0][3] * p_phi; // dt/dλ
        dydt[1] = g_inv[1][1] * p_r;                       // dr/dλ
        dydt[2] = g_inv[2][2] * p_theta;                   // dθ/dλ
        dydt[3] = g_inv[3][0] * p_t + g_inv[3][3] * p_phi; // dφ/dλ

        // Get metric derivatives
        double dg_dr[4][4], dg_dtheta[4][4];
        metric_derivatives(r, theta, dg_dr, dg_dtheta);

        // Momentum derivatives: dp_μ/dλ = -(1/2)(∂g^νσ/∂x^μ)p_ν p_σ
        // dp_t/dλ = 0 (metric is time-independent)
        dydt[4] = 0.0;

        // dp_r/dλ = -(1/2)(∂g^νσ/∂r)p_ν p_σ
        dydt[5] = 0.0;
        for (int nu = 0; nu < 4; ++nu) {
            for (int sigma = 0; sigma < 4; ++sigma) {
                double p_nu = (nu == 0 ? p_t : (nu == 1 ? p_r : (nu == 2 ? p_theta : p_phi)));
                double p_sigma =
                    (sigma == 0 ? p_t : (sigma == 1 ? p_r : (sigma == 2 ? p_theta : p_phi)));
                dydt[5] -= 0.5 * dg_dr[nu][sigma] * p_nu * p_sigma;
            }
        }

        // dp_θ/dλ = -(1/2)(∂g^νσ/∂θ)p_ν p_σ
        dydt[6] = 0.0;
        for (int nu = 0; nu < 4; ++nu) {
            for (int sigma = 0; sigma < 4; ++sigma) {
                double p_nu = (nu == 0 ? p_t : (nu == 1 ? p_r : (nu == 2 ? p_theta : p_phi)));
                double p_sigma =
                    (sigma == 0 ? p_t : (sigma == 1 ? p_r : (sigma == 2 ? p_theta : p_phi)));
                dydt[6] -= 0.5 * dg_dtheta[nu][sigma] * p_nu * p_sigma;
            }
        }

        // dp_φ/dλ = 0 (metric is φ-independent / axisymmetric)
        dydt[7] = 0.0;

        return dydt;
    }

    /**
     * Check if ray has been captured by black hole (crossed event horizon)
     */
    bool is_captured(double r) const { return r < event_horizon() * 1.01; }

    /**
     * Check if ray has escaped to infinity
     */
    bool has_escaped(double r) const { return r > Physics::ESCAPE_RADIUS; }

    /**
     * Get disk 4-velocity for circular orbit at radius r (Novikov-Thorne model)
     * This is needed for calculating the relativistic redshift factor g.
     *
     * Returns normalized 4-velocity u^μ = (u^t, 0, 0, u^φ)
     * for a particle in a circular, equatorial, prograde Keplerian orbit.
     */
    std::array<double, 4> disk_four_velocity(double r) const {
        std::array<double, 4> u;

        double r2 = r * r;
        double Delta_val = Delta(r);

        // Angular velocity Ω for Keplerian orbit
        double Omega = M / (r * std::sqrt(r * M));

        // u^t and u^φ from normalization condition g_μν u^μ u^ν = -1
        double A = -(1.0 - 2.0 * M / r);
        double B = -4.0 * M * a / r;
        double C = (r2 + a * a + 2.0 * M * a * a / r);

        double u_phi = Omega;
        double u_t = std::sqrt(-(A + 2.0 * B * Omega + C * Omega * Omega)) / (A + B * Omega);

        u[0] = u_t;
        u[1] = 0.0; // Circular orbit: no radial motion
        u[2] = 0.0; // Equatorial: no polar motion
        u[3] = u_phi;

        return u;
    }
};

#endif // KERR_H
