#ifndef ADAPTIVE_INTEGRATOR_H
#define ADAPTIVE_INTEGRATOR_H

#include "kerr.h"

#include <array>
#include <cmath>
#include <functional>

/**
 * Runge-Kutta-Fehlberg (RKF45) Adaptive Step-Size Integrator
 *
 * This is the CORRECT integrator for stiff geodesic equations near a black hole.
 * A fixed-step integrator (like RK4) is inadequate because:
 *
 * - Far from BH: Spacetime is nearly flat → can use large steps (fast)
 * - Near horizon: Curvature changes rapidly → requires tiny steps (accurate)
 *
 * RKF45 automatically adapts the step size based on local truncation error,
 * providing both:
 * 1. High performance (60+ FPS) by taking large steps when safe
 * 2. Physical accuracy by taking tiny steps near the event horizon
 *
 * Algorithm:
 * - Computes TWO solutions per step: 4th-order and 5th-order
 * - Error estimate: Δ = |Y₅ - Y₄|
 * - If Δ ≤ TOL: Accept step, increase h
 * - If Δ > TOL: Reject step, decrease h, retry
 *
 * This is the method used by professional GRRT codes (GYOTO, etc.)
 */
class RKF45Integrator {
  public:
    using State = Kerr::State; // 8-component [x^μ, p_μ]
    using DerivativeFunc = std::function<State(const State&)>;

  private:
    DerivativeFunc derivatives;
    double tolerance;      // Error tolerance for adaptive stepping
    double h_min, h_max;   // Min/max allowed step sizes
    double safety_factor;  // Safety margin for step size adjustment

    // Butcher tableau coefficients for RKF45
    // https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method
    static constexpr double a2 = 1.0 / 4.0;
    static constexpr double a3 = 3.0 / 8.0;
    static constexpr double a4 = 12.0 / 13.0;
    static constexpr double a5 = 1.0;
    static constexpr double a6 = 1.0 / 2.0;

    static constexpr double b21 = 1.0 / 4.0;
    static constexpr double b31 = 3.0 / 32.0;
    static constexpr double b32 = 9.0 / 32.0;
    static constexpr double b41 = 1932.0 / 2197.0;
    static constexpr double b42 = -7200.0 / 2197.0;
    static constexpr double b43 = 7296.0 / 2197.0;
    static constexpr double b51 = 439.0 / 216.0;
    static constexpr double b52 = -8.0;
    static constexpr double b53 = 3680.0 / 513.0;
    static constexpr double b54 = -845.0 / 4104.0;
    static constexpr double b61 = -8.0 / 27.0;
    static constexpr double b62 = 2.0;
    static constexpr double b63 = -3544.0 / 2565.0;
    static constexpr double b64 = 1859.0 / 4104.0;
    static constexpr double b65 = -11.0 / 40.0;

    // 4th-order solution coefficients
    static constexpr double c1 = 25.0 / 216.0;
    static constexpr double c2 = 0.0;
    static constexpr double c3 = 1408.0 / 2565.0;
    static constexpr double c4 = 2197.0 / 4104.0;
    static constexpr double c5 = -1.0 / 5.0;
    static constexpr double c6 = 0.0;

    // 5th-order solution coefficients
    static constexpr double d1 = 16.0 / 135.0;
    static constexpr double d2 = 0.0;
    static constexpr double d3 = 6656.0 / 12825.0;
    static constexpr double d4 = 28561.0 / 56430.0;
    static constexpr double d5 = -9.0 / 50.0;
    static constexpr double d6 = 2.0 / 55.0;

    /**
     * Vector addition
     */
    static State add(const State& a, const State& b) {
        State result;
        for (size_t i = 0; i < 8; ++i)
            result[i] = a[i] + b[i];
        return result;
    }

    /**
     * Scalar multiplication
     */
    static State multiply(const State& a, double scalar) {
        State result;
        for (size_t i = 0; i < 8; ++i)
            result[i] = a[i] * scalar;
        return result;
    }

    /**
     * L2 norm of state vector (for error estimation)
     */
    static double norm(const State& a) {
        double sum = 0.0;
        for (size_t i = 0; i < 8; ++i)
            sum += a[i] * a[i];
        return std::sqrt(sum);
    }

  public:
    RKF45Integrator(DerivativeFunc f, double tol = 1e-6, double h_min_val = 1e-6,
                    double h_max_val = 1.0, double safety = 0.9)
        : derivatives(f), tolerance(tol), h_min(h_min_val), h_max(h_max_val),
          safety_factor(safety) {}

    /**
     * Perform one adaptive RKF45 step
     *
     * @param y Current state (will be updated)
     * @param h Current step size (will be updated)
     * @return true if step was accepted, false if rejected (caller should retry)
     */
    bool step(State& y, double& h) {
        // Compute the six stages k1, k2, ..., k6
        State k1 = derivatives(y);

        State k2 = derivatives(add(y, multiply(k1, h * b21)));

        State k3 = derivatives(
            add(add(y, multiply(k1, h * b31)), multiply(k2, h * b32)));

        State k4 = derivatives(
            add(add(add(y, multiply(k1, h * b41)), multiply(k2, h * b42)),
                multiply(k3, h * b43)));

        State k5 = derivatives(add(
            add(add(add(y, multiply(k1, h * b51)), multiply(k2, h * b52)),
                multiply(k3, h * b53)),
            multiply(k4, h * b54)));

        State k6 = derivatives(
            add(add(add(add(add(y, multiply(k1, h * b61)), multiply(k2, h * b62)),
                        multiply(k3, h * b63)),
                    multiply(k4, h * b64)),
                multiply(k5, h * b65)));

        // Compute 4th-order solution
        State y4 = y;
        for (size_t i = 0; i < 8; ++i) {
            y4[i] += h * (c1 * k1[i] + c2 * k2[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] +
                          c6 * k6[i]);
        }

        // Compute 5th-order solution
        State y5 = y;
        for (size_t i = 0; i < 8; ++i) {
            y5[i] += h * (d1 * k1[i] + d2 * k2[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] +
                          d6 * k6[i]);
        }

        // Estimate error
        State error;
        for (size_t i = 0; i < 8; ++i)
            error[i] = y5[i] - y4[i];

        double error_norm = norm(error);

        // Compute optimal step size for next iteration
        double h_optimal;
        if (error_norm < 1e-12) {
            // Error is negligible - increase step size significantly
            h_optimal = h_max;
        } else {
            // Standard step size control formula
            // h_new = safety × h × (TOL / error)^(1/5)
            // Exponent is 1/(order+1) = 1/5 for 4th-order method
            h_optimal = safety_factor * h * std::pow(tolerance / error_norm, 0.20);
        }

        // Clamp to allowed range and apply growth limiters
        h_optimal = std::max(h_min, std::min(h_max, h_optimal));
        h_optimal = std::min(h_optimal, 5.0 * h);  // Don't grow too fast
        h_optimal = std::max(h_optimal, 0.1 * h);  // Don't shrink too fast

        // Check if step is accepted
        if (error_norm <= tolerance) {
            // SUCCESS: Accept the 5th-order solution (more accurate)
            y = y5;
            h = h_optimal;
            return true;
        } else {
            // FAILURE: Reject step and retry with smaller h
            h = h_optimal;
            return false;
        }
    }

    /**
     * Integrate until a termination condition is met
     *
     * @param y0 Initial state
     * @param h0 Initial step size
     * @param max_steps Maximum integration steps
     * @param should_terminate Function that returns true when integration should stop
     * @return Final state
     */
    State integrate(const State& y0, double h0, int max_steps,
                    std::function<bool(const State&)> should_terminate) {
        State y = y0;
        double h = h0;

        for (int i = 0; i < max_steps; ++i) {
            // Check termination condition
            if (should_terminate(y)) {
                break;
            }

            // Attempt step (may be rejected and retried)
            int retry_count = 0;
            while (!step(y, h) && retry_count < 10) {
                retry_count++;
            }

            if (retry_count >= 10) {
                // Failed to find acceptable step - likely near singularity
                break;
            }
        }

        return y;
    }

    // Getters/setters
    void set_tolerance(double tol) { tolerance = tol; }
    double get_tolerance() const { return tolerance; }

    void set_step_limits(double h_min_val, double h_max_val) {
        h_min = h_min_val;
        h_max = h_max_val;
    }
};

/**
 * Geodesic integrator for Kerr metric with adaptive stepping
 */
class KerrGeodesicIntegrator {
  private:
    Kerr metric;
    RKF45Integrator integrator;

  public:
    using State = Kerr::State;

    KerrGeodesicIntegrator(const Kerr& kerr_metric, double tolerance = 1e-6)
        : metric(kerr_metric),
          integrator(
              [this](const State& y) { return this->metric.hamiltonian_derivatives(y); },
              tolerance) {}

    /**
     * Perform one adaptive integration step
     */
    bool step(State& y, double& h) { return integrator.step(y, h); }

    /**
     * Check if geodesic has terminated (captured or escaped)
     */
    bool is_terminated(const State& y) const {
        double r = y[1];
        return metric.is_captured(r) || metric.has_escaped(r);
    }

    const Kerr& get_metric() const { return metric; }

    void set_tolerance(double tol) { integrator.set_tolerance(tol); }
};

#endif // ADAPTIVE_INTEGRATOR_H
