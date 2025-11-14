#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "schwarzschild.h"
#include <functional>
#include <array>

/**
 * 4th-Order Runge-Kutta (RK4) integrator for geodesic equations
 *
 * Solves the system: dY/dλ = f(λ, Y)
 * where Y = [t, r, θ, φ, v_t, v_r, v_θ, v_φ] is the 8-component state vector
 *
 * RK4 algorithm:
 * k₁ = h·f(λₙ, Yₙ)
 * k₂ = h·f(λₙ + h/2, Yₙ + k₁/2)
 * k₃ = h·f(λₙ + h/2, Yₙ + k₂/2)
 * k₄ = h·f(λₙ + h, Yₙ + k₃)
 * Yₙ₊₁ = Yₙ + (k₁ + 2k₂ + 2k₃ + k₄)/6
 */
class RK4Integrator {
public:
    using State = Schwarzschild::State;
    using DerivativeFunc = std::function<State(double, const State&)>;

private:
    DerivativeFunc derivatives;
    double step_size;

    /**
     * Vector addition for State
     */
    static State add(const State& a, const State& b) {
        State result;
        for (size_t i = 0; i < 8; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    /**
     * Scalar multiplication for State
     */
    static State multiply(const State& a, double scalar) {
        State result;
        for (size_t i = 0; i < 8; ++i) {
            result[i] = a[i] * scalar;
        }
        return result;
    }

public:
    RK4Integrator(DerivativeFunc f, double h = Physics::STEP_SIZE)
        : derivatives(f), step_size(h) {}

    /**
     * Perform one RK4 step
     * @param lambda Current affine parameter value
     * @param y Current state
     * @return New state after one step
     */
    State step(double lambda, const State& y) {
        // k1 = h * f(λ, Y)
        State k1 = multiply(derivatives(lambda, y), step_size);

        // k2 = h * f(λ + h/2, Y + k1/2)
        State k2 = multiply(derivatives(lambda + step_size/2.0,
                                       add(y, multiply(k1, 0.5))),
                           step_size);

        // k3 = h * f(λ + h/2, Y + k2/2)
        State k3 = multiply(derivatives(lambda + step_size/2.0,
                                       add(y, multiply(k2, 0.5))),
                           step_size);

        // k4 = h * f(λ + h, Y + k3)
        State k4 = multiply(derivatives(lambda + step_size,
                                       add(y, k3)),
                           step_size);

        // Yₙ₊₁ = Yₙ + (k₁ + 2k₂ + 2k₃ + k₄)/6
        State result = y;
        for (size_t i = 0; i < 8; ++i) {
            result[i] += (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        }

        return result;
    }

    /**
     * Integrate from lambda0 to lambda1
     * @param lambda0 Starting affine parameter
     * @param lambda1 Ending affine parameter
     * @param y0 Initial state
     * @return Final state
     */
    State integrate(double lambda0, double lambda1, const State& y0) {
        State y = y0;
        double lambda = lambda0;

        int num_steps = static_cast<int>((lambda1 - lambda0) / step_size);

        for (int i = 0; i < num_steps; ++i) {
            y = step(lambda, y);
            lambda += step_size;
        }

        return y;
    }

    void set_step_size(double h) { step_size = h; }
    double get_step_size() const { return step_size; }
};

/**
 * Geodesic integrator specifically for Schwarzschild metric
 */
class GeodesicIntegrator {
private:
    Schwarzschild metric;
    RK4Integrator integrator;

public:
    using State = Schwarzschild::State;

    GeodesicIntegrator(const Schwarzschild& schw, double step_size = Physics::STEP_SIZE)
        : metric(schw),
          integrator([this](double lambda, const State& y) {
              return this->geodesic_derivatives(lambda, y);
          }, step_size) {}

    /**
     * Compute derivatives for geodesic equation
     * Returns dY/dλ where Y = [t, r, θ, φ, v_t, v_r, v_θ, v_φ]
     *
     * First 4 components: velocities (dx^μ/dλ)
     * Last 4 components: accelerations (d²x^μ/dλ²)
     */
    State geodesic_derivatives(double lambda, const State& y) {
        State dydt;

        // Position derivatives = velocities
        dydt[0] = y[4];  // dt/dλ = v_t
        dydt[1] = y[5];  // dr/dλ = v_r
        dydt[2] = y[6];  // dθ/dλ = v_θ
        dydt[3] = y[7];  // dφ/dλ = v_φ

        // Velocity derivatives = accelerations (from geodesic equation)
        auto accel = metric.geodesic_acceleration(y);
        dydt[4] = accel[0];  // d²t/dλ²
        dydt[5] = accel[1];  // d²r/dλ²
        dydt[6] = accel[2];  // d²θ/dλ²
        dydt[7] = accel[3];  // d²φ/dλ²

        return dydt;
    }

    /**
     * Perform one integration step
     */
    State step(double lambda, const State& y) {
        return integrator.step(lambda, y);
    }

    /**
     * Check if geodesic has terminated (captured or escaped)
     */
    bool is_terminated(const State& y) const {
        double r = y[1];
        return metric.is_captured(r) || metric.has_escaped(r);
    }

    const Schwarzschild& get_metric() const { return metric; }
};

#endif // INTEGRATOR_H
