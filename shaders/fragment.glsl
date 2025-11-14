#version 330 core

/**
 * General Relativistic Ray Tracing Fragment Shader
 * Schwarzschild Black Hole Simulation
 *
 * This shader implements the complete GRRT pipeline:
 * 1. Initialize null geodesic for each pixel
 * 2. Integrate using RK4
 * 3. Check termination conditions (captured, escaped, disk intersection)
 * 4. Return final color
 */

out vec4 FragColor;
in vec2 TexCoord;

// Uniforms
uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;  // Camera position (r, theta, phi)
uniform float u_fov;

// Physical constants (in natural units where G = M = c = 1)
const float r_s = 2.0;           // Schwarzschild radius
const float r_horizon = 2.0;     // Event horizon
const float r_photon = 3.0;      // Photon sphere
const float r_ISCO = 6.0;        // Innermost stable circular orbit
const float r_disk_inner = 6.0;  // Accretion disk inner edge
const float r_disk_outer = 20.0; // Accretion disk outer edge
const float escape_radius = 100.0;
const int MAX_STEPS = 500;       // Reduced for real-time performance
const float STEP_SIZE = 0.15;

const float PI = 3.14159265359;
const float EPSILON = 1e-6;

// State vector: [t, r, theta, phi, v_t, v_r, v_theta, v_phi]
struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Add two states
State state_add(State a, State b) {
    State result;
    result.t = a.t + b.t;
    result.r = a.r + b.r;
    result.theta = a.theta + b.theta;
    result.phi = a.phi + b.phi;
    result.v_t = a.v_t + b.v_t;
    result.v_r = a.v_r + b.v_r;
    result.v_theta = a.v_theta + b.v_theta;
    result.v_phi = a.v_phi + b.v_phi;
    return result;
}

// Multiply state by scalar
State state_mul(State a, float s) {
    State result;
    result.t = a.t * s;
    result.r = a.r * s;
    result.theta = a.theta * s;
    result.phi = a.phi * s;
    result.v_t = a.v_t * s;
    result.v_r = a.v_r * s;
    result.v_theta = a.v_theta * s;
    result.v_phi = a.v_phi * s;
    return result;
}

// Compute geodesic derivatives (equations of motion)
State geodesic_derivatives(State y) {
    State dydt;

    float r = y.r;
    float theta = y.theta;
    float v_t = y.v_t;
    float v_r = y.v_r;
    float v_theta = y.v_theta;
    float v_phi = y.v_phi;

    // Avoid singularities
    if (r < r_s + EPSILON) r = r_s + EPSILON;
    if (theta < EPSILON) theta = EPSILON;
    if (theta > PI - EPSILON) theta = PI - EPSILON;

    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    if (abs(sin_theta) < EPSILON) sin_theta = EPSILON;

    // Position derivatives = velocities
    dydt.t = v_t;
    dydt.r = v_r;
    dydt.theta = v_theta;
    dydt.phi = v_phi;

    // Christoffel symbols
    float Gamma_r_tt = (r_s * (r - r_s)) / (2.0 * r * r * r);
    float Gamma_t_tr = r_s / (2.0 * r * (r - r_s));
    float Gamma_r_rr = -r_s / (2.0 * r * (r - r_s));
    float Gamma_r_thetatheta = -(r - r_s);
    float Gamma_r_phiphi = -(r - r_s) * sin_theta * sin_theta;
    float Gamma_theta_rtheta = 1.0 / r;
    float Gamma_theta_phiphi = -sin_theta * cos_theta;
    float Gamma_phi_rphi = 1.0 / r;
    float Gamma_phi_thetaphi = cos_theta / sin_theta;

    // Velocity derivatives = accelerations (geodesic equation)
    dydt.v_t = -2.0 * Gamma_t_tr * v_t * v_r;

    dydt.v_r = -Gamma_r_tt * v_t * v_t
               - Gamma_r_rr * v_r * v_r
               - Gamma_r_thetatheta * v_theta * v_theta
               - Gamma_r_phiphi * v_phi * v_phi;

    dydt.v_theta = -2.0 * Gamma_theta_rtheta * v_r * v_theta
                   - Gamma_theta_phiphi * v_phi * v_phi;

    dydt.v_phi = -2.0 * Gamma_phi_rphi * v_r * v_phi
                 - 2.0 * Gamma_phi_thetaphi * v_theta * v_phi;

    return dydt;
}

// 4th-order Runge-Kutta step
State rk4_step(State y, float h) {
    State k1 = state_mul(geodesic_derivatives(y), h);

    State y2 = state_add(y, state_mul(k1, 0.5));
    State k2 = state_mul(geodesic_derivatives(y2), h);

    State y3 = state_add(y, state_mul(k2, 0.5));
    State k3 = state_mul(geodesic_derivatives(y3), h);

    State y4 = state_add(y, k3);
    State k4 = state_mul(geodesic_derivatives(y4), h);

    // Combine: y_next = y + (k1 + 2*k2 + 2*k3 + k4) / 6
    State result = y;
    State increment = state_add(
        state_add(k1, state_mul(k2, 2.0)),
        state_add(state_mul(k3, 2.0), k4)
    );
    increment = state_mul(increment, 1.0 / 6.0);
    result = state_add(result, increment);

    return result;
}

// Initialize ray for given pixel
State initialize_ray(vec2 pixel_coord) {
    State state;

    // Camera position
    float r_cam = u_camera_pos.x;
    float theta_cam = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    // Initial position
    state.t = 0.0;
    state.r = r_cam;
    state.theta = theta_cam;
    state.phi = phi_cam;

    // Ray direction based on pixel
    float tan_fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    float offset_x = pixel_coord.x * tan_fov * aspect;
    float offset_y = pixel_coord.y * tan_fov;

    // Initial 4-velocity (photon)
    float v_r = -1.0;  // Toward black hole
    float v_theta = offset_y * 0.5;
    float v_phi = offset_x * 0.5;

    // Solve for v_t using null condition: g_μν v^μ v^ν = 0
    float f = 1.0 - r_s / r_cam;
    float sin_theta = sin(theta_cam);

    float spatial_term = (v_r * v_r) / f
                       + r_cam * r_cam * (v_theta * v_theta
                       + sin_theta * sin_theta * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial_term / f));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

// Background color (checkerboard)
vec3 get_background_color(float theta, float phi) {
    int checker_x = int(floor(phi * 4.0 / PI));
    int checker_y = int(floor(theta * 4.0 / PI));

    bool is_white = ((checker_x + checker_y) & 1) == 0;

    if (is_white) {
        return vec3(0.8, 0.8, 0.9);
    } else {
        return vec3(0.1, 0.1, 0.3);
    }
}

// Accretion disk color
vec3 get_disk_color(float r) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature profile: T ∝ r^(-3/4)
    float temp_ratio = pow(r_disk_inner / r, 0.75);
    float t_norm = min(1.0, temp_ratio);

    // Blackbody color approximation
    vec3 color = vec3(1.0, 0.5 + 0.5 * t_norm, t_norm);

    // Brightness falloff
    float brightness = pow(r_disk_inner / r, 2.0);
    brightness = min(1.0, brightness);

    return color * brightness;
}

// Check disk intersection
bool check_disk_intersection(State s_old, State s_new, out float r_intersect) {
    float pi_2 = PI / 2.0;

    // Check if crossed equatorial plane
    bool crossed = (s_old.theta - pi_2) * (s_new.theta - pi_2) < 0.0;

    if (crossed) {
        // Linear interpolation
        float t = (pi_2 - s_old.theta) / (s_new.theta - s_old.theta);
        r_intersect = s_old.r + t * (s_new.r - s_old.r);

        return r_intersect >= r_disk_inner && r_intersect <= r_disk_outer;
    }

    return false;
}

// Main ray tracing function
vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;
        float theta = state.theta;
        float phi = state.phi;

        // Condition A: Captured by black hole
        if (r < r_horizon * 1.01) {
            return vec3(0.0);  // Black
        }

        // Condition B: Escaped to infinity
        if (r > escape_radius) {
            return get_background_color(theta, phi);
        }

        // Condition C: Hit accretion disk
        float r_intersect;
        if (check_disk_intersection(state_old, state, r_intersect)) {
            return get_disk_color(r_intersect);
        }
    }

    // Max steps reached
    return vec3(0.0);
}

void main() {
    // Normalize pixel coordinates to [-1, 1]
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;  // Flip Y

    // Trace ray
    vec3 color = trace_ray(pixel_coord);

    // Gamma correction
    color = pow(color, vec3(1.0 / 2.2));

    FragColor = vec4(color, 1.0);
}
