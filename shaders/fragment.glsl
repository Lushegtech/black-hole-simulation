#version 330 core

/**
 * Real Black Hole Simulation - Schwarzschild Geodesic Ray Tracing
 * Actual General Relativity physics with optimizations for GTX 1050
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float r_s = 2.0;
const float r_horizon = 2.0;
const float r_disk_inner = 6.0;
const float r_disk_outer = 20.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 300;
const float STEP_SIZE = 0.2;
const float PI = 3.14159265359;

struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Hash for stars
float hash12(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

// Simple starfield
vec3 get_starfield(vec3 dir) {
    float theta = acos(dir.y);
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.0);

    // Stars
    for (int i = 0; i < 80; i++) {
        vec2 starPos = vec2(hash12(vec2(i, 0)), hash12(vec2(i, 1)));
        float dist = length(uv - starPos);

        if (dist < 0.003) {
            float brightness = 1.0 - (dist / 0.003);
            brightness = pow(brightness, 2.0);
            color += vec3(1.0, 0.95, 0.9) * brightness * 0.8;
        }
    }

    return color;
}

// Accretion disk color with temperature gradient
vec3 get_disk_color(float r, float phi) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature: hotter closer to black hole
    float temp_ratio = pow(r_disk_inner / r, 0.75);

    // Rotation effect
    float rotation_phase = phi + u_time * 0.3 / r;
    float brightness = 1.0 + 0.3 * sin(rotation_phase * 8.0);

    // Color based on temperature
    vec3 color;
    if (temp_ratio > 1.5) {
        color = vec3(1.0, 0.95, 0.8);  // White hot
    } else if (temp_ratio > 1.0) {
        color = vec3(1.0, 0.8, 0.4);   // Yellow-orange
    } else {
        color = vec3(1.0, 0.5, 0.2);   // Orange-red
    }

    // Brightness falloff
    float disk_brightness = pow(r_disk_inner / r, 2.0) * brightness;
    disk_brightness = clamp(disk_brightness, 0.0, 3.0);

    return color * disk_brightness;
}

// State operations
State state_add(State a, State b) {
    return State(
        a.t + b.t, a.r + b.r, a.theta + b.theta, a.phi + b.phi,
        a.v_t + b.v_t, a.v_r + b.v_r, a.v_theta + b.v_theta, a.v_phi + b.v_phi
    );
}

State state_mul(State s, float k) {
    return State(
        s.t * k, s.r * k, s.theta * k, s.phi * k,
        s.v_t * k, s.v_r * k, s.v_theta * k, s.v_phi * k
    );
}

// Geodesic equation - Schwarzschild metric
State geodesic_derivatives(State s) {
    State dydt;

    float r = s.r;
    float theta = s.theta;
    float v_t = s.v_t;
    float v_r = s.v_r;
    float v_theta = s.v_theta;
    float v_phi = s.v_phi;

    // Position derivatives
    dydt.t = v_t;
    dydt.r = v_r;
    dydt.theta = v_theta;
    dydt.phi = v_phi;

    // Christoffel symbols and accelerations
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    float Gamma_r_tt = r_s * (r - r_s) / (2.0 * r * r * r);
    float Gamma_t_tr = r_s / (2.0 * r * (r - r_s));
    float Gamma_r_rr = -r_s / (2.0 * r * (r - r_s));
    float Gamma_r_thetatheta = -(r - r_s);
    float Gamma_r_phiphi = -(r - r_s) * sin_theta * sin_theta;
    float Gamma_theta_rtheta = 1.0 / r;
    float Gamma_theta_phiphi = -sin_theta * cos_theta;
    float Gamma_phi_rphi = 1.0 / r;
    float Gamma_phi_thetaphi = cos_theta / (sin_theta + 0.0001);

    dydt.v_t = -2.0 * Gamma_t_tr * v_t * v_r;
    dydt.v_r = -Gamma_r_tt * v_t * v_t - Gamma_r_rr * v_r * v_r
               - Gamma_r_thetatheta * v_theta * v_theta
               - Gamma_r_phiphi * v_phi * v_phi;
    dydt.v_theta = -2.0 * Gamma_theta_rtheta * v_r * v_theta
                   - Gamma_theta_phiphi * v_phi * v_phi;
    dydt.v_phi = -2.0 * Gamma_phi_rphi * v_r * v_phi
                 - 2.0 * Gamma_phi_thetaphi * v_theta * v_phi;

    return dydt;
}

// RK4 integration
State rk4_step(State y, float h) {
    State k1 = state_mul(geodesic_derivatives(y), h);
    State y2 = state_add(y, state_mul(k1, 0.5));
    State k2 = state_mul(geodesic_derivatives(y2), h);
    State y3 = state_add(y, state_mul(k2, 0.5));
    State k3 = state_mul(geodesic_derivatives(y3), h);
    State y4 = state_add(y, k3);
    State k4 = state_mul(geodesic_derivatives(y4), h);

    State increment = state_add(
        state_add(k1, state_mul(k2, 2.0)),
        state_add(state_mul(k3, 2.0), k4)
    );
    increment = state_mul(increment, 1.0 / 6.0);

    return state_add(y, increment);
}

// Initialize ray
State initialize_ray(vec2 pixel_coord) {
    State state;

    float r_cam = u_camera_pos.x;
    float theta_cam = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    state.t = 0.0;
    state.r = r_cam;
    state.theta = theta_cam;
    state.phi = phi_cam;

    // Ray direction
    float tan_fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;
    float offset_x = pixel_coord.x * tan_fov * aspect;
    float offset_y = pixel_coord.y * tan_fov;

    float v_r = -1.0;
    float v_theta = offset_y * 0.5;
    float v_phi = offset_x * 0.5;

    // Normalize to null geodesic
    float f = 1.0 - r_s / r_cam;
    float sin_theta = sin(theta_cam);
    float spatial_term = (v_r * v_r) / f + r_cam * r_cam *
                        (v_theta * v_theta + sin_theta * sin_theta * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial_term / f));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

// Check disk intersection
bool check_disk_intersection(State s_old, State s_new, out float r_intersect, out float phi_intersect) {
    float pi_2 = PI / 2.0;
    bool crossed = (s_old.theta - pi_2) * (s_new.theta - pi_2) < 0.0;

    if (crossed) {
        float t = (pi_2 - s_old.theta) / (s_new.theta - s_old.theta);
        r_intersect = s_old.r + t * (s_new.r - s_old.r);
        phi_intersect = s_old.phi + t * (s_new.phi - s_old.phi);
        return r_intersect >= r_disk_inner && r_intersect <= r_disk_outer;
    }

    return false;
}

// Main ray tracing
vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;

        // Captured by black hole
        if (r < r_horizon * 1.01) {
            return vec3(0.0);
        }

        // Escaped to infinity
        if (r > escape_radius) {
            vec3 ray_dir = vec3(
                sin(state.theta) * cos(state.phi),
                cos(state.theta),
                sin(state.theta) * sin(state.phi)
            );
            return get_starfield(ray_dir);
        }

        // Check disk intersection
        float r_intersect, phi_intersect;
        if (check_disk_intersection(state_old, state, r_intersect, phi_intersect)) {
            return get_disk_color(r_intersect, phi_intersect);
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    vec3 color = trace_ray(pixel_coord);

    // Tone mapping
    color = color / (color + vec3(1.0));

    // Gamma correction
    color = pow(color, vec3(0.85));

    // Vignette
    float dist = length(pixel_coord);
    float vignette = 1.0 - smoothstep(0.7, 1.5, dist);
    color *= vignette * 0.5 + 0.5;

    FragColor = vec4(color, 1.0);
}
