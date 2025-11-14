#version 330 core

/**
 * CINEMATIC Black Hole Fragment Shader
 * Inspired by Interstellar (2014) - DNGR renderer
 *
 * Features:
 * - Procedural starfield background
 * - Hot plasma accretion disk
 * - Relativistic Doppler shift and beaming
 * - Realistic colors and lighting
 */

out vec4 FragColor;
in vec2 TexCoord;

// Uniforms
uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;  // Camera position (r, theta, phi)
uniform float u_fov;

// Physical constants
const float r_s = 2.0;
const float r_horizon = 2.0;
const float r_photon = 3.0;
const float r_ISCO = 6.0;
const float r_disk_inner = 6.0;
const float r_disk_outer = 20.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 500;
const float STEP_SIZE = 0.15;
const float PI = 3.14159265359;
const float EPSILON = 1e-6;

// State vector for geodesic integration
struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Hash function for procedural generation
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453123);
}

// Procedural starfield
vec3 get_starfield(vec3 dir) {
    // Convert direction to spherical coordinates
    float theta = acos(dir.y);
    float phi = atan(dir.z, dir.x);

    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    // Multiple layers of stars at different scales
    vec3 color = vec3(0.0);

    // Bright stars
    for (int i = 0; i < 100; i++) {
        vec2 starPos = vec2(hash(vec2(float(i), 0.0)), hash(vec2(float(i), 1.0)));
        float dist = length(uv - starPos);
        if (dist < 0.002) {
            float brightness = 1.0 - (dist / 0.002);
            brightness = pow(brightness, 3.0);
            // Star color variation
            vec3 starColor = mix(vec3(0.8, 0.9, 1.0), vec3(1.0, 0.95, 0.8), hash(vec2(float(i), 2.0)));
            color += starColor * brightness;
        }
    }

    // Dim stars
    for (int i = 100; i < 500; i++) {
        vec2 starPos = vec2(hash(vec2(float(i), 0.0)), hash(vec2(float(i), 1.0)));
        float dist = length(uv - starPos);
        if (dist < 0.001) {
            float brightness = 0.3 * (1.0 - (dist / 0.001));
            color += vec3(brightness);
        }
    }

    // Milky Way-like glow
    float galaxy = pow(abs(sin(uv.x * PI * 3.0)), 2.0) * pow(abs(sin(uv.y * PI)), 4.0);
    color += vec3(0.15, 0.1, 0.2) * galaxy * 0.3;

    // Background space color
    color += vec3(0.01, 0.01, 0.015);

    return clamp(color, 0.0, 1.0);
}

// Blackbody radiation color (Planck's law approximation)
vec3 blackbody_color(float temperature) {
    // Simplified blackbody color for temperatures 1000K - 20000K
    float t = clamp(temperature / 10000.0, 0.1, 2.0);

    vec3 color;
    if (t < 0.5) {
        // Red to orange
        color = mix(vec3(1.0, 0.2, 0.0), vec3(1.0, 0.6, 0.2), t * 2.0);
    } else {
        // Orange to yellow-white
        color = mix(vec3(1.0, 0.6, 0.2), vec3(1.0, 0.95, 0.8), (t - 0.5) * 2.0);
    }

    return color;
}

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

State geodesic_derivatives(State y) {
    State dydt;
    float r = y.r;
    float theta = y.theta;
    float v_t = y.v_t;
    float v_r = y.v_r;
    float v_theta = y.v_theta;
    float v_phi = y.v_phi;

    if (r < r_s + EPSILON) r = r_s + EPSILON;
    if (theta < EPSILON) theta = EPSILON;
    if (theta > PI - EPSILON) theta = PI - EPSILON;

    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    if (abs(sin_theta) < EPSILON) sin_theta = EPSILON;

    dydt.t = v_t;
    dydt.r = v_r;
    dydt.theta = v_theta;
    dydt.phi = v_phi;

    float Gamma_r_tt = (r_s * (r - r_s)) / (2.0 * r * r * r);
    float Gamma_t_tr = r_s / (2.0 * r * (r - r_s));
    float Gamma_r_rr = -r_s / (2.0 * r * (r - r_s));
    float Gamma_r_thetatheta = -(r - r_s);
    float Gamma_r_phiphi = -(r - r_s) * sin_theta * sin_theta;
    float Gamma_theta_rtheta = 1.0 / r;
    float Gamma_theta_phiphi = -sin_theta * cos_theta;
    float Gamma_phi_rphi = 1.0 / r;
    float Gamma_phi_thetaphi = cos_theta / sin_theta;

    dydt.v_t = -2.0 * Gamma_t_tr * v_t * v_r;
    dydt.v_r = -Gamma_r_tt * v_t * v_t - Gamma_r_rr * v_r * v_r
               - Gamma_r_thetatheta * v_theta * v_theta - Gamma_r_phiphi * v_phi * v_phi;
    dydt.v_theta = -2.0 * Gamma_theta_rtheta * v_r * v_theta - Gamma_theta_phiphi * v_phi * v_phi;
    dydt.v_phi = -2.0 * Gamma_phi_rphi * v_r * v_phi - 2.0 * Gamma_phi_thetaphi * v_theta * v_phi;

    return dydt;
}

State rk4_step(State y, float h) {
    State k1 = state_mul(geodesic_derivatives(y), h);
    State y2 = state_add(y, state_mul(k1, 0.5));
    State k2 = state_mul(geodesic_derivatives(y2), h);
    State y3 = state_add(y, state_mul(k2, 0.5));
    State k3 = state_mul(geodesic_derivatives(y3), h);
    State y4 = state_add(y, k3);
    State k4 = state_mul(geodesic_derivatives(y4), h);

    State increment = state_add(state_add(k1, state_mul(k2, 2.0)),
                                state_add(state_mul(k3, 2.0), k4));
    increment = state_mul(increment, 1.0 / 6.0);
    return state_add(y, increment);
}

State initialize_ray(vec2 pixel_coord) {
    State state;
    float r_cam = u_camera_pos.x;
    float theta_cam = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    state.t = 0.0;
    state.r = r_cam;
    state.theta = theta_cam;
    state.phi = phi_cam;

    float tan_fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;
    float offset_x = pixel_coord.x * tan_fov * aspect;
    float offset_y = pixel_coord.y * tan_fov;

    float v_r = -1.0;
    float v_theta = offset_y * 0.5;
    float v_phi = offset_x * 0.5;

    float f = 1.0 - r_s / r_cam;
    float sin_theta = sin(theta_cam);
    float spatial_term = (v_r * v_r) / f + r_cam * r_cam * (v_theta * v_theta
                       + sin_theta * sin_theta * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial_term / f));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

vec3 get_disk_color_cinematic(float r, float v_phi_disk, float v_phi_photon) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature profile: T ∝ r^(-3/4)
    float temp_ratio = pow(r_disk_inner / r, 0.75);
    float base_temp = 8000.0;  // Kelvin
    float temperature = base_temp * temp_ratio;

    // Relativistic Doppler shift
    // Disk rotates, so we compare angular velocities
    float doppler_factor = 1.0 + (v_phi_disk - v_phi_photon) * 0.3;
    float observed_temp = temperature * doppler_factor;

    // Blackbody color based on temperature
    vec3 color = blackbody_color(observed_temp);

    // Relativistic beaming: brightness ∝ (doppler_factor)^3
    float beaming = pow(max(doppler_factor, 0.1), 3.0);

    // Brightness falloff with radius
    float brightness = pow(r_disk_inner / r, 2.5) * beaming;
    brightness = clamp(brightness, 0.0, 3.0);

    return color * brightness;
}

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

vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;
        float theta = state.theta;
        float phi = state.phi;

        // Captured by black hole
        if (r < r_horizon * 1.01) {
            return vec3(0.0);
        }

        // Escaped to infinity - show starfield
        if (r > escape_radius) {
            vec3 dir = vec3(
                sin(theta) * cos(phi),
                cos(theta),
                sin(theta) * sin(phi)
            );
            return get_starfield(dir);
        }

        // Hit accretion disk
        float r_intersect, phi_intersect;
        if (check_disk_intersection(state_old, state, r_intersect, phi_intersect)) {
            // Calculate disk rotation velocity (Keplerian orbit)
            float v_phi_disk = sqrt(1.0 / r_intersect) * 0.5;
            return get_disk_color_cinematic(r_intersect, v_phi_disk, state.v_phi);
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    vec3 color = trace_ray(pixel_coord);

    // Tone mapping and color grading
    color = color / (color + vec3(1.0));  // Reinhard tone mapping
    color = pow(color, vec3(0.8));  // Slight gamma adjustment for space

    // Add slight blue tint to space
    color = mix(color, color * vec3(0.9, 0.95, 1.0), 0.1);

    FragColor = vec4(color, 1.0);
}
