#version 330 core

/**
 * SAGITTARIUS A* - SIMPLIFIED & WORKING
 * Clear black hole shadow with visible accretion disk
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float M = 1.0;
const float a = 0.7;
const float r_s = 2.0 * M;
const float r_horizon = 1.7;
const float r_photon = 2.6;
const float r_isco = 2.3;
const float r_disk_inner = 2.5;
const float r_disk_outer = 15.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 200;
const float STEP_SIZE = 0.15;
const float PI = 3.14159265359;

struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Simple hash for stars
float hash(vec2 p) {
    return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453);
}

// Simple starfield
vec3 get_starfield(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 color = vec3(0.01, 0.01, 0.02);  // Dark background

    // Stars
    for (int i = 0; i < 50; i++) {
        vec2 starPos = vec2(hash(vec2(float(i), 0.0)), hash(vec2(float(i), 1.0)));
        float dist = length(uv - starPos);
        if (dist < 0.003) {
            color += vec3(1.0, 0.95, 0.9) * (1.0 - dist / 0.003);
        }
    }

    return color;
}

// Bright orange disk
vec3 get_disk_color(float r, float phi) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Brightness falls off with radius
    float brightness = pow(r_disk_inner / r, 1.2);
    brightness = clamp(brightness, 0.0, 3.0);

    // Rotation pattern
    float pattern = 0.8 + 0.2 * sin(phi * 10.0 + u_time * 0.5 / sqrt(r));

    // Orange/yellow EHT colors
    vec3 color_inner = vec3(1.0, 0.85, 0.5);   // Bright yellow
    vec3 color_outer = vec3(1.0, 0.4, 0.1);    // Deep orange

    float mix_factor = (r - r_disk_inner) / (r_disk_outer - r_disk_inner);
    vec3 color = mix(color_inner, color_outer, mix_factor);

    // Doppler boost (simple approximation)
    float doppler = 1.0 + 0.3 * sin(phi);
    brightness *= pow(max(0.5, doppler), 3.0);

    return color * brightness * pattern * 2.0;
}

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

// Simplified Kerr-like geodesics
State geodesic_derivatives(State s) {
    State dydt;
    float r = s.r;
    float theta = s.theta;

    dydt.t = s.v_t;
    dydt.r = s.v_r;
    dydt.theta = s.v_theta;
    dydt.phi = s.v_phi;

    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    // Simplified gravity (works well visually)
    float grav_pull = M / max(r * r, 0.1);

    dydt.v_t = 0.0;
    dydt.v_r = -grav_pull * s.v_t * s.v_t + (r - r_s) * s.v_theta * s.v_theta / max(r * r * r, 0.1);
    dydt.v_r += (r - r_s) * sin_theta * sin_theta * s.v_phi * s.v_phi / max(r * r * r, 0.1);
    dydt.v_theta = -2.0 * s.v_r * s.v_theta / max(r, 0.1);
    dydt.v_theta += sin_theta * cos_theta * s.v_phi * s.v_phi;
    dydt.v_phi = -2.0 * (s.v_r * s.v_phi / max(r, 0.1) + s.v_theta * s.v_phi * cos_theta / max(sin_theta, 0.01));

    return dydt;
}

State rk4_step(State y, float h) {
    State k1 = state_mul(geodesic_derivatives(y), h);
    State k2 = state_mul(geodesic_derivatives(state_add(y, state_mul(k1, 0.5))), h);
    State k3 = state_mul(geodesic_derivatives(state_add(y, state_mul(k2, 0.5))), h);
    State k4 = state_mul(geodesic_derivatives(state_add(y, k3)), h);

    return state_add(y, state_mul(state_add(state_add(k1, state_mul(k2, 2.0)), state_add(state_mul(k3, 2.0), k4)), 1.0 / 6.0));
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

    state.v_r = -1.0;
    state.v_theta = pixel_coord.y * tan_fov * 0.5;
    state.v_phi = pixel_coord.x * tan_fov * aspect * 0.5;

    // Normalize
    float speed = sqrt(state.v_r * state.v_r + state.v_theta * state.v_theta + state.v_phi * state.v_phi);
    state.v_t = speed;

    return state;
}

bool check_disk_hit(State s_old, State s_new, out float r_hit, out float phi_hit) {
    float pi_2 = PI / 2.0;
    bool crossed = (s_old.theta - pi_2) * (s_new.theta - pi_2) < 0.0;

    if (crossed) {
        float t = (pi_2 - s_old.theta) / (s_new.theta - s_old.theta);
        r_hit = s_old.r + t * (s_new.r - s_old.r);
        phi_hit = s_old.phi + t * (s_new.phi - s_old.phi);
        return r_hit >= r_disk_inner && r_hit <= r_disk_outer;
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

        // Hit event horizon - BLACK
        if (r < r_horizon * 1.05) {
            return vec3(0.0);
        }

        // Escaped to infinity - STARS
        if (r > escape_radius) {
            vec3 ray_dir = normalize(vec3(
                sin(state.theta) * cos(state.phi),
                cos(state.theta),
                sin(state.theta) * sin(state.phi)
            ));
            return get_starfield(ray_dir);
        }

        // Hit disk - BRIGHT ORANGE
        float r_hit, phi_hit;
        if (check_disk_hit(state_old, state, r_hit, phi_hit)) {
            return get_disk_color(r_hit, phi_hit);
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    vec3 color = trace_ray(pixel_coord);

    // Simple tone mapping
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(0.9));

    FragColor = vec4(color, 1.0);
}
