#version 330 core

/**
 * ULTIMATE CINEMATIC Black Hole Shader
 * Full animated simulation with particles, bloom, and turbulence
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float r_s = 2.0;
const float r_horizon = 2.0;
const float r_ISCO = 6.0;
const float r_disk_inner = 6.0;
const float r_disk_outer = 20.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 200;  // Reduced from 500 for Intel HD Graphics
const float STEP_SIZE = 0.25;  // Increased for fewer steps
const float PI = 3.14159265359;
const float EPSILON = 1e-6;

// Particle system constants
const int NUM_PARTICLES = 20;  // Reduced from 50 for performance
const float PARTICLE_SIZE = 0.3;

struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Improved hash functions for better randomness
float hash11(float p) {
    p = fract(p * 0.1031);
    p *= p + 33.33;
    p *= p + p;
    return fract(p);
}

float hash12(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

vec2 hash22(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

// Noise function for turbulence
float noise(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    f = f * f * (3.0 - 2.0 * f);

    float a = hash12(i);
    float b = hash12(i + vec2(1.0, 0.0));
    float c = hash12(i + vec2(0.0, 1.0));
    float d = hash12(i + vec2(1.0, 1.0));

    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

// Fractal Brownian Motion for realistic turbulence
float fbm(vec2 p) {
    float value = 0.0;
    float amplitude = 0.5;
    float frequency = 1.0;

    // Reduced from 5 to 3 octaves for performance on Intel HD Graphics
    for (int i = 0; i < 3; i++) {
        value += amplitude * noise(p * frequency);
        frequency *= 2.0;
        amplitude *= 0.5;
    }

    return value;
}

// Enhanced starfield with twinkling
vec3 get_starfield(vec3 dir) {
    float theta = acos(dir.y);
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.0);

    // Large bright stars with twinkling (reduced for performance)
    for (int i = 0; i < 40; i++) {
        vec2 starPos = hash22(vec2(float(i), 0.0));
        float dist = length(uv - starPos);

        if (dist < 0.003) {
            float brightness = 1.0 - (dist / 0.003);
            brightness = pow(brightness, 3.0);

            // Twinkling effect
            float twinkle = 0.7 + 0.3 * sin(u_time * 2.0 + float(i));
            brightness *= twinkle;

            // Star color
            float hue = hash11(float(i) * 7.123);
            vec3 starColor = mix(vec3(0.7, 0.8, 1.0), vec3(1.0, 0.95, 0.7), hue);
            color += starColor * brightness;
        }
    }

    // Medium stars (reduced for performance)
    for (int i = 40; i < 100; i++) {
        vec2 starPos = hash22(vec2(float(i), 1.0));
        float dist = length(uv - starPos);

        if (dist < 0.0015) {
            float brightness = 0.5 * (1.0 - (dist / 0.0015));
            color += vec3(brightness);
        }
    }

    // Dust and nebula
    float nebula = fbm(uv * 5.0 + u_time * 0.02);
    color += vec3(0.05, 0.02, 0.08) * nebula * 0.5;

    // Milky Way
    float galaxy = pow(abs(sin(uv.x * PI * 3.0)), 2.0) * pow(abs(sin(uv.y * PI)), 4.0);
    color += vec3(0.2, 0.15, 0.25) * galaxy * 0.4;

    // Deep space
    color += vec3(0.01, 0.01, 0.02);

    return clamp(color, 0.0, 1.0);
}

// Blackbody color with HDR
vec3 blackbody_color(float temperature) {
    float t = clamp(temperature / 10000.0, 0.1, 3.0);

    vec3 color;
    if (t < 0.5) {
        color = mix(vec3(1.0, 0.15, 0.0), vec3(1.0, 0.5, 0.1), t * 2.0);
    } else if (t < 1.0) {
        color = mix(vec3(1.0, 0.5, 0.1), vec3(1.0, 0.9, 0.6), (t - 0.5) * 2.0);
    } else {
        color = mix(vec3(1.0, 0.9, 0.6), vec3(0.9, 0.95, 1.0), (t - 1.0) / 2.0);
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
    float v_t = y.v_t, v_r = y.v_r, v_theta = y.v_theta, v_phi = y.v_phi;

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

// Animated accretion disk with turbulence
vec3 get_disk_color_animated(float r, float phi, float v_phi_photon) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Rotating disk - angular position changes with time
    float rotation_speed = 0.3;
    float rotated_phi = phi + u_time * rotation_speed / r;

    // Temperature with turbulence
    float temp_ratio = pow(r_disk_inner / r, 0.75);
    float base_temp = 8000.0;

    // Add turbulent hot spots
    vec2 turb_coord = vec2(rotated_phi * 5.0, r * 0.5);
    float turbulence = fbm(turb_coord + u_time * 0.5);
    turbulence = pow(turbulence, 2.0);

    // Hot spots create temperature variations
    float temperature = base_temp * temp_ratio * (1.0 + turbulence * 0.5);

    // Keplerian velocity
    float v_phi_disk = sqrt(1.0 / r) * 0.5;

    // Doppler shift + rotation
    float doppler_factor = 1.0 + (v_phi_disk * cos(rotated_phi) - v_phi_photon) * 0.4;
    float observed_temp = temperature * doppler_factor;

    vec3 color = blackbody_color(observed_temp);

    // Relativistic beaming
    float beaming = pow(max(0.15, doppler_factor), 3.5);

    // Brightness with turbulent variations
    float brightness = pow(r_disk_inner / r, 2.5) * beaming;
    brightness *= (1.0 + turbulence * 0.3);
    brightness = clamp(brightness, 0.0, 5.0);

    return color * brightness;
}

// Render accretion particles
vec3 get_particle_contribution(vec3 ray_dir, float ray_r) {
    vec3 particle_color = vec3(0.0);

    for (int i = 0; i < NUM_PARTICLES; i++) {
        float seed = float(i) * 13.7;

        // Particle orbit parameters (initialized with hash)
        float p_r = r_disk_inner + hash11(seed) * (r_disk_outer - r_disk_inner);
        float p_phase = hash11(seed + 1.0) * 2.0 * PI;

        // Orbital motion
        float orbit_speed = 1.0 / sqrt(p_r) * 0.5;
        float current_phi = p_phase + u_time * orbit_speed;

        // Spiraling inward slowly
        p_r -= u_time * 0.02 * hash11(seed + 2.0);
        if (p_r < r_disk_inner) continue;

        // 3D position
        vec3 p_pos = vec3(
            p_r * sin(PI/2.0) * cos(current_phi),
            0.0,
            p_r * sin(PI/2.0) * sin(current_phi)
        );

        // Distance to ray (simplified)
        float dist_to_particle = length(cross(ray_dir, p_pos)) / length(ray_dir);

        if (dist_to_particle < PARTICLE_SIZE && abs(ray_r - p_r) < 2.0) {
            float particle_brightness = 1.0 - (dist_to_particle / PARTICLE_SIZE);
            particle_brightness = pow(particle_brightness, 3.0);

            // Hot particle color
            vec3 p_color = vec3(1.0, 0.7, 0.3);
            particle_color += p_color * particle_brightness * 0.5;
        }
    }

    return particle_color;
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

    vec3 accumulated_color = vec3(0.0);

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;
        float theta = state.theta;
        float phi = state.phi;

        // Particle contribution along ray
        vec3 ray_dir = vec3(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
        accumulated_color += get_particle_contribution(ray_dir, r) * 0.02;

        // Captured
        if (r < r_horizon * 1.01) {
            return accumulated_color;
        }

        // Escaped - starfield
        if (r > escape_radius) {
            return accumulated_color + get_starfield(ray_dir);
        }

        // Disk intersection
        float r_intersect, phi_intersect;
        if (check_disk_intersection(state_old, state, r_intersect, phi_intersect)) {
            vec3 disk_color = get_disk_color_animated(r_intersect, phi_intersect, state.v_phi);
            return accumulated_color + disk_color;
        }
    }

    return accumulated_color;
}

// Simple bloom effect
vec3 apply_bloom(vec3 color, vec2 uv) {
    vec3 bloom = vec3(0.0);

    // Sample surrounding pixels
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            vec2 offset = vec2(float(x), float(y)) / u_resolution * 3.0;
            vec3 sample_col = trace_ray(uv + offset);

            // Only bloom bright areas
            float brightness = max(max(sample_col.r, sample_col.g), sample_col.b);
            if (brightness > 1.0) {
                bloom += sample_col * 0.04;
            }
        }
    }

    return color + bloom * 0.5;
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    // Main ray trace
    vec3 color = trace_ray(pixel_coord);

    // Tone mapping
    color = color / (color + vec3(1.0));

    // Simplified color grading (removed bloom for performance)
    color = pow(color, vec3(0.85));

    // Vignette
    float dist = length(pixel_coord);
    float vignette = 1.0 - smoothstep(0.7, 1.5, dist);
    color *= vignette * 0.3 + 0.7;

    FragColor = vec4(color, 1.0);
}
