#version 330 core

/**
 * Optimized Black Hole - Faster with better visuals
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
const int MAX_STEPS = 150;  // Reduced for speed
const float STEP_SIZE = 0.3;  // Larger steps = faster
const float PI = 3.14159265359;

struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Better hash
float hash12(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

vec2 hash22(vec2 p) {
    vec3 p3 = fract(vec3(p.xyx) * vec3(0.1031, 0.1030, 0.0973));
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.xx + p3.yz) * p3.zy);
}

// Noise
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

// Beautiful starfield
vec3 get_starfield(vec3 dir) {
    float theta = acos(dir.y);
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.0);

    // Bright stars
    for (int i = 0; i < 60; i++) {
        vec2 starPos = hash22(vec2(float(i), 0.0));
        float dist = length(uv - starPos);

        if (dist < 0.004) {
            float brightness = 1.0 - (dist / 0.004);
            brightness = pow(brightness, 2.0);
            
            // Color variation
            float hue = hash12(vec2(float(i), 1.0));
            vec3 starColor = mix(vec3(0.8, 0.9, 1.0), vec3(1.0, 0.9, 0.7), hue);
            color += starColor * brightness * 1.5;
        }
    }

    // Nebula glow
    float nebula = noise(uv * 3.0 + u_time * 0.01);
    nebula += 0.5 * noise(uv * 6.0);
    color += vec3(0.1, 0.05, 0.15) * nebula * 0.3;

    return color;
}

// DRAMATIC disk with better colors
vec3 get_disk_color(float r, float phi) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Distance from inner edge
    float norm_r = (r - r_disk_inner) / (r_disk_outer - r_disk_inner);

    // Rotation
    float rotation = phi + u_time * 0.5 / (r * 0.5);
    
    // Turbulence
    vec2 turb_uv = vec2(r * 0.3, rotation);
    float turbulence = noise(turb_uv * 5.0);
    turbulence += 0.5 * noise(turb_uv * 10.0);

    // Temperature gradient - MUCH hotter near black hole
    float temp = pow(1.0 - norm_r, 1.5);
    
    // Color based on temperature
    vec3 color;
    if (temp > 0.7) {
        // Inner: white-hot to blue-white
        color = mix(vec3(1.0, 0.95, 0.8), vec3(0.9, 0.95, 1.0), (temp - 0.7) / 0.3);
    } else if (temp > 0.4) {
        // Middle: yellow-white
        color = mix(vec3(1.0, 0.7, 0.3), vec3(1.0, 0.95, 0.8), (temp - 0.4) / 0.3);
    } else {
        // Outer: orange to red
        color = mix(vec3(1.0, 0.4, 0.1), vec3(1.0, 0.7, 0.3), temp / 0.4);
    }

    // Brightness with turbulence
    float brightness = pow(1.0 - norm_r, 2.5) * (1.0 + turbulence * 0.5);
    brightness *= 1.0 + 0.3 * sin(rotation * 12.0);
    brightness = clamp(brightness, 0.0, 5.0);

    return color * brightness;
}

// State math
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

// Geodesic equations
State geodesic_derivatives(State s) {
    State dydt;
    float r = s.r;
    float theta = s.theta;
    float v_t = s.v_t;
    float v_r = s.v_r;
    float v_theta = s.v_theta;
    float v_phi = s.v_phi;

    dydt.t = v_t;
    dydt.r = v_r;
    dydt.theta = v_theta;
    dydt.phi = v_phi;

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
    float Gamma_phi_thetaphi = cos_theta / max(sin_theta, 0.001);

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

// RK4
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

// Initialize
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
    float spatial_term = (v_r * v_r) / f + r_cam * r_cam *
                        (v_theta * v_theta + sin_theta * sin_theta * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial_term / f));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

// Disk intersection
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

// Main tracing
vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    vec3 accumulated = vec3(0.0);

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;

        // Captured
        if (r < r_horizon * 1.01) {
            return accumulated;
        }

        // Escaped
        if (r > escape_radius) {
            vec3 ray_dir = vec3(
                sin(state.theta) * cos(state.phi),
                cos(state.theta),
                sin(state.theta) * sin(state.phi)
            );
            return accumulated + get_starfield(ray_dir);
        }

        // Disk
        float r_intersect, phi_intersect;
        if (check_disk_intersection(state_old, state, r_intersect, phi_intersect)) {
            return accumulated + get_disk_color(r_intersect, phi_intersect);
        }
    }

    return accumulated;
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    vec3 color = trace_ray(pixel_coord);

    // HDR tone mapping
    color = color / (color + vec3(1.0));
    
    // Gamma
    color = pow(color, vec3(0.8));

    // Boost contrast
    color = (color - 0.5) * 1.2 + 0.5;
    color = clamp(color, 0.0, 1.0);

    // Subtle vignette
    float dist = length(pixel_coord);
    float vignette = 1.0 - smoothstep(0.8, 1.8, dist);
    color *= vignette * 0.4 + 0.6;

    FragColor = vec4(color, 1.0);
}
