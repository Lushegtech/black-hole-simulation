#version 330 core

/**
 * KERR BLACK HOLE - Rotating spacetime with frame dragging
 * Maximum visual impact - Interstellar-grade simulation
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Physical constants (natural units: G=M=c=1)
const float M = 1.0;                      // Black hole mass
const float a = 0.95;                     // Spin parameter (0.95 = fast spin!)
const float r_s = 2.0 * M;                // Schwarzschild radius
const float r_isco = M * (3.0 + sqrt(9.0 - 8.0 * a * a / (M * M)));  // Kerr ISCO
const float r_disk_inner = r_isco;
const float r_disk_outer = 20.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 300;
const float STEP_SIZE = 0.2;
const float PI = 3.14159265359;

struct State {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Hash functions
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

// Noise for turbulence
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

// Fractal Brownian Motion
float fbm(vec2 p) {
    float value = 0.0;
    float amplitude = 0.5;
    for (int i = 0; i < 4; i++) {
        value += amplitude * noise(p);
        p *= 2.0;
        amplitude *= 0.5;
    }
    return value;
}

// EPIC starfield with nebula
vec3 get_starfield(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.0);

    // Bright stars
    for (int i = 0; i < 120; i++) {
        vec2 starPos = hash22(vec2(float(i), 0.0));
        float dist = length(uv - starPos);
        
        if (dist < 0.006) {
            float brightness = pow(1.0 - (dist / 0.006), 2.5);
            float hue = hash12(vec2(float(i), 1.0));
            vec3 starColor = mix(vec3(0.6, 0.7, 1.0), vec3(1.0, 0.8, 0.5), hue);
            color += starColor * brightness * 1.2;
        }
    }

    // Nebula clouds
    float nebula = fbm(uv * 4.0 + u_time * 0.005);
    nebula += 0.5 * fbm(uv * 8.0);
    color += vec3(0.15, 0.08, 0.2) * nebula * 0.4;

    // Galactic band
    float band = pow(abs(sin(theta * 3.0)), 5.0) * 0.25;
    color += vec3(0.3, 0.2, 0.35) * band;

    return color;
}

// Proper blackbody with blue-white hot inner regions
vec3 blackbody_color(float temp) {
    temp = clamp(temp, 0.0, 1.5);
    
    if (temp > 1.2) {
        // Ultra-hot: blue-white
        return mix(vec3(0.8, 0.9, 1.0), vec3(0.6, 0.8, 1.0), (temp - 1.2) / 0.3);
    } else if (temp > 0.9) {
        // Very hot: white
        return mix(vec3(1.0, 0.95, 0.85), vec3(0.8, 0.9, 1.0), (temp - 0.9) / 0.3);
    } else if (temp > 0.6) {
        // Hot: yellow-white
        return mix(vec3(1.0, 0.75, 0.3), vec3(1.0, 0.95, 0.85), (temp - 0.6) / 0.3);
    } else if (temp > 0.3) {
        // Medium: orange
        return mix(vec3(1.0, 0.5, 0.15), vec3(1.0, 0.75, 0.3), (temp - 0.3) / 0.3);
    } else {
        // Cool: red-orange
        return mix(vec3(0.7, 0.2, 0.05), vec3(1.0, 0.5, 0.15), temp / 0.3);
    }
}

// KERR DISK: With frame dragging, Doppler, and beaming
vec3 get_disk_color(float r, float phi, float v_phi_photon) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature: T ∝ r^(-3/4) but BOOSTED near ISCO for Kerr
    float temp_ratio = pow(r_disk_inner / r, 0.8);  // Steeper for Kerr
    
    // Keplerian velocity WITH frame dragging effect
    float v_phi_disk = sqrt(M / r) * (1.0 + a * 0.3 / r);  // Frame dragging boost
    
    // RELATIVISTIC DOPPLER
    float doppler_factor = 1.0 + (v_phi_disk - v_phi_photon) * 0.6;
    
    // RELATIVISTIC BEAMING: (1+z)^4
    float beaming = pow(max(0.05, doppler_factor), 4.5);  // Even MORE dramatic
    
    // Observed temperature
    float temp_observed = temp_ratio / doppler_factor;
    temp_observed = clamp(temp_observed, 0.0, 1.8);
    
    // Turbulence
    vec2 turb_uv = vec2(r * 0.4, phi * 3.0);
    float turbulence = fbm(turb_uv + u_time * 0.1);
    
    // Rotation with frame dragging
    float rotation = phi + u_time * (0.5 / r + a * 0.2);  // Faster near hole!
    float pattern = 0.7 + 0.3 * sin(rotation * 10.0);
    pattern *= (1.0 + turbulence * 0.4);
    
    // Color
    vec3 color = blackbody_color(temp_observed);
    
    // EPIC brightness with beaming
    float brightness = beaming * pattern * temp_ratio * 3.0;
    brightness *= (1.0 + turbulence * 0.5);
    brightness = clamp(brightness, 0.0, 15.0);
    
    return color * brightness;
}

// State math
State state_add(State s1, State s2) {
    return State(
        s1.t + s2.t, s1.r + s2.r, s1.theta + s2.theta, s1.phi + s2.phi,
        s1.v_t + s2.v_t, s1.v_r + s2.v_r, s1.v_theta + s2.v_theta, s1.v_phi + s2.v_phi
    );
}

State state_mul(State s, float k) {
    return State(
        s.t * k, s.r * k, s.theta * k, s.phi * k,
        s.v_t * k, s.v_r * k, s.v_theta * k, s.v_phi * k
    );
}

// KERR geodesic equations (simplified Boyer-Lindquist)
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
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;
    
    // Kerr metric functions
    float rho2 = r * r + a * a * cos2;
    float delta = r * r - r_s * r + a * a;
    
    // Simplified Kerr Christoffel symbols (approximate)
    float Gamma_r_tt = (r_s * (r * r - a * a * cos2)) / (2.0 * rho2 * rho2);
    float Gamma_t_tr = (r_s * r) / (2.0 * rho2 * delta);
    float Gamma_r_rr = (r - M) / delta - r / rho2;
    float Gamma_r_thetatheta = -r * delta / rho2;
    float Gamma_r_phiphi = -sin2 * (r * delta - a * a * r_s * r * sin2) / rho2;
    float Gamma_theta_rtheta = 1.0 / r;
    float Gamma_theta_phiphi = -sin_theta * cos_theta * (1.0 + (a * a * r_s * r) / (rho2 * rho2));
    float Gamma_phi_rphi = (1.0 / r) * (1.0 + (a * a * r_s * r * sin2) / (rho2 * rho2));
    float Gamma_phi_thetaphi = cos_theta / max(sin_theta, 0.001) * (1.0 - (a * a * r_s * r * cos2) / (rho2 * rho2));

    // Geodesic accelerations
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
    return state_add(y, state_mul(increment, 1.0 / 6.0));
}

// Initialize photon
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

    // Null geodesic normalization (approximate for Kerr)
    float rho2 = r_cam * r_cam + a * a * cos(theta_cam) * cos(theta_cam);
    float delta = r_cam * r_cam - r_s * r_cam + a * a;
    float f = delta / rho2;
    
    float sin_theta_cam = sin(theta_cam);
    float spatial = (v_r * v_r) / f + r_cam * r_cam *
                   (v_theta * v_theta + sin_theta_cam * sin_theta_cam * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial * rho2 / delta));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

// Disk intersection
bool check_disk_intersection(State s_old, State s_new, out float r_int, out float phi_int, out float v_phi_int) {
    float pi_2 = PI / 2.0;
    bool crossed = (s_old.theta - pi_2) * (s_new.theta - pi_2) < 0.0;

    if (crossed) {
        float t = (pi_2 - s_old.theta) / (s_new.theta - s_old.theta);
        r_int = s_old.r + t * (s_new.r - s_old.r);
        phi_int = s_old.phi + t * (s_new.phi - s_old.phi);
        v_phi_int = s_old.v_phi + t * (s_new.v_phi - s_old.v_phi);
        return r_int >= r_disk_inner && r_int <= r_disk_outer;
    }
    return false;
}

// Main GRRT
vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;

        // Event horizon (Kerr)
        float r_horizon = M + sqrt(M * M - a * a);
        if (r < r_horizon * 1.05) {
            return vec3(0.0);
        }

        // Escaped
        if (r > escape_radius) {
            vec3 ray_dir = normalize(vec3(
                sin(state.theta) * cos(state.phi),
                cos(state.theta),
                sin(state.theta) * sin(state.phi)
            ));
            return get_starfield(ray_dir);
        }

        // Disk hit
        float r_int, phi_int, v_phi_photon;
        if (check_disk_intersection(state_old, state, r_int, phi_int, v_phi_photon)) {
            return get_disk_color(r_int, phi_int, v_phi_photon);
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    vec3 color = trace_ray(pixel_coord);

    // HDR tone mapping
    color = color / (color + vec3(0.8));
    
    // Contrast boost
    color = pow(color, vec3(0.82));
    color = (color - 0.5) * 1.15 + 0.5;
    color = clamp(color, 0.0, 1.0);

    // Cinematic vignette
    float dist = length(pixel_coord);
    float vignette = 1.0 - smoothstep(0.9, 1.9, dist);
    color *= vignette * 0.35 + 0.65;

    FragColor = vec4(color, 1.0);
}
