#version 330 core

/**
 * Proper General Relativistic Black Hole Simulation
 * Following the technical guide exactly
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Physical constants (natural units: G=M=c=1)
const float r_s = 2.0;                    // Schwarzschild radius
const float r_horizon = 2.0;              // Event horizon
const float r_photon_sphere = 3.0;        // Photon sphere at 1.5*r_s
const float r_isco = 6.0;                 // ISCO at 3*r_s
const float r_disk_inner = 6.0;           // Disk starts at ISCO
const float r_disk_outer = 20.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 250;
const float STEP_SIZE = 0.25;
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

// Starfield with proper sphere mapping
vec3 get_starfield(vec3 dir) {
    // Convert direction to spherical coordinates
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.0);

    // Dense starfield
    for (int i = 0; i < 100; i++) {
        vec2 starPos = hash22(vec2(float(i), 0.0));
        float dist = length(uv - starPos);
        
        if (dist < 0.005) {
            float brightness = pow(1.0 - (dist / 0.005), 2.0);
            vec3 starColor = mix(vec3(0.7, 0.8, 1.0), vec3(1.0, 0.9, 0.6), 
                                hash12(vec2(float(i), 1.0)));
            color += starColor * brightness;
        }
    }

    // Milky Way glow
    float glow = pow(abs(sin(theta * 2.0)), 3.0) * 0.15;
    color += vec3(0.2, 0.15, 0.25) * glow;

    return color;
}

// Blackbody color from temperature (simplified)
vec3 blackbody_color(float temp) {
    // temp is normalized 0-1
    if (temp > 0.8) {
        return mix(vec3(1.0, 0.9, 0.7), vec3(0.7, 0.9, 1.0), (temp - 0.8) / 0.2);
    } else if (temp > 0.5) {
        return mix(vec3(1.0, 0.6, 0.2), vec3(1.0, 0.9, 0.7), (temp - 0.5) / 0.3);
    } else {
        return mix(vec3(0.8, 0.2, 0.1), vec3(1.0, 0.6, 0.2), temp / 0.5);
    }
}

// CRITICAL: Relativistic disk color with proper (1+z)^4 beaming
vec3 get_disk_color(float r, float phi, float v_phi_photon) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature profile: T ∝ r^(-3/4)
    float temp_ratio = pow(r_disk_inner / r, 0.75);
    
    // Disk fluid velocity (Keplerian orbit)
    float v_phi_disk = sqrt(1.0 / r) * 0.5;  // Simplified Keplerian
    
    // RELATIVISTIC DOPPLER SHIFT
    // The key to the lopsided appearance!
    float doppler_factor = 1.0 + (v_phi_disk - v_phi_photon) * 0.5;
    
    // RELATIVISTIC BEAMING: Brightness ∝ (doppler_factor)^4
    // This is THE critical effect from the guide!
    float beaming = pow(max(0.1, doppler_factor), 4.0);
    
    // Observed temperature (redshifted/blueshifted)
    float temp_observed = temp_ratio / doppler_factor;
    temp_observed = clamp(temp_observed, 0.0, 1.5);
    
    // Get color from temperature
    vec3 color = blackbody_color(temp_observed);
    
    // Rotation pattern
    float rotation = phi + u_time * 0.4 / r;
    float pattern = 0.8 + 0.2 * sin(rotation * 8.0);
    
    // Final brightness with beaming (THIS IS KEY!)
    float brightness = beaming * pattern * temp_ratio * 2.0;
    brightness = clamp(brightness, 0.0, 10.0);
    
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

// EXACT Schwarzschild geodesic equations from the guide
State geodesic_derivatives(State s) {
    State dydt;
    float r = s.r;
    float theta = s.theta;
    float v_t = s.v_t;
    float v_r = s.v_r;
    float v_theta = s.v_theta;
    float v_phi = s.v_phi;

    // Position derivatives = velocities
    dydt.t = v_t;
    dydt.r = v_r;
    dydt.theta = v_theta;
    dydt.phi = v_phi;

    // Christoffel symbols (exact from guide)
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

    // Velocity derivatives = accelerations (geodesic equation)
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

// RK4 integrator
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

// Initialize null geodesic (photon path)
State initialize_ray(vec2 pixel_coord) {
    State state;
    float r_cam = u_camera_pos.x;
    float theta_cam = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    state.t = 0.0;
    state.r = r_cam;
    state.theta = theta_cam;
    state.phi = phi_cam;

    // Ray direction based on pixel
    float tan_fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;
    float offset_x = pixel_coord.x * tan_fov * aspect;
    float offset_y = pixel_coord.y * tan_fov;

    float v_r = -1.0;
    float v_theta = offset_y * 0.5;
    float v_phi = offset_x * 0.5;

    // Normalize to null geodesic (photon constraint: ds^2 = 0)
    float f = 1.0 - r_s / r_cam;
    float sin_theta_cam = sin(theta_cam);
    float spatial = (v_r * v_r) / f + r_cam * r_cam *
                   (v_theta * v_theta + sin_theta_cam * sin_theta_cam * v_phi * v_phi);

    state.v_t = sqrt(max(0.0, spatial / f));
    state.v_r = v_r;
    state.v_theta = v_theta;
    state.v_phi = v_phi;

    return state;
}

// Check if ray crosses accretion disk
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

// Main ray tracing (GRRT algorithm from guide)
vec3 trace_ray(vec2 pixel_coord) {
    State state = initialize_ray(pixel_coord);
    State state_old = state;

    for (int step = 0; step < MAX_STEPS; step++) {
        state_old = state;
        state = rk4_step(state, STEP_SIZE);

        float r = state.r;

        // Condition A: Captured by event horizon
        if (r < r_horizon * 1.01) {
            return vec3(0.0);  // Black
        }

        // Condition B: Escaped to infinity
        if (r > escape_radius) {
            vec3 ray_dir = normalize(vec3(
                sin(state.theta) * cos(state.phi),
                cos(state.theta),
                sin(state.theta) * sin(state.phi)
            ));
            return get_starfield(ray_dir);
        }

        // Condition C: Hit accretion disk
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

    // Tone mapping for HDR
    color = color / (color + vec3(1.0));
    
    // Gamma correction
    color = pow(color, vec3(0.85));

    // Subtle vignette
    float dist = length(pixel_coord);
    float vignette = 1.0 - smoothstep(1.0, 2.0, dist);
    color *= vignette * 0.3 + 0.7;

    FragColor = vec4(color, 1.0);
}
