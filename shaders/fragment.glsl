#version 330 core

/**
 * SAGITTARIUS A* - SUPERMASSIVE BLACK HOLE SIMULATION
 * Based on Event Horizon Telescope (EHT) observations (2022)
 *
 * Physical Properties:
 * - Mass: 4.3 million solar masses
 * - Distance: 27,000 light-years from Earth
 * - Rotating black hole (Kerr metric)
 * - First direct image captured by EHT in 2022
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Physical constants (natural units: G=M=c=1)
const float M = 1.0;                      // Black hole mass (normalized)
const float a = 0.7;                      // Kerr spin parameter (0.7 = moderate-high rotation)
const float r_s = 2.0 * M;                // Schwarzschild radius
const float r_horizon = M * (1.0 + sqrt(1.0 - a*a/(M*M)));  // Kerr event horizon
const float r_photon_sphere = 3.0 * M;    // Approximate for Kerr
const float r_isco = M * (3.0 + sqrt(9.0 - 8.0 * a*a/(M*M)));  // Kerr ISCO
const float r_disk_inner = r_isco;        // Disk starts at ISCO
const float r_disk_outer = 25.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 280;                // Optimized for performance
const float STEP_SIZE = 0.22;
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

// Milky Way starfield (view from Sgr A* at galactic center)
vec3 get_starfield(vec3 dir) {
    // Convert direction to spherical coordinates
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI) + 0.5, theta / PI);

    vec3 color = vec3(0.02, 0.01, 0.03);  // Deep space background

    // Realistic starfield (view from galactic center)
    for (int i = 0; i < 80; i++) {
        vec2 starPos = hash22(vec2(float(i), 0.0));
        float dist = length(uv - starPos);

        if (dist < 0.004) {
            float brightness = pow(1.0 - (dist / 0.004), 3.0);
            // Realistic star colors (blue-white to orange)
            vec3 starColor = mix(vec3(0.8, 0.85, 1.0), vec3(1.0, 0.8, 0.5),
                                hash12(vec2(float(i), 1.0)));
            color += starColor * brightness * 0.8;
        }
    }

    // Galactic plane glow (we're at the center of the Milky Way!)
    float galactic_glow = pow(abs(sin(theta * 1.8)), 4.0) * 0.12;
    color += vec3(0.15, 0.12, 0.18) * galactic_glow;

    // Distant nebulae hints
    float nebula = hash12(uv * 3.0) * 0.03;
    color += vec3(0.15, 0.08, 0.12) * nebula;

    return color;
}

// EHT-style orange/yellow glow (based on 2022 Sgr A* image)
vec3 eht_color(float temp, float brightness) {
    // The EHT image shows predominantly orange/yellow/amber tones
    // Not blue-white like stellar objects

    // Base color: warm orange/amber
    vec3 base = vec3(1.0, 0.45, 0.1);  // Deep orange
    vec3 mid = vec3(1.0, 0.7, 0.25);   // Amber
    vec3 hot = vec3(1.0, 0.85, 0.5);   // Yellow-white

    vec3 color;
    if (temp > 0.65) {
        color = mix(mid, hot, (temp - 0.65) / 0.35);
    } else {
        color = mix(base, mid, temp / 0.65);
    }

    // Add subtle color variation based on brightness
    color = mix(color, vec3(1.0, 0.6, 0.2), brightness * 0.2);

    return color;
}

// Turbulence function for realistic disk structure
float turbulence(float r, float phi, float time) {
    // EHT observations show material swirls rapidly around Sgr A*
    float angle = phi + time * 0.3 / sqrt(r);

    // Multi-scale turbulent pattern
    float turb = 0.0;
    turb += 0.5 * sin(angle * 5.0 + r * 2.0);
    turb += 0.3 * sin(angle * 11.0 - r * 3.5 + time * 0.5);
    turb += 0.2 * sin(angle * 17.0 + r * 1.3 - time * 0.8);

    // Radial variation
    float radial = sin(r * 4.0 - time * 0.2) * 0.3;

    return 0.7 + (turb + radial) * 0.3;
}

// Sgr A* accretion disk with EHT-style appearance
vec3 get_disk_color(float r, float phi, float v_phi_photon) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Temperature profile: T ∝ r^(-3/4) (standard accretion disk model)
    float temp_ratio = pow(r_disk_inner / r, 0.75);

    // Kerr disk velocity with frame dragging
    float v_phi_disk = sqrt(M / r) * (1.0 + a * 0.35 / r);  // Frame dragging boost

    // RELATIVISTIC DOPPLER SHIFT
    // Critical for the lopsided appearance seen in EHT image
    float doppler_factor = 1.0 + (v_phi_disk - v_phi_photon) * 0.55;

    // RELATIVISTIC BEAMING: Brightness ∝ (doppler_factor)^4
    // This creates the bright crescent in the EHT image
    float beaming = pow(max(0.08, doppler_factor), 4.2);

    // Observed temperature (gravitationally redshifted/blueshifted)
    float temp_observed = temp_ratio / doppler_factor;
    temp_observed = clamp(temp_observed, 0.0, 1.0);

    // Turbulent structure (Sgr A* is highly variable)
    float turb = turbulence(r, phi, u_time);

    // Enhanced brightness near ISCO (inner edge glow)
    float isco_boost = 1.0 + 1.5 * exp(-3.0 * (r - r_disk_inner));

    // Final brightness combining all effects
    float brightness = beaming * turb * temp_ratio * isco_boost * 2.5;
    brightness = clamp(brightness, 0.0, 12.0);

    // Get EHT-style orange/amber color
    vec3 color = eht_color(temp_observed, brightness / 12.0);

    // Add photon ring glow (concentrated near photon sphere)
    float ring_factor = exp(-5.0 * abs(r - r_photon_sphere) / r_photon_sphere);
    color += vec3(1.0, 0.75, 0.4) * ring_factor * 2.0;

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

// KERR METRIC geodesic equations (rotating black hole)
// Implements frame dragging and ergosphere effects
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

    // Kerr metric functions
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);
    float sin2 = sin_theta * sin_theta;
    float cos2 = cos_theta * cos_theta;

    float rho2 = r * r + a * a * cos2;
    float delta = r * r - r_s * r + a * a;
    float sigma = (r * r + a * a) * (r * r + a * a) - a * a * delta * sin2;

    // Simplified Kerr Christoffel symbols (Boyer-Lindquist coordinates)
    // These capture the essential physics: frame dragging and rotation

    // Time component (frame dragging effect)
    float Gamma_t_rphi = (2.0 * a * r * r_s) / (rho2 * rho2);
    float Gamma_phi_rt = Gamma_t_rphi;

    // Radial components (modified by rotation)
    float Gamma_r_tt = (r_s * (r * r - a * a * cos2)) / (2.0 * rho2 * rho2);
    float Gamma_r_rr = (r - r_s / 2.0) / (delta * rho2);
    float Gamma_r_thetatheta = -r * delta / rho2;
    float Gamma_r_phiphi = -(delta * sin2 / rho2) * (r - (r_s * r * sin2 * (r * r + a * a)) / rho2);

    // Angular components
    float Gamma_theta_rtheta = r / rho2;
    float Gamma_theta_phiphi = -(sin_theta * cos_theta / rho2) * ((r * r + a * a) + (r_s * r * a * a * sin2) / rho2);
    float Gamma_phi_rphi = (r * (r * r + a * a) - r_s * r * a * a * sin2 / rho2) / (rho2 * sigma);
    float Gamma_phi_thetaphi = (cos_theta / sin_theta) * (1.0 + (r_s * r * a * a) / (rho2 * rho2));

    // Geodesic equations with Kerr corrections
    dydt.v_t = -2.0 * Gamma_t_rphi * v_r * v_phi;
    dydt.v_r = -Gamma_r_tt * v_t * v_t - Gamma_r_rr * v_r * v_r
               - Gamma_r_thetatheta * v_theta * v_theta
               - Gamma_r_phiphi * v_phi * v_phi;
    dydt.v_theta = -2.0 * Gamma_theta_rtheta * v_r * v_theta
                   - Gamma_theta_phiphi * v_phi * v_phi;
    dydt.v_phi = -2.0 * Gamma_phi_rphi * v_r * v_phi
                 - 2.0 * Gamma_phi_thetaphi * v_theta * v_phi
                 - 2.0 * Gamma_phi_rt * v_t * v_r;  // Frame dragging!

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
