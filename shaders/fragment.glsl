#version 330 core

/**
 * SAGITTARIUS A* - THE REAL DEAL
 * Visually stunning black hole with proper gravitational lensing
 */

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float M = 1.0;
const float r_s = 2.0;
const float r_horizon = 1.8;
const float r_photon = 2.6;
const float r_isco = 2.8;
const float r_disk_inner = 3.0;
const float r_disk_outer = 12.0;
const float escape_radius = 100.0;
const int MAX_STEPS = 250;
const float STEP_SIZE = 0.12;
const float PI = 3.14159265359;

struct Ray {
    float t, r, theta, phi;
    float v_t, v_r, v_theta, v_phi;
};

// Better hash
float hash(vec2 p) {
    p = fract(p * vec2(123.45, 678.90));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

// Starfield with variety
vec3 starfield(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 col = vec3(0.005, 0.008, 0.015);  // Deep space blue

    // Multiple star layers
    for (int i = 0; i < 80; i++) {
        vec2 starPos = vec2(hash(vec2(float(i), 0.0)), hash(vec2(float(i), 1.0)));
        float dist = length(uv - starPos);
        float size = 0.002 + hash(vec2(float(i), 2.0)) * 0.002;

        if (dist < size) {
            float bright = pow(1.0 - dist / size, 2.5);
            float temp = hash(vec2(float(i), 3.0));
            vec3 starCol = mix(vec3(0.7, 0.8, 1.0), vec3(1.0, 0.9, 0.7), temp);
            col += starCol * bright * 1.5;
        }
    }

    // Subtle nebula
    float nebula = hash(uv * 2.0) * 0.03;
    col += vec3(0.15, 0.08, 0.2) * nebula;

    return col;
}

// BRIGHT and dramatic accretion disk
vec3 disk_color(float r, float phi) {
    if (r < r_disk_inner || r > r_disk_outer) {
        return vec3(0.0);
    }

    // Super bright inner regions
    float brightness = pow(r_disk_inner / r, 1.5);
    brightness = clamp(brightness, 0.0, 5.0);

    // Spiral structure
    float spiral = sin(phi * 8.0 - u_time * 0.3 + log(r) * 5.0) * 0.5 + 0.5;
    float rotation = sin(phi * 15.0 + u_time * 0.8 / sqrt(r)) * 0.3 + 0.7;

    // Temperature gradient - VIVID colors
    float temp = (r - r_disk_inner) / (r_disk_outer - r_disk_inner);
    vec3 color;

    if (temp < 0.3) {
        // Ultra hot: blue-white to white
        color = mix(vec3(0.9, 0.95, 1.0), vec3(1.0, 1.0, 1.0), temp / 0.3);
    } else if (temp < 0.6) {
        // Hot: white to yellow
        color = mix(vec3(1.0, 1.0, 1.0), vec3(1.0, 0.9, 0.5), (temp - 0.3) / 0.3);
    } else {
        // Warm: yellow to orange to red
        color = mix(vec3(1.0, 0.9, 0.5), vec3(1.0, 0.3, 0.1), (temp - 0.6) / 0.4);
    }

    // Doppler boost - DRAMATIC asymmetry
    float v_disk = sqrt(M / r);
    float doppler = 1.0 + v_disk * sin(phi) * 0.8;
    float beaming = pow(max(0.3, doppler), 4.5);

    // Final brightness
    float final_brightness = brightness * spiral * rotation * beaming * 4.0;
    final_brightness = clamp(final_brightness, 0.0, 15.0);

    return color * final_brightness;
}

Ray ray_add(Ray a, Ray b) {
    return Ray(
        a.t + b.t, a.r + b.r, a.theta + b.theta, a.phi + b.phi,
        a.v_t + b.v_t, a.v_r + b.v_r, a.v_theta + b.v_theta, a.v_phi + b.v_phi
    );
}

Ray ray_mul(Ray r, float k) {
    return Ray(
        r.t * k, r.r * k, r.theta * k, r.phi * k,
        r.v_t * k, r.v_r * k, r.v_theta * k, r.v_phi * k
    );
}

// Schwarzschild geodesics with proper bending
Ray geodesic(Ray s) {
    Ray d;
    float r = s.r;
    float theta = s.theta;

    d.t = s.v_t;
    d.r = s.v_r;
    d.theta = s.v_theta;
    d.phi = s.v_phi;

    float sin_t = sin(theta);
    float cos_t = cos(theta);
    float sin2 = sin_t * sin_t;

    // Schwarzschild Christoffel symbols
    float f = 1.0 - r_s / r;

    d.v_t = r_s * s.v_t * s.v_r / (r * r * f);
    d.v_r = r_s * s.v_t * s.v_t / (2.0 * r * r) - r_s * s.v_r * s.v_r / (2.0 * r * r * f * f);
    d.v_r += f * (s.v_theta * s.v_theta + sin2 * s.v_phi * s.v_phi);
    d.v_theta = -2.0 * s.v_r * s.v_theta / r + sin_t * cos_t * s.v_phi * s.v_phi;
    d.v_phi = -2.0 * s.v_r * s.v_phi / r - 2.0 * cos_t * s.v_theta * s.v_phi / max(sin_t, 0.001);

    return d;
}

Ray rk4(Ray y, float h) {
    Ray k1 = ray_mul(geodesic(y), h);
    Ray k2 = ray_mul(geodesic(ray_add(y, ray_mul(k1, 0.5))), h);
    Ray k3 = ray_mul(geodesic(ray_add(y, ray_mul(k2, 0.5))), h);
    Ray k4 = ray_mul(geodesic(ray_add(y, k3)), h);

    return ray_add(y, ray_mul(ray_add(ray_add(k1, ray_mul(k2, 2.0)), ray_add(ray_mul(k3, 2.0), k4)), 1.0 / 6.0));
}

Ray init_ray(vec2 pixel) {
    Ray r;

    r.t = 0.0;
    r.r = u_camera_pos.x;
    r.theta = u_camera_pos.y;
    r.phi = u_camera_pos.z;

    float fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    r.v_r = -1.0;
    r.v_theta = pixel.y * fov;
    r.v_phi = pixel.x * fov * aspect;

    // Normalize for null geodesic
    float f = 1.0 - r_s / r.r;
    float sin_cam = sin(r.theta);
    float norm = (r.v_r * r.v_r / f) + r.r * r.r * (r.v_theta * r.v_theta + sin_cam * sin_cam * r.v_phi * r.v_phi);
    r.v_t = sqrt(max(0.0, norm / f));

    return r;
}

bool hit_disk(Ray old, Ray new, out float r_hit, out float phi_hit) {
    float eq = PI / 2.0;
    bool cross = (old.theta - eq) * (new.theta - eq) < 0.0;

    if (cross) {
        float t = (eq - old.theta) / (new.theta - old.theta);
        r_hit = old.r + t * (new.r - old.r);
        phi_hit = old.phi + t * (new.phi - old.phi);
        return r_hit >= r_disk_inner && r_hit <= r_disk_outer;
    }
    return false;
}

vec3 trace(vec2 pixel) {
    Ray ray = init_ray(pixel);
    Ray old = ray;

    for (int i = 0; i < MAX_STEPS; i++) {
        old = ray;
        ray = rk4(ray, STEP_SIZE);

        // Captured by black hole
        if (ray.r < r_horizon * 1.1) {
            return vec3(0.0);
        }

        // Escaped to infinity
        if (ray.r > escape_radius) {
            vec3 dir = normalize(vec3(
                sin(ray.theta) * cos(ray.phi),
                cos(ray.theta),
                sin(ray.theta) * sin(ray.phi)
            ));
            return starfield(dir);
        }

        // Hit accretion disk
        float r_hit, phi_hit;
        if (hit_disk(old, ray, r_hit, phi_hit)) {
            vec3 disk = disk_color(r_hit, phi_hit);

            // Add photon ring glow
            if (ray.r > r_photon * 0.8 && ray.r < r_photon * 1.2) {
                float ring_glow = exp(-10.0 * abs(ray.r - r_photon));
                disk += vec3(1.0, 0.8, 0.5) * ring_glow * 3.0;
            }

            return disk;
        }
    }

    return vec3(0.0);
}

void main() {
    vec2 uv = (2.0 * TexCoord - 1.0);
    uv.y = -uv.y;

    vec3 color = trace(uv);

    // HDR tone mapping
    color = color / (color + vec3(0.5));

    // Enhance contrast
    color = pow(color, vec3(0.85));
    color = (color - 0.5) * 1.2 + 0.5;
    color = clamp(color, 0.0, 1.0);

    // Subtle vignette
    float dist = length(uv);
    float vignette = 1.0 - smoothstep(0.8, 1.5, dist);
    color *= vignette * 0.4 + 0.6;

    FragColor = vec4(color, 1.0);
}
