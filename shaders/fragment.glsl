#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Physical constants - keeping it simple
const float rs = 2.0;              // Schwarzschild radius
const float r_inner = 2.2 * rs;    // Disk inner (like the reference)
const float r_outer = 8.0 * rs;    // Disk outer
const float MAX_R = 100.0;
const int MAX_ITER = 100;          // MUCH fewer iterations
const float DL = 0.25;             // Step size
const float PI = 3.14159265359;

vec3 stars(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 col = vec3(0.0);

    // Quick stars
    for(int i = 0; i < 30; i++) {
        vec2 p = vec2(fract(sin(float(i)*12.9898)*43758.5453),
                     fract(sin(float(i)*78.233)*43758.5453));
        float d = length(uv - p);
        if(d < 0.005) col += vec3(1.0) * (1.0 - d/0.005);
    }

    return col;
}

vec3 diskColor(float r, float phi) {
    // Brightness falloff
    float bright = 5.0 * pow(r_inner / r, 1.5);
    bright = clamp(bright, 0.0, 10.0);

    // Simple animated pattern
    float anim = sin(phi * 10.0 - u_time * 0.5 + r * 0.5) * 0.5 + 0.5;

    // Orange/yellow gradient
    float t = (r - r_inner) / (r_outer - r_inner);
    vec3 hot = vec3(1.0, 0.9, 0.6);   // Yellow-white
    vec3 cool = vec3(1.0, 0.4, 0.1);  // Orange
    vec3 col = mix(hot, cool, t);

    // Doppler boost (one side brighter)
    float doppler = 1.0 + 0.6 * cos(phi);
    bright *= pow(doppler, 2.5);

    return col * bright * anim;
}

void main() {
    vec2 uv = (TexCoord * 2.0 - 1.0);
    uv.y = -uv.y;

    // Camera setup
    float r = u_camera_pos.x;
    float theta = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    // Ray direction
    float fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    float dr = -1.0;
    float dtheta = uv.y * fov;
    float dphi = uv.x * fov * aspect;

    // Ray trace
    float theta_ray = theta;
    float phi_ray = phi_cam;

    vec3 color = vec3(0.0);
    bool hit = false;

    for(int i = 0; i < MAX_ITER && !hit; i++) {
        // Simple Schwarzschild geodesic
        float f = 1.0 - rs / r;
        float sin_t = sin(theta_ray);
        float cos_t = cos(theta_ray);

        // Accelerations (simplified but correct)
        float acc_r = rs / (2.0 * r * r) - f * (dtheta * dtheta + sin_t * sin_t * dphi * dphi);
        float acc_theta = -2.0 * dr * dtheta / r + sin_t * cos_t * dphi * dphi;
        float acc_phi = -2.0 * dr * dphi / r - 2.0 * cos_t / max(sin_t, 0.01) * dtheta * dphi;

        // Update velocities
        dr += acc_r * DL;
        dtheta += acc_theta * DL;
        dphi += acc_phi * DL;

        // Update positions
        r += dr * DL;
        theta_ray += dtheta * DL;
        phi_ray += dphi * DL;

        // Check if hit event horizon
        if(r <= rs * 1.05) {
            color = vec3(0.0);  // BLACK
            hit = true;
        }

        // Check if escaped
        else if(r > MAX_R) {
            vec3 dir = normalize(vec3(
                sin(theta_ray) * cos(phi_ray),
                cos(theta_ray),
                sin(theta_ray) * sin(phi_ray)
            ));
            color = stars(dir);
            hit = true;
        }

        // Check if hit disk (at equator)
        else if(abs(theta_ray - PI/2.0) < 0.1) {
            if(r >= r_inner && r <= r_outer) {
                color = diskColor(r, phi_ray);
                hit = true;
            }
        }
    }

    // Simple tone mapping
    color = color / (color + 0.5);
    color = pow(color, vec3(0.9));

    FragColor = vec4(color, 1.0);
}
