#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Physical constants
const float rs = 2.0;
const float r_inner = 2.5 * rs;
const float r_outer = 10.0 * rs;
const float r_photon = 1.5 * rs;  // Photon sphere
const float MAX_R = 100.0;
const int MAX_ITER = 120;  // Slightly more for better lensing
const float DL = 0.22;
const float PI = 3.14159265359;

// Better hash for stars
float hash(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

vec3 stars(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 col = vec3(0.005, 0.008, 0.012);  // Deep space

    // Varied stars
    for(int i = 0; i < 50; i++) {
        vec2 p = vec2(hash(vec2(float(i), 0.0)), hash(vec2(float(i), 1.0)));
        float d = length(uv - p);
        float size = 0.002 + hash(vec2(float(i), 2.0)) * 0.002;

        if(d < size) {
            float bright = pow(1.0 - d/size, 2.0);
            // Star colors: blue-white to yellow-orange
            float temp = hash(vec2(float(i), 3.0));
            vec3 starCol = mix(vec3(0.8, 0.9, 1.0), vec3(1.0, 0.9, 0.6), temp);
            col += starCol * bright;
        }
    }

    // Subtle galactic glow
    float glow = pow(abs(sin(theta * 2.0)), 4.0) * 0.08;
    col += vec3(0.15, 0.12, 0.2) * glow;

    return col;
}

vec3 diskColor(float r, float phi, bool isSecondary) {
    // Temperature-based color (realistic blackbody)
    float temp = pow(r_inner / r, 0.75);  // T ∝ r^(-3/4)

    vec3 color;
    if(temp > 1.2) {
        // Ultra hot: blue-white
        color = vec3(0.9, 0.95, 1.0);
    } else if(temp > 0.8) {
        // Very hot: white to yellow
        color = mix(vec3(1.0, 1.0, 1.0), vec3(1.0, 0.9, 0.6), (temp - 0.8) / 0.4);
    } else if(temp > 0.5) {
        // Hot: yellow to orange
        color = mix(vec3(1.0, 0.9, 0.6), vec3(1.0, 0.6, 0.3), (temp - 0.5) / 0.3);
    } else {
        // Warm: orange to red
        color = mix(vec3(1.0, 0.6, 0.3), vec3(0.9, 0.3, 0.1), temp / 0.5);
    }

    // Brightness with realistic falloff
    float bright = 6.0 * pow(r_inner / r, 1.4);
    bright = clamp(bright, 0.0, 12.0);

    // Turbulent pattern (spiral structure)
    float spiral = sin(phi * 8.0 - log(r) * 5.0 + u_time * 0.3) * 0.5 + 0.5;
    float rotation = sin(phi * 12.0 - u_time * 0.5 / sqrt(r)) * 0.3 + 0.7;

    // RELATIVISTIC EFFECTS
    // Doppler shift + beaming
    float v = sqrt(1.0 / r);  // Keplerian velocity
    float doppler = 1.0 + v * sin(phi) * 0.7;
    float beaming = pow(max(0.3, doppler), 3.5);  // Relativistic beaming

    // Apply effects
    bright *= spiral * rotation * beaming;

    // Secondary image is dimmer (went around the back)
    if(isSecondary) {
        bright *= 0.4;
    }

    return color * bright;
}

void main() {
    vec2 uv = (TexCoord * 2.0 - 1.0);
    uv.y = -uv.y;

    // Camera
    float r = u_camera_pos.x;
    float theta = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    // Ray direction
    float fov = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    float dr = -1.0;
    float dtheta = uv.y * fov;
    float dphi = uv.x * fov * aspect;

    // Ray position
    float theta_ray = theta;
    float phi_ray = phi_cam;

    vec3 color = vec3(0.0);
    vec3 diskContribution = vec3(0.0);
    int diskHits = 0;
    bool capturedByHole = false;

    for(int i = 0; i < MAX_ITER; i++) {
        // Schwarzschild geodesic
        float f = 1.0 - rs / r;
        float sin_t = sin(theta_ray);
        float cos_t = cos(theta_ray);

        // Accelerations
        float acc_r = rs / (2.0 * r * r) - f * (dtheta * dtheta + sin_t * sin_t * dphi * dphi);
        float acc_theta = -2.0 * dr * dtheta / r + sin_t * cos_t * dphi * dphi;
        float acc_phi = -2.0 * dr * dphi / r - 2.0 * cos_t / max(sin_t, 0.01) * dtheta * dphi;

        // Update
        dr += acc_r * DL;
        dtheta += acc_theta * DL;
        dphi += acc_phi * DL;

        float r_old = r;
        float theta_old = theta_ray;

        r += dr * DL;
        theta_ray += dtheta * DL;
        phi_ray += dphi * DL;

        // Check event horizon
        if(r <= rs * 1.05) {
            capturedByHole = true;
            break;
        }

        // Check escape
        if(r > MAX_R) {
            vec3 dir = normalize(vec3(
                sin(theta_ray) * cos(phi_ray),
                cos(theta_ray),
                sin(theta_ray) * sin(phi_ray)
            ));
            color = stars(dir);
            break;
        }

        // Check disk intersection (multiple passes for lensing!)
        float eq = PI / 2.0;
        if((theta_old - eq) * (theta_ray - eq) < 0.0) {
            float t = (eq - theta_old) / (theta_ray - theta_old);
            float r_hit = r_old + t * (r - r_old);

            if(r_hit >= r_inner && r_hit <= r_outer && diskHits < 2) {
                float phi_hit = phi_ray;
                bool isSecondary = (diskHits == 1);  // Second crossing is dimmer
                diskContribution += diskColor(r_hit, phi_hit, isSecondary);
                diskHits++;
            }
        }

        // Add photon ring glow
        if(r > r_photon * 0.9 && r < r_photon * 1.15) {
            float ringGlow = exp(-15.0 * abs(r - r_photon) / r_photon);
            diskContribution += vec3(1.0, 0.85, 0.6) * ringGlow * 1.5;
        }
    }

    // Combine contributions
    if(capturedByHole) {
        color = vec3(0.0);  // Pure black
    }

    color += diskContribution;

    // HDR tone mapping
    color = color / (color + vec3(0.6));

    // Contrast enhancement
    color = pow(color, vec3(0.88));
    color = (color - 0.5) * 1.15 + 0.5;
    color = clamp(color, 0.0, 1.0);

    // Subtle vignette
    float dist = length(uv);
    float vignette = 1.0 - smoothstep(0.9, 1.6, dist);
    color *= vignette * 0.5 + 0.5;

    FragColor = vec4(color, 1.0);
}
