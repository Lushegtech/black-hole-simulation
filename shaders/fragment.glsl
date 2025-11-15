#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Sagittarius A* parameters (simplified but accurate)
const float rs = 2.0;                    // Schwarzschild radius
const float r_inner = 2.2 * rs;          // Like the reference!
const float r_outer = 5.5 * rs;          // Smaller disk
const float r_photon = 1.5 * rs;         // Photon sphere
const float MAX_R = 100.0;
const int MAX_ITER = 150;                // Good balance
const float DL = 0.18;                   // Smaller steps for better accuracy
const float PI = 3.14159265359;

// Better hash
float hash(vec2 p) {
    p = fract(p * 0.3183099 + 0.1);
    p *= 17.0;
    return fract(p.x * p.y * (p.x + p.y));
}

// Simple but beautiful starfield
vec3 stars(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 col = vec3(0.0, 0.0, 0.02);  // Deep space

    // Stars with variety
    for(int i = 0; i < 60; i++) {
        float fi = float(i);
        vec2 p = vec2(hash(vec2(fi * 0.1, 0.0)), hash(vec2(fi * 0.1, 1.0)));
        float d = length(uv - p);
        float size = 0.002 + hash(vec2(fi, 2.0)) * 0.0015;

        if(d < size) {
            float bright = pow(1.0 - d/size, 2.5);
            // Star color temperature
            float temp = hash(vec2(fi, 3.0));
            vec3 starCol = mix(vec3(0.7, 0.8, 1.0), vec3(1.0, 0.85, 0.6), temp);
            col += starCol * bright * 1.2;
        }
    }

    return col;
}

// Disk color - SIMPLE like the reference!
vec3 diskColor(float r, float phi) {
    // Normalized radius (0 to 1)
    float t = (r - r_inner) / (r_outer - r_inner);

    // Simple orange gradient (like reference: vec3(1.0, r, 0.2))
    vec3 color = vec3(1.0, 0.3 + t * 0.6, 0.1 + t * 0.3);

    // Brightness falloff
    float bright = 7.0 * pow(r_inner / r, 1.6);
    bright = clamp(bright, 0.0, 15.0);

    // Animated spiral
    float spiral = sin(phi * 10.0 - log(r + 1.0) * 6.0 + u_time * 0.4) * 0.3 + 0.7;

    // Doppler boost
    float v = sqrt(1.0 / r);
    float doppler = 1.0 + v * sin(phi) * 0.8;
    bright *= pow(max(0.2, doppler), 3.8);

    return color * bright * spiral;
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
    vec3 diskAcc = vec3(0.0);
    int hits = 0;
    bool done = false;

    // Ray trace with proper RK4
    for(int i = 0; i < MAX_ITER && !done; i++) {
        // Schwarzschild geodesic equations
        float f = 1.0 - rs / r;
        float sin_t = sin(theta_ray);
        float cos_t = cos(theta_ray);
        float sin2 = sin_t * sin_t;

        // RK4 step for accuracy
        // k1
        float k1_dr = dr;
        float k1_dtheta = dtheta;
        float k1_dphi = dphi;
        float k1_acc_r = rs/(2.0*r*r) - f*(dtheta*dtheta + sin2*dphi*dphi);
        float k1_acc_theta = -2.0*dr*dtheta/r + sin_t*cos_t*dphi*dphi;
        float k1_acc_phi = -2.0*dr*dphi/r - 2.0*cos_t*dtheta*dphi/max(sin_t, 0.01);

        // k2 (simplified - half step)
        float r2 = r + k1_dr * DL * 0.5;
        float dr2 = dr + k1_acc_r * DL * 0.5;
        float k2_acc_r = rs/(2.0*r2*r2) - (1.0-rs/r2)*(dtheta*dtheta + sin2*dphi*dphi);

        // Average (simplified RK4)
        float acc_r = (k1_acc_r + k2_acc_r) * 0.5;
        float acc_theta = k1_acc_theta;
        float acc_phi = k1_acc_phi;

        // Update velocities
        dr += acc_r * DL;
        dtheta += acc_theta * DL;
        dphi += acc_phi * DL;

        float r_old = r;
        float theta_old = theta_ray;

        // Update positions
        r += dr * DL;
        theta_ray += dtheta * DL;
        phi_ray += dphi * DL;

        // Event horizon
        if(r <= rs * 1.05) {
            color = vec3(0.0);
            done = true;
        }

        // Escaped
        else if(r > MAX_R) {
            vec3 dir = normalize(vec3(
                sin(theta_ray) * cos(phi_ray),
                cos(theta_ray),
                sin(theta_ray) * sin(phi_ray)
            ));
            color = stars(dir);
            done = true;
        }

        // Disk intersection
        else {
            float eq = PI / 2.0;
            if((theta_old - eq) * (theta_ray - eq) < 0.0 && hits < 3) {
                float t = (eq - theta_old) / (theta_ray - theta_old);
                float r_hit = r_old + t * (r - r_old);

                if(r_hit >= r_inner && r_hit <= r_outer) {
                    float dimming = 1.0 / (1.0 + float(hits) * 0.7);
                    diskAcc += diskColor(r_hit, phi_ray) * dimming;
                    hits++;
                }
            }

            // Photon ring glow
            if(r > r_photon * 0.88 && r < r_photon * 1.18) {
                float ringGlow = exp(-12.0 * abs(r - r_photon) / r_photon);
                diskAcc += vec3(1.0, 0.75, 0.4) * ringGlow * 2.5;
            }
        }
    }

    color += diskAcc;

    // Tone mapping
    color = color / (color + vec3(0.55));
    color = pow(color, vec3(0.87));

    // Contrast
    color = (color - 0.5) * 1.18 + 0.5;
    color = clamp(color, 0.0, 1.0);

    // Vignette
    float dist = length(uv);
    float vignette = 1.0 - smoothstep(0.85, 1.7, dist);
    color *= vignette * 0.45 + 0.55;

    FragColor = vec4(color, 1.0);
}
