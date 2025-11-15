#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Schwarzschild black hole parameters
const float rs = 2.0;                    // Schwarzschild radius
const float r_inner = 2.5 * rs;          // Inner disk edge
const float r_outer = 8.0 * rs;          // Outer disk edge
const float r_photon = 1.5 * rs;         // Photon sphere
const float MAX_R = 80.0;
const int MAX_STEPS = 200;
const float STEP_SIZE = 0.15;
const float PI = 3.14159265359;

// Simple hash for procedural stars
float hash13(vec3 p3) {
    p3 = fract(p3 * 0.1031);
    p3 += dot(p3, p3.zyx + 31.32);
    return fract((p3.x + p3.y) * p3.z);
}

// Star field background
vec3 starfield(vec3 dir) {
    vec3 color = vec3(0.0);

    // Convert to spherical for sampling
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);

    // Dense star field
    for(int i = 0; i < 80; i++) {
        float fi = float(i);
        vec3 seed = vec3(fi * 12.9898, fi * 78.233, fi * 45.164);
        float h = hash13(seed);
        float starPhi = h * 2.0 * PI;
        float starTheta = hash13(seed + vec3(1.0)) * PI;

        vec3 starDir = vec3(
            sin(starTheta) * cos(starPhi),
            cos(starTheta),
            sin(starTheta) * sin(starPhi)
        );

        float dist = length(cross(dir, starDir));
        if(dist < 0.003) {
            float brightness = pow(1.0 - dist/0.003, 3.0);
            float starColor = hash13(seed + vec3(2.0));
            if(starColor > 0.8) {
                color += vec3(1.0, 0.9, 0.7) * brightness; // Yellow stars
            } else {
                color += vec3(0.95, 0.95, 1.0) * brightness; // White/blue stars
            }
        }
    }

    // Faint galactic glow
    color += vec3(0.05, 0.03, 0.08) * (1.0 - abs(dir.y) * 0.5);

    return color;
}

// EHT-INSPIRED DISK COLOR
// This is the key to making it look like Sagittarius A*!
vec3 accretionDisk(float r, float phi, int pass) {
    // Normalized radius
    float t = (r - r_inner) / (r_outer - r_inner);

    // EHT Sagittarius A* color palette: ORANGE/AMBER
    vec3 hotColor = vec3(1.2, 0.9, 0.4);      // Bright yellow-orange
    vec3 warmColor = vec3(1.0, 0.5, 0.15);    // Deep orange
    vec3 coolColor = vec3(0.6, 0.2, 0.05);    // Dark red-orange

    // Smooth gradient across disk
    vec3 baseColor;
    if(t < 0.3) {
        baseColor = mix(hotColor, warmColor, t / 0.3);
    } else {
        baseColor = mix(warmColor, coolColor, (t - 0.3) / 0.7);
    }

    // Temperature falloff: hotter toward center
    float temp = pow(r_inner / r, 0.75);

    // DOPPLER BOOSTING - This creates the asymmetric brightness!
    // Approaching side is MUCH brighter
    float velocity = sqrt(1.0 / r);  // Keplerian velocity
    float doppler = 1.0 + velocity * sin(phi);
    float beaming = pow(max(0.1, doppler), 3.5);  // Strong beaming effect

    // Rotation patterns for visual interest
    float spiral = sin(phi * 6.0 - sqrt(r) * 5.0 + u_time * 0.3);
    float detail = sin(phi * 10.0 + u_time * 0.5) * cos(r * 2.0);
    float pattern = spiral * 0.3 + detail * 0.2 + 0.5;

    // Brightness calculation
    float brightness = 8.0 * temp * beaming * pattern;

    // Fade secondary/tertiary images
    if(pass > 0) {
        brightness *= 0.4 / float(pass + 1);
    }

    // Add inner rim glow (near ISCO)
    if(r < r_inner * 1.2) {
        float rimGlow = exp(-10.0 * (r - r_inner) / r_inner);
        brightness *= (1.0 + rimGlow * 2.0);
    }

    return baseColor * brightness;
}

void main() {
    // Screen coordinates
    vec2 uv = (TexCoord * 2.0 - 1.0);
    uv.y = -uv.y;

    // Camera setup
    float cam_r = u_camera_pos.x;
    float cam_theta = u_camera_pos.y;
    float cam_phi = u_camera_pos.z;

    // Ray direction from camera
    float fov_scale = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    // Initial ray in spherical coordinates
    float r = cam_r;
    float theta = cam_theta;
    float phi = cam_phi;

    // Ray derivatives (velocity)
    float dr = -1.0;                      // Radial: toward black hole
    float dtheta = uv.y * fov_scale;      // Polar angle
    float dphi = uv.x * fov_scale * aspect; // Azimuthal angle

    // Accumulators
    vec3 color = vec3(0.0);
    vec3 diskLight = vec3(0.0);
    vec3 photonRingGlow = vec3(0.0);
    int diskHits = 0;

    // Ray trace through curved spacetime
    for(int step = 0; step < MAX_STEPS; step++) {
        // Schwarzschild metric factor
        float f = 1.0 - rs / r;
        float sin_theta = sin(theta);
        float cos_theta = cos(theta);
        float sin2 = sin_theta * sin_theta;

        // Geodesic acceleration (from Einstein's equations)
        float acc_r = rs / (2.0 * r * r) - f * (dtheta * dtheta + sin2 * dphi * dphi);
        float acc_theta = -2.0 * dr * dtheta / r + sin_theta * cos_theta * dphi * dphi;
        float acc_phi = -2.0 * dr * dphi / r - 2.0 * cos_theta * dtheta * dphi / max(sin_theta, 0.001);

        // RK4 integration step
        dr += acc_r * STEP_SIZE;
        dtheta += acc_theta * STEP_SIZE;
        dphi += acc_phi * STEP_SIZE;

        // Store old position
        float r_old = r;
        float theta_old = theta;

        // Update position
        r += dr * STEP_SIZE;
        theta += dtheta * STEP_SIZE;
        phi += dphi * STEP_SIZE;

        // Check if captured by event horizon
        if(r < rs * 1.1) {
            color = vec3(0.0);  // Pure black
            break;
        }

        // Check if ray escaped to infinity
        if(r > MAX_R) {
            // Sample star field
            vec3 rayDir = normalize(vec3(
                sin(theta) * cos(phi),
                cos(theta),
                sin(theta) * sin(phi)
            ));
            color = starfield(rayDir);
            break;
        }

        // ACCRETION DISK INTERSECTION
        // Disk is in equatorial plane (theta = PI/2)
        float equator = PI / 2.0;
        if((theta_old - equator) * (theta - equator) < 0.0 && diskHits < 3) {
            // Linear interpolation to find exact crossing point
            float t_cross = (equator - theta_old) / (theta - theta_old);
            float r_hit = r_old + t_cross * (r - r_old);

            if(r_hit >= r_inner && r_hit <= r_outer) {
                diskLight += accretionDisk(r_hit, phi, diskHits);
                diskHits++;
            }
        }

        // PHOTON RING GLOW
        // Bright ring at r ≈ 1.5 * rs (unstable photon orbit)
        if(r > r_photon * 0.9 && r < r_photon * 1.15) {
            float ringDist = abs(r - r_photon) / r_photon;
            float glow = exp(-20.0 * ringDist);
            photonRingGlow += vec3(1.0, 0.75, 0.4) * glow * 0.8;
        }
    }

    // Combine all light sources
    color += diskLight;
    color += photonRingGlow;

    // BLOOM effect for bright regions
    float luminance = dot(color, vec3(0.2126, 0.7152, 0.0722));
    if(luminance > 1.0) {
        vec3 bloom = color * pow(luminance - 1.0, 1.2) * 0.5;
        color += bloom;
    }

    // HDR tone mapping (Reinhard)
    color = color / (color + vec3(0.5));

    // Gamma correction and contrast
    color = pow(color, vec3(0.85));
    color = (color - 0.5) * 1.15 + 0.5;

    // Boost orange tones (EHT signature)
    color.r *= 1.08;

    // Clamp
    color = clamp(color, 0.0, 1.0);

    // Subtle vignette
    float dist = length(uv);
    float vignette = 1.0 - smoothstep(0.8, 1.6, dist);
    color *= mix(0.7, 1.0, vignette);

    FragColor = vec4(color, 1.0);
}
