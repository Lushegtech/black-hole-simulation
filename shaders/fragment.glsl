#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

// Schwarzschild black hole
const float rs = 2.0;                    // Schwarzschild radius
const float r_isco = 3.0 * rs;           // ISCO for Schwarzschild
const float r_inner = 2.8 * rs;          // Inner disk edge
const float r_outer = 7.0 * rs;          // Outer disk edge
const float r_photon = 1.5 * rs;         // Photon sphere
const float MAX_R = 100.0;
const int MAX_STEPS = 250;
const float STEP_SIZE = 0.12;
const float PI = 3.14159265359;

// Procedural noise for subtle texture
float hash(vec2 p) {
    p = fract(p * vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x * p.y);
}

// Star field
vec3 starfield(vec3 dir) {
    vec3 color = vec3(0.0);

    // Simplified procedural stars
    vec2 starCoord = vec2(atan(dir.z, dir.x), acos(clamp(dir.y, -1.0, 1.0)));
    starCoord *= 100.0;

    for(int i = 0; i < 60; i++) {
        vec2 offset = vec2(float(i) * 12.34, float(i) * 56.78);
        vec2 cellCoord = floor(starCoord + offset);
        float h = hash(cellCoord);

        if(h > 0.98) {
            vec2 starPos = cellCoord + vec2(hash(cellCoord + vec2(1.0)), hash(cellCoord + vec2(2.0)));
            float dist = length(starCoord + offset - starPos);
            if(dist < 0.5) {
                float brightness = pow(1.0 - dist/0.5, 4.0) * 0.8;
                color += vec3(1.0, 0.95, 0.9) * brightness;
            }
        }
    }

    // Subtle background glow
    color += vec3(0.02, 0.015, 0.025);

    return color;
}

// EHT SAGITTARIUS A* ACCRETION DISK
vec3 diskEmission(float r, float phi, int pass) {
    // Distance from inner edge
    float t = (r - r_inner) / (r_outer - r_inner);

    // EHT Color Palette: ORANGE/AMBER gradient
    vec3 innerColor = vec3(1.5, 1.0, 0.5);     // Bright orange-yellow
    vec3 midColor = vec3(1.2, 0.6, 0.2);       // Deep orange
    vec3 outerColor = vec3(0.8, 0.3, 0.1);     // Dark red-orange

    vec3 diskColor;
    if(t < 0.4) {
        diskColor = mix(innerColor, midColor, t / 0.4);
    } else {
        diskColor = mix(midColor, outerColor, (t - 0.4) / 0.6);
    }

    // Temperature profile: T ∝ r^(-3/4)
    float temperature = pow(r_inner / r, 0.75);

    // DOPPLER BOOSTING - key to asymmetry!
    // Disk rotating counter-clockwise as seen from above
    float omega = sqrt(1.0 / (r * r * r));  // Keplerian angular velocity
    float velocity = r * omega;              // Linear velocity

    // Doppler factor: approaching side is brighter
    float cosPhi = cos(phi);
    float doppler = 1.0 + velocity * cosPhi;

    // Relativistic beaming: I' = I * doppler^(3+α), α ~ 0 for thermal emission
    float beaming = pow(max(0.05, doppler), 4.0);

    // Add turbulent structure (not too much - EHT is relatively smooth)
    float turbulence = 0.5 + 0.5 * sin(phi * 5.0 + r * 3.0 - u_time * 0.2);
    turbulence *= 0.5 + 0.5 * sin(phi * 11.0 - r * 7.0 + u_time * 0.15);
    turbulence = 0.6 + turbulence * 0.4;  // Keep it subtle

    // Base intensity
    float intensity = temperature * beaming * turbulence;

    // Brightness multiplier
    intensity *= 12.0;

    // Fade secondary images (gravitationally lensed paths)
    if(pass > 0) {
        intensity *= 0.3 / float(pass + 1);
    }

    // Hot inner edge glow
    if(r < r_inner * 1.15) {
        float edgeGlow = exp(-8.0 * (r - r_inner) / r_inner);
        intensity *= (1.0 + edgeGlow * 3.0);
    }

    return diskColor * intensity;
}

void main() {
    // Normalized screen coordinates
    vec2 uv = (TexCoord * 2.0 - 1.0);
    uv.y = -uv.y;

    // Camera position in spherical coordinates
    float cam_r = u_camera_pos.x;
    float cam_theta = u_camera_pos.y;
    float cam_phi = u_camera_pos.z;

    // FOV and aspect ratio
    float fov_scale = tan(u_fov * 0.5);
    float aspect = u_resolution.x / u_resolution.y;

    // Initialize ray in spherical coordinates
    float r = cam_r;
    float theta = cam_theta;
    float phi = cam_phi;

    // Ray direction (velocity in phase space)
    float dr = -1.0;                          // Radial velocity (toward BH)
    float dtheta = uv.y * fov_scale;          // Polar angle velocity
    float dphi = uv.x * fov_scale * aspect;   // Azimuthal velocity

    // Conserved quantities for better integration
    float E = sqrt(1.0 - rs / r);             // Energy at infinity
    float L = r * r * dphi;                   // Angular momentum

    // Light accumulators
    vec3 finalColor = vec3(0.0);
    vec3 accumulatedDisk = vec3(0.0);
    vec3 photonRing = vec3(0.0);
    int diskCrossings = 0;
    bool hitHorizon = false;

    // GEODESIC INTEGRATION
    for(int step = 0; step < MAX_STEPS; step++) {
        // Metric function
        float f = 1.0 - rs / r;

        // Trig values
        float sin_theta = sin(theta);
        float cos_theta = cos(theta);
        float sin2_theta = sin_theta * sin_theta;

        if(sin_theta < 0.001) sin_theta = 0.001; // Avoid singularity

        // Geodesic equations in Schwarzschild spacetime
        // d²r/dλ²
        float d2r = -rs / (2.0 * r * r) * f + f * (dtheta * dtheta + sin2_theta * dphi * dphi);

        // d²θ/dλ²
        float d2theta = -2.0 * dr * dtheta / r + sin_theta * cos_theta * dphi * dphi;

        // d²φ/dλ²
        float d2phi = -2.0 * dr * dphi / r - 2.0 * cos_theta * dtheta * dphi / sin_theta;

        // Save old position
        float r_old = r;
        float theta_old = theta;

        // Euler step (simple but fast)
        dr += d2r * STEP_SIZE;
        dtheta += d2theta * STEP_SIZE;
        dphi += d2phi * STEP_SIZE;

        r += dr * STEP_SIZE;
        theta += dtheta * STEP_SIZE;
        phi += dphi * STEP_SIZE;

        // Check if ray hit event horizon
        if(r < rs * 1.05) {
            hitHorizon = true;
            break;
        }

        // Check if ray escaped to infinity
        if(r > MAX_R) {
            // Sample background stars
            vec3 rayDir = vec3(
                sin(theta) * cos(phi),
                cos(theta),
                sin(theta) * sin(phi)
            );
            finalColor = starfield(rayDir);
            break;
        }

        // ACCRETION DISK INTERSECTION
        // Disk lies in equatorial plane (theta = PI/2)
        float equator = PI / 2.0;
        float crossing = (theta_old - equator) * (theta - equator);

        if(crossing < 0.0 && diskCrossings < 3) {
            // Interpolate to find exact crossing point
            float t = (equator - theta_old) / (theta - theta_old);
            float r_cross = r_old + t * (r - r_old);

            // Check if within disk radii
            if(r_cross >= r_inner && r_cross <= r_outer) {
                vec3 emission = diskEmission(r_cross, phi, diskCrossings);
                accumulatedDisk += emission;
                diskCrossings++;
            }
        }

        // PHOTON RING visualization
        // Bright ring at photon sphere radius
        if(r > r_photon * 0.95 && r < r_photon * 1.1) {
            float ringDist = abs(r - r_photon) / (r_photon * 0.1);
            float ringGlow = exp(-15.0 * ringDist);
            photonRing += vec3(1.2, 0.8, 0.4) * ringGlow * 0.6;
        }
    }

    // Combine all components
    vec3 color = finalColor;

    // Add disk emission
    color += accumulatedDisk;

    // Add photon ring
    color += photonRing;

    // If hit horizon, ensure it's pure black
    if(hitHorizon) {
        color = vec3(0.0);
    }

    // POST-PROCESSING for EHT look

    // Strong bloom on bright areas
    float lum = dot(color, vec3(0.2126, 0.7152, 0.0722));
    if(lum > 1.5) {
        vec3 bloom = color * pow(lum - 1.5, 1.5) * 0.6;
        color += bloom;
    }

    // Tone mapping: Reinhard with slight adjustment
    color = color / (color + vec3(0.8));

    // Gamma correction
    color = pow(color, vec3(0.9));

    // Boost contrast
    color = (color - 0.5) * 1.2 + 0.5;

    // Emphasize orange (EHT signature color)
    color.r *= 1.1;
    color.g *= 0.95;

    // Clamp
    color = clamp(color, 0.0, 1.0);

    // Vignette
    float vignette = 1.0 - 0.3 * pow(length(uv), 2.0);
    color *= vignette;

    FragColor = vec4(color, 1.0);
}
