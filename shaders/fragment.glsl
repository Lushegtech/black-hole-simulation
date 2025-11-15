#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float rs = 2.0;
const float r_inner = 2.2 * rs;
const float r_outer = 6.0 * rs;
const float r_photon = 1.5 * rs;
const float MAX_R = 100.0;
const int MAX_ITER = 140;
const float DL = 0.2;
const float PI = 3.14159265359;

// Better noise for dust/turbulence
float hash(vec2 p) {
    p = fract(p * 0.3183099 + 0.1);
    p *= 17.0;
    return fract(p.x * p.y * (p.x + p.y));
}

float noise(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));
    return mix(mix(a, b, f.x), mix(c, d, f.x), f.y);
}

// Fractal noise for turbulent dust
float fbm(vec2 p) {
    float value = 0.0;
    float amplitude = 0.5;
    for(int i = 0; i < 4; i++) {
        value += amplitude * noise(p);
        p *= 2.0;
        amplitude *= 0.5;
    }
    return value;
}

vec3 stars(vec3 dir) {
    float theta = acos(clamp(dir.y, -1.0, 1.0));
    float phi = atan(dir.z, dir.x);
    vec2 uv = vec2(phi / (2.0 * PI), theta / PI);

    vec3 col = vec3(0.0);

    for(int i = 0; i < 40; i++) {
        float fi = float(i);
        vec2 p = vec2(hash(vec2(fi * 0.1, 0.0)), hash(vec2(fi * 0.1, 1.0)));
        float d = length(uv - p);
        if(d < 0.004) {
            float bright = pow(1.0 - d/0.004, 2.0);
            col += vec3(1.0, 0.95, 0.9) * bright;
        }
    }

    return col;
}

// STUNNING disk with glow, dust, and atmosphere
vec3 diskColor(float r, float phi, int passNum) {
    // Normalized distance
    float t = (r - r_inner) / (r_outer - r_inner);

    // EHT-STYLE ORANGE GLOW
    // Sagittarius A* appears orange/amber in EHT images
    vec3 innerColor = vec3(1.0, 0.85, 0.5);  // Bright yellow-orange
    vec3 midColor = vec3(1.0, 0.6, 0.2);      // Orange
    vec3 outerColor = vec3(0.9, 0.3, 0.1);    // Deep orange-red

    vec3 color;
    if(t < 0.4) {
        color = mix(innerColor, midColor, t / 0.4);
    } else {
        color = mix(midColor, outerColor, (t - 0.4) / 0.6);
    }

    // TURBULENT DUST/GAS
    vec2 turbCoord = vec2(phi * 3.0, log(r) * 2.0);
    float turbulence = fbm(turbCoord + u_time * 0.1);

    // Add dust variation to color
    color = mix(color, color * 1.3, turbulence * 0.4);

    // BRIGHTNESS with realistic falloff
    float brightness = 12.0 * pow(r_inner / r, 1.8);
    brightness = clamp(brightness, 0.0, 20.0);

    // ANIMATED ROTATION with turbulence
    float rotation = sin(phi * 8.0 - sqrt(r) * 4.0 + u_time * 0.4) * 0.4 + 0.6;
    float swirl = sin(phi * 12.0 - log(r) * 5.0 + u_time * 0.3) * 0.3 + 0.7;

    // DRAMATIC DOPPLER BOOSTING
    // This is what makes one side SUPER bright in EHT image!
    float v = sqrt(1.0 / r);
    float doppler = 1.0 + v * sin(phi) * 0.9;
    float beaming = pow(max(0.15, doppler), 4.5);  // VERY strong!

    // Combine all effects
    brightness *= rotation * swirl * beaming * turbulence;

    // Add GLOW for dusty atmosphere
    float glow = exp(-2.0 * t);  // Brighter inner regions glow more
    brightness *= (1.0 + glow * 0.8);

    // Dimming for secondary/tertiary passes (lensed images)
    if(passNum > 0) {
        brightness *= 0.5 / float(passNum + 1);
    }

    return color * brightness;
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

    float theta_ray = theta;
    float phi_ray = phi_cam;

    vec3 color = vec3(0.0);
    vec3 diskAcc = vec3(0.0);
    int diskHits = 0;
    float glowAcc = 0.0;

    for(int i = 0; i < MAX_ITER; i++) {
        float f = 1.0 - rs / r;
        float sin_t = sin(theta_ray);
        float cos_t = cos(theta_ray);
        float sin2 = sin_t * sin_t;

        // Geodesic
        float acc_r = rs/(2.0*r*r) - f*(dtheta*dtheta + sin2*dphi*dphi);
        float acc_theta = -2.0*dr*dtheta/r + sin_t*cos_t*dphi*dphi;
        float acc_phi = -2.0*dr*dphi/r - 2.0*cos_t*dtheta*dphi/max(sin_t, 0.01);

        dr += acc_r * DL;
        dtheta += acc_theta * DL;
        dphi += acc_phi * DL;

        float r_old = r;
        float theta_old = theta_ray;

        r += dr * DL;
        theta_ray += dtheta * DL;
        phi_ray += dphi * DL;

        // Event horizon
        if(r <= rs * 1.05) {
            break;
        }

        // Escaped
        if(r > MAX_R) {
            vec3 dir = normalize(vec3(
                sin(theta_ray) * cos(phi_ray),
                cos(theta_ray),
                sin(theta_ray) * sin(phi_ray)
            ));
            color = stars(dir);
            break;
        }

        // Disk intersection with MULTIPLE passes
        float eq = PI / 2.0;
        if((theta_old - eq) * (theta_ray - eq) < 0.0 && diskHits < 4) {
            float t = (eq - theta_old) / (theta_ray - theta_old);
            float r_hit = r_old + t * (r - r_old);

            if(r_hit >= r_inner && r_hit <= r_outer) {
                diskAcc += diskColor(r_hit, phi_ray, diskHits);
                diskHits++;
            }
        }

        // ATMOSPHERIC GLOW near disk
        if(r > r_inner * 0.9 && r < r_outer * 1.3) {
            float dist_to_eq = abs(theta_ray - PI/2.0);
            if(dist_to_eq < 0.3) {
                float glow = exp(-dist_to_eq * 8.0) * exp(-abs(r - r_inner*1.5)/r_outer);
                glowAcc += glow * 0.15;
            }
        }

        // PHOTON RING GLOW
        if(r > r_photon * 0.85 && r < r_photon * 1.2) {
            float ringGlow = exp(-15.0 * abs(r - r_photon) / r_photon);
            diskAcc += vec3(1.0, 0.8, 0.5) * ringGlow * 3.0;
        }
    }

    color += diskAcc;

    // Add atmospheric glow (dusty halo)
    color += vec3(1.0, 0.65, 0.3) * glowAcc;

    // BLOOM effect for bright areas
    float avgBright = (color.r + color.g + color.b) / 3.0;
    if(avgBright > 1.0) {
        vec3 bloom = color * pow(avgBright - 1.0, 1.5) * 0.4;
        color += bloom;
    }

    // HDR tone mapping
    color = color / (color + vec3(0.4));

    // Enhanced contrast and saturation
    color = pow(color, vec3(0.85));
    color = (color - 0.5) * 1.25 + 0.5;

    // Boost orange tones (EHT style)
    color.r *= 1.1;
    color.g *= 1.05;

    color = clamp(color, 0.0, 1.0);

    // Vignette
    float dist = length(uv);
    float vignette = 1.0 - smoothstep(0.7, 1.8, dist);
    color *= vignette * 0.4 + 0.6;

    FragColor = vec4(color, 1.0);
}
