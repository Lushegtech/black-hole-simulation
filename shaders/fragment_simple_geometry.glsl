#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float PI = 3.14159265359;

// Simple hash for stars
float hash12(vec2 p) {
    vec3 p3  = fract(vec3(p.xyx) * 0.1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

void main() {
    vec2 uv = (2.0 * TexCoord - 1.0);
    uv.y = -uv.y;
    uv.x *= u_resolution.x / u_resolution.y;

    float dist = length(uv);
    vec3 color = vec3(0.0);

    // Stars background
    for (int i = 0; i < 100; i++) {
        vec2 starPos = vec2(hash12(vec2(i, 0)), hash12(vec2(i, 1)));
        starPos = starPos * 4.0 - 2.0;
        float d = length(uv - starPos);
        if (d < 0.01) {
            float brightness = 1.0 - (d / 0.01);
            color += vec3(1.0, 0.9, 0.8) * brightness * 0.5;
        }
    }

    // Accretion disk (simple ring)
    float diskDist = abs(dist - 0.4);
    if (diskDist < 0.15) {
        float diskBrightness = 1.0 - (diskDist / 0.15);

        // Rotating effect
        float angle = atan(uv.y, uv.x);
        float rotation = angle + u_time * 0.5;

        // Color gradient
        vec3 diskColor = vec3(1.0, 0.6, 0.2);
        diskColor *= diskBrightness * (0.7 + 0.3 * sin(rotation * 5.0));

        color += diskColor;
    }

    // Black hole shadow
    if (dist < 0.25) {
        color = vec3(0.0);
    }

    // Vignette
    float vignette = 1.0 - smoothstep(0.5, 1.5, dist);
    color *= vignette * 0.5 + 0.5;

    FragColor = vec4(color, 1.0);
}
