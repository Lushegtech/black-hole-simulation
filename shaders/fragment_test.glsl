#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;

void main() {
    // Simple test: show a color gradient
    vec2 uv = TexCoord;
    vec3 color = vec3(uv.x, uv.y, 0.5 + 0.5 * sin(u_time));
    FragColor = vec4(color, 1.0);
}
