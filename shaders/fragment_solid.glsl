#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

void main() {
    // Absolute simplest shader - just output solid red
    // If this doesn't work, something is fundamentally broken
    FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}
