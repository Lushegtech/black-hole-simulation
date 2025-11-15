#version 330 core

out vec4 FragColor;
in vec2 TexCoord;

uniform float u_time;
uniform vec2 u_resolution;
uniform vec3 u_camera_pos;
uniform float u_fov;

const float r_s = 2.0;
const float r_disk_inner = 6.0;
const float r_disk_outer = 20.0;
const float PI = 3.14159265359;

void main() {
    // Camera setup
    float r_cam = u_camera_pos.x;
    float theta_cam = u_camera_pos.y;
    float phi_cam = u_camera_pos.z;

    // Ray direction (simplified - no geodesics yet)
    vec2 pixel_coord = (2.0 * TexCoord - 1.0);
    pixel_coord.y = -pixel_coord.y;

    float aspect = u_resolution.x / u_resolution.y;
    vec3 forward = vec3(0, 0, -1);
    vec3 right = vec3(aspect, 0, 0);
    vec3 up = vec3(0, 1, 0);

    vec3 ray_dir = normalize(forward + pixel_coord.x * right * 0.5 + pixel_coord.y * up * 0.5);

    // Simple sphere test - show orange circle where black hole should be
    vec3 cam_pos = vec3(0, 0, r_cam);

    // Ray-sphere intersection for black hole shadow
    vec3 oc = cam_pos;
    float a = dot(ray_dir, ray_dir);
    float b = 2.0 * dot(oc, ray_dir);
    float c = dot(oc, oc) - r_s * r_s * 4.0; // Make it visible
    float discriminant = b*b - 4.0*a*c;

    vec3 color = vec3(0.0);

    if (discriminant > 0.0) {
        // Hit the black hole region - show orange
        color = vec3(1.0, 0.5, 0.0);
    } else {
        // Miss - show blue background with gradient
        color = vec3(0.1, 0.2, 0.4) + vec3(0.0, 0.0, pixel_coord.y * 0.2);
    }

    // Add some animation to verify it's updating
    color += vec3(0.1) * sin(u_time);

    FragColor = vec4(color, 1.0);
}
