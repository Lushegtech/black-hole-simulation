/**
 * Sagittarius A* Black Hole Simulator
 *
 * Physically accurate General Relativistic Ray Tracing using:
 * - Kerr metric (rotating black hole)
 * - Adaptive RKF45 integration
 * - Relativistic beaming and redshift
 * - Novikov-Thorne accretion disk model
 *
 * Based on observations from the Event Horizon Telescope (EHT)
 */

#include "shader.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Simple vector/matrix types (to avoid GLM dependency)
struct vec2 { float x, y; };
struct vec3 { float x, y, z; };
struct vec4 { float x, y, z, w; };
struct mat4 {
    float data[16];
    mat4() { for (int i = 0; i < 16; i++) data[i] = (i % 5 == 0) ? 1.0f : 0.0f; }
};

// ============================================================================
// SAGITTARIUS A* PARAMETERS (from EHT observations)
// ============================================================================

namespace SgrA {
// Black hole properties
constexpr float MASS = 4.15e6;           // Solar masses (Ghez et al. 2008, EHT 2022)
constexpr float SPIN = 0.90f;            // Dimensionless spin a* ∈ [0.0, 0.94] (high spin favored by EHT)
constexpr float DISTANCE = 8.2e3;        // Parsecs from Earth
constexpr float INCLINATION = 17.0f;     // Degrees (face-on to ~30°, model dependent)

// Accretion disk properties
constexpr float DISK_INNER_MULT = 1.0f;  // Inner edge at ISCO
constexpr float DISK_OUTER = 15.0f;      // Outer radius in M units

// Simulation properties
constexpr float OBSERVER_DISTANCE = 100.0f;  // Camera distance in M units
} // namespace SgrA

// ============================================================================
// WINDOW & CAMERA STATE
// ============================================================================

const unsigned int SCR_WIDTH = 1920;
const unsigned int SCR_HEIGHT = 1080;

struct Camera {
    float distance = SgrA::OBSERVER_DISTANCE;
    float theta = M_PI / 2.0f;  // Equatorial
    float phi = 0.0f;
    float fov = 60.0f * M_PI / 180.0f;
} camera;

struct SceneParameters {
    mat4 viewMatrix;
    vec3 cameraPos;
    float blackHoleSpin;
    float blackHoleMass;
    float time;
    float inclination;
    float diskInnerRadiusMult;
    float diskOuterRadius;
    vec2 resolution;
    float fov;
    float _pad1, _pad2;  // Padding for std140 alignment
};

// Mouse state
bool firstMouse = true;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;

// Timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// Adjustable parameters
float spin_parameter = SgrA::SPIN;
float inclination_angle = SgrA::INCLINATION * M_PI / 180.0f;

// ============================================================================
// INPUT HANDLING
// ============================================================================

void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    float camera_speed = 10.0f * deltaTime;
    float angle_speed = 0.5f * deltaTime;

    // Camera distance
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.distance -= camera_speed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.distance += camera_speed;

    // Camera rotation
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.phi -= angle_speed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.phi += angle_speed;
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        camera.theta -= angle_speed;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        camera.theta += angle_speed;

    // Black hole spin adjustment (1-9 keys)
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
        spin_parameter = 0.0f;  // Schwarzschild
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
        spin_parameter = 0.2f;
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
        spin_parameter = 0.4f;
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
        spin_parameter = 0.6f;
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
        spin_parameter = 0.8f;
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
        spin_parameter = 0.90f;  // Sgr A* (nominal)
    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
        spin_parameter = 0.94f;  // Sgr A* (high spin)
    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
        spin_parameter = 0.98f;  // Near-extremal
    if (glfwGetKey(window, GLFW_KEY_9) == GLFW_PRESS)
        spin_parameter = 0.998f; // Extremal Kerr

    // Inclination adjustment (I/K keys)
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS)
        inclination_angle += angle_speed * 0.5f;
    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
        inclination_angle -= angle_speed * 0.5f;

    // Clamp values
    camera.distance = std::max(3.0f, std::min(200.0f, camera.distance));
    camera.theta = std::max(0.1f, std::min((float)M_PI - 0.1f, camera.theta));
    spin_parameter = std::max(0.0f, std::min(0.998f, spin_parameter));
    inclination_angle = std::max(0.0f, std::min((float)M_PI / 2.0f, inclination_angle));
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    (void)window;
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.002f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        camera.phi += xoffset;
        camera.theta += yoffset;
        camera.theta = std::max(0.1f, std::min((float)M_PI - 0.1f, camera.theta));
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    (void)window;
    (void)xoffset;
    camera.fov -= (float)yoffset * 0.05f;
    camera.fov = std::max(10.0f * (float)M_PI / 180.0f,
                          std::min(120.0f * (float)M_PI / 180.0f, camera.fov));
}

// ============================================================================
// COMPUTE SHADER UTILITIES
// ============================================================================

/**
 * Load and compile compute shader
 */
GLuint loadComputeShader(const char* path) {
    std::string code;
    std::ifstream shaderFile;
    shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        shaderFile.open(path);
        std::stringstream shaderStream;
        shaderStream << shaderFile.rdbuf();
        shaderFile.close();
        code = shaderStream.str();
    } catch (std::ifstream::failure& e) {
        std::cerr << "ERROR: Failed to read compute shader file: " << path << std::endl;
        return 0;
    }

    const char* shaderCode = code.c_str();
    GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    glShaderSource(shader, 1, &shaderCode, NULL);
    glCompileShader(shader);

    // Check compilation
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        GLchar infoLog[1024];
        glGetShaderInfoLog(shader, 1024, NULL, infoLog);
        std::cerr << "ERROR: Compute shader compilation failed:\n" << infoLog << std::endl;
        return 0;
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, shader);
    glLinkProgram(program);

    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        GLchar infoLog[1024];
        glGetProgramInfoLog(program, 1024, NULL, infoLog);
        std::cerr << "ERROR: Compute shader linking failed:\n" << infoLog << std::endl;
        return 0;
    }

    glDeleteShader(shader);
    return program;
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Require OpenGL 4.3+ for compute shaders
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT,
                                          "Sagittarius A* - Black Hole Simulator", NULL, NULL);
    if (window == NULL) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSwapInterval(0);  // Disable VSync for max FPS

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    // Print info
    std::cout << "========================================\n";
    std::cout << "  SAGITTARIUS A* BLACK HOLE SIMULATOR\n";
    std::cout << "========================================\n\n";
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << "\n";
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << "\n";
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << "\n\n";

    std::cout << "Astrophysical Parameters (Sgr A*):\n";
    std::cout << "  Mass: " << SgrA::MASS << " M☉\n";
    std::cout << "  Spin (a*): " << SgrA::SPIN << "\n";
    std::cout << "  Distance: " << SgrA::DISTANCE << " pc\n";
    std::cout << "  Inclination: " << SgrA::INCLINATION << "°\n\n";

    std::cout << "Controls:\n";
    std::cout << "  W/S: Camera distance\n";
    std::cout << "  A/D: Orbit left/right\n";
    std::cout << "  Q/E: Move up/down\n";
    std::cout << "  Mouse drag: Free rotation\n";
    std::cout << "  Mouse wheel: Zoom (FOV)\n";
    std::cout << "  1-9: Change black hole spin (1=0.0, 6=0.9, 9=0.998)\n";
    std::cout << "  I/K: Adjust inclination angle\n";
    std::cout << "  ESC: Exit\n\n";

    // Load compute shader
    const char* compute_paths[] = {"shaders/geodesic.comp", "build/shaders/geodesic.comp",
                                   "../shaders/geodesic.comp", nullptr};

    const char* compute_path = nullptr;
    for (int i = 0; compute_paths[i] != nullptr; i++) {
        std::ifstream test(compute_paths[i]);
        if (test.good()) {
            compute_path = compute_paths[i];
            break;
        }
    }

    if (!compute_path) {
        std::cerr << "ERROR: Could not find geodesic.comp shader!\n";
        glfwTerminate();
        return -1;
    }

    std::cout << "Loading compute shader: " << compute_path << "\n";
    GLuint computeShader = loadComputeShader(compute_path);
    if (computeShader == 0) {
        glfwTerminate();
        return -1;
    }
    std::cout << "Compute shader compiled successfully!\n\n";

    // Create output texture (HDR)
    GLuint outputTexture;
    glGenTextures(1, &outputTexture);
    glBindTexture(GL_TEXTURE_2D, outputTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindImageTexture(0, outputTexture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

    // Create UBO for scene parameters
    GLuint uboHandle;
    glGenBuffers(1, &uboHandle);
    glBindBuffer(GL_UNIFORM_BUFFER, uboHandle);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(SceneParameters), NULL, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, uboHandle);

    // Create fullscreen quad for display
    float quadVertices[] = {-1.0f, 1.0f,  0.0f, 1.0f, -1.0f, -1.0f, 0.0f, 0.0f,
                            1.0f,  -1.0f, 1.0f, 0.0f, 1.0f,  1.0f,  1.0f, 1.0f};
    unsigned int quadIndices[] = {0, 1, 2, 0, 2, 3};

    GLuint quadVAO, quadVBO, quadEBO;
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glGenBuffers(1, &quadEBO);

    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, quadEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIndices), quadIndices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // Simple display shader (vertex + fragment)
    const char* vertexShaderSource = R"(
        #version 430 core
        layout (location = 0) in vec2 aPos;
        layout (location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);
            TexCoord = aTexCoord;
        }
    )";

    const char* fragmentShaderSource = R"(
        #version 430 core
        out vec4 FragColor;
        in vec2 TexCoord;
        uniform sampler2D screenTexture;
        void main() {
            FragColor = texture(screenTexture, TexCoord);
        }
    )";

    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);

    GLuint displayShader = glCreateProgram();
    glAttachShader(displayShader, vertexShader);
    glAttachShader(displayShader, fragmentShader);
    glLinkProgram(displayShader);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    std::cout << "Starting real-time rendering...\n\n";

    // Render loop
    int frame_count = 0;
    float fps_timer = 0.0f;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        fps_timer += deltaTime;
        frame_count++;
        if (fps_timer >= 1.0f) {
            std::cout << "FPS: " << frame_count << " | Spin: a*=" << spin_parameter
                      << " | Incl: " << (inclination_angle * 180.0f / M_PI) << "° | "
                      << "Camera: r=" << camera.distance << "    \r" << std::flush;
            fps_timer = 0.0f;
            frame_count = 0;
        }

        processInput(window);

        // Update UBO
        SceneParameters params;
        params.viewMatrix = mat4();  // Identity matrix
        params.cameraPos = {camera.distance, camera.theta, camera.phi};
        params.blackHoleSpin = spin_parameter;
        params.blackHoleMass = 1.0f;  // Geometric units
        params.time = currentFrame;
        params.inclination = inclination_angle;
        params.diskInnerRadiusMult = SgrA::DISK_INNER_MULT;
        params.diskOuterRadius = SgrA::DISK_OUTER;
        params.resolution = {(float)SCR_WIDTH, (float)SCR_HEIGHT};
        params.fov = camera.fov;

        glBindBuffer(GL_UNIFORM_BUFFER, uboHandle);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(SceneParameters), &params);

        // Dispatch compute shader
        glUseProgram(computeShader);
        glDispatchCompute((GLuint)std::ceil(SCR_WIDTH / 8.0),
                         (GLuint)std::ceil(SCR_HEIGHT / 8.0), 1);

        // Wait for compute to finish
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        // Draw result to screen
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(displayShader);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, outputTexture);
        glBindVertexArray(quadVAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        glfwSwapBuffers(window);
    }

    // Cleanup
    glDeleteProgram(computeShader);
    glDeleteProgram(displayShader);
    glDeleteTextures(1, &outputTexture);
    glDeleteBuffers(1, &uboHandle);
    glDeleteVertexArrays(1, &quadVAO);
    glDeleteBuffers(1, &quadVBO);
    glDeleteBuffers(1, &quadEBO);

    glfwTerminate();
    std::cout << "\n\nShutdown complete.\n";
    return 0;
}
