/**
 * Black Hole GPU Renderer - Real-time visualization
 *
 * GPU-accelerated General Relativistic Ray Tracing using OpenGL/GLSL
 * This implementation uses fragment shaders to trace geodesics in parallel
 * for real-time rendering of a Schwarzschild black hole.
 */

#include "shader.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <cmath>
#include <iostream>

// Window settings
const unsigned int SCR_WIDTH = 1280;
const unsigned int SCR_HEIGHT = 720;

// Camera settings
float camera_distance = 20.0f;
float camera_theta = M_PI / 2.0f; // Equatorial plane
float camera_phi = 0.0f;
float camera_fov = 60.0f * M_PI / 180.0f;

// Mouse state
bool firstMouse = true;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;

// Timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

/**
 * Process keyboard input
 */
void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    float camera_speed = 5.0f * deltaTime;

    // Camera distance controls
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera_distance -= camera_speed;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera_distance += camera_speed;

    // Camera angle controls
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera_phi -= camera_speed * 0.1f;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera_phi += camera_speed * 0.1f;

    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        camera_theta -= camera_speed * 0.1f;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
        camera_theta += camera_speed * 0.1f;

    // Clamp values
    camera_distance = std::max(3.0f, std::min(100.0f, camera_distance));
    camera_theta = std::max(0.1f, std::min((float)M_PI - 0.1f, camera_theta));
}

/**
 * Framebuffer size callback
 */
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    (void)window;  // Unused parameter (required by GLFW callback signature)
    glViewport(0, 0, width, height);
}

/**
 * Mouse callback for camera control
 */
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // Reversed: y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.001f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    // Only rotate if mouse button is pressed
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        camera_phi += xoffset;
        camera_theta += yoffset;
        camera_theta = std::max(0.1f, std::min((float)M_PI - 0.1f, camera_theta));
    }
}

/**
 * Scroll callback for FOV control
 */
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    (void)window;   // Unused parameter
    (void)xoffset;  // Unused parameter
    camera_fov -= (float)yoffset * 0.05f;
    camera_fov = std::max(10.0f * (float)M_PI / 180.0f,
                          std::min(120.0f * (float)M_PI / 180.0f, camera_fov));
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create window
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT,
                                          "Black Hole Simulation - GPU", NULL, NULL);
    if (window == NULL) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        return -1;
    }

    // Print OpenGL info
    std::cout << "========================================\n";
    std::cout << "  Black Hole GPU Renderer\n";
    std::cout << "========================================\n\n";
    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << "\n";
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << "\n";
    std::cout << "Renderer: " << glGetString(GL_RENDERER) << "\n\n";

    std::cout << "Controls:\n";
    std::cout << "  W/S: Move camera closer/farther\n";
    std::cout << "  A/D: Rotate around black hole\n";
    std::cout << "  Q/E: Move up/down\n";
    std::cout << "  Mouse drag: Rotate view\n";
    std::cout << "  Mouse wheel: Zoom (FOV)\n";
    std::cout << "  ESC: Exit\n\n";

    // Build and compile shader program
    // Shaders are copied to build directory by CMake
    Shader shader("shaders/vertex.glsl", "shaders/fragment.glsl");

    // Set up vertex data for fullscreen quad
    float vertices[] = {// positions        // texture coords
                        -1.0f, 1.0f,  0.0f, 0.0f, 1.0f, -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
                        1.0f,  -1.0f, 0.0f, 1.0f, 0.0f, 1.0f,  1.0f,  0.0f, 1.0f, 1.0f};

    unsigned int indices[] = {0, 1, 2, 0, 2, 3};

    // Create VAO, VBO, EBO
    unsigned int VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // Texture coord attribute
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float),
                          (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    std::cout << "Starting real-time rendering...\n";

    // Render loop
    int frame_count = 0;
    float fps_timer = 0.0f;

    while (!glfwWindowShouldClose(window)) {
        // Per-frame time logic
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // FPS counter
        fps_timer += deltaTime;
        frame_count++;
        if (fps_timer >= 1.0f) {
            std::cout << "FPS: " << frame_count << " | Camera: r=" << camera_distance
                      << ", θ=" << camera_theta << ", φ=" << camera_phi << "    \r"
                      << std::flush;
            fps_timer = 0.0f;
            frame_count = 0;
        }

        // Input
        processInput(window);

        // Render
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Use shader and set uniforms
        shader.use();
        shader.setFloat("u_time", currentFrame);
        shader.setVec2("u_resolution", SCR_WIDTH, SCR_HEIGHT);
        shader.setVec3("u_camera_pos", camera_distance, camera_theta, camera_phi);
        shader.setFloat("u_fov", camera_fov);

        // Draw fullscreen quad
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);

    glfwTerminate();

    std::cout << "\n\nShutdown complete.\n";
    return 0;
}
