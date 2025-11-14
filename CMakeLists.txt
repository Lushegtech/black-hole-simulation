cmake_minimum_required(VERSION 3.10)
project(BlackHoleRaytracer VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    # Use -O3 for release, but skip -march=native for portability in CI
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Find OpenGL and related libraries (optional for GPU version)
find_package(OpenGL)
find_package(glfw3)
find_package(GLEW)

# CPU-based ray tracer (always build)
set(CPU_SOURCES
    src/main.cpp
)

add_executable(black-hole-raytracer ${CPU_SOURCES})

# Link math library (Unix-like systems only)
if(UNIX)
    target_link_libraries(black-hole-raytracer m)
endif()

# GPU-based real-time renderer (build if OpenGL found)
if(OpenGL_FOUND AND glfw3_FOUND AND GLEW_FOUND)
    message(STATUS "OpenGL found - building GPU renderer")

    set(GPU_SOURCES
        src/main_gpu.cpp
        src/shader.cpp
    )

    add_executable(black-hole-gpu ${GPU_SOURCES})
    target_link_libraries(black-hole-gpu
        OpenGL::GL
        glfw
        GLEW::GLEW
    )

    # Link math library (Unix-like systems only)
    if(UNIX)
        target_link_libraries(black-hole-gpu m)
    endif()
else()
    message(STATUS "OpenGL not found - skipping GPU renderer")
    message(STATUS "  Install: sudo apt-get install libglfw3-dev libglew-dev")
endif()

# Install targets
install(TARGETS black-hole-raytracer DESTINATION bin)
if(TARGET black-hole-gpu)
    install(TARGETS black-hole-gpu DESTINATION bin)
endif()

# Testing
enable_testing()

# Basic smoke test - verify the executable runs and produces output
add_test(NAME BasicRenderTest
    COMMAND black-hole-raytracer test_output.ppm 100 75
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Set a timeout for the test (10 seconds should be plenty)
set_tests_properties(BasicRenderTest PROPERTIES TIMEOUT 10)

# Print configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
