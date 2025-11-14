cmake_minimum_required(VERSION 3.10)
project(BlackHoleRaytracer VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/main.cpp
    src/raytracer.cpp
    src/vector3d.cpp
)

# Create executable
add_executable(black-hole-raytracer ${SOURCES})

# Install targets
install(TARGETS black-hole-raytracer DESTINATION bin)

# Testing
enable_testing()

# Optional: Add tests if test files exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vector3d.cpp")
    add_executable(test_vector3d tests/test_vector3d.cpp src/vector3d.cpp)
    add_test(NAME Vector3DTest COMMAND test_vector3d)
endif()

# Print configuration
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
