#!/bin/bash

# Build script for Black Hole Simulation

set -e  # Exit on error

echo "========================================="
echo " Black Hole Simulation - Build Script"
echo "========================================="
echo ""

# Check for dependencies
echo "Checking dependencies..."

# Check for compiler
if ! command -v g++ &> /dev/null; then
    echo "Error: g++ not found. Please install build-essential."
    echo "  sudo apt-get install build-essential"
    exit 1
fi

# Check for CMake
if ! command -v cmake &> /dev/null; then
    echo "Error: cmake not found. Please install CMake."
    echo "  sudo apt-get install cmake"
    exit 1
fi

echo "✓ Compiler found: $(g++ --version | head -n1)"
echo "✓ CMake found: $(cmake --version | head -n1)"

# Check for OpenGL libraries (optional)
echo ""
echo "Checking optional GPU dependencies..."

if pkg-config --exists glfw3; then
    echo "✓ GLFW3 found"
    HAS_GLFW=true
else
    echo "✗ GLFW3 not found (GPU renderer will be skipped)"
    HAS_GLFW=false
fi

if pkg-config --exists glew; then
    echo "✓ GLEW found"
    HAS_GLEW=true
else
    echo "✗ GLEW not found (GPU renderer will be skipped)"
    HAS_GLEW=false
fi

if [ "$HAS_GLFW" = false ] || [ "$HAS_GLEW" = false ]; then
    echo ""
    echo "To install GPU dependencies:"
    echo "  sudo apt-get install libglfw3-dev libglew-dev libgl1-mesa-dev"
fi

# Create build directory
echo ""
echo "Creating build directory..."
mkdir -p build
mkdir -p output

# Configure with CMake
echo ""
echo "Configuring with CMake..."
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build
echo ""
echo "Building..."
cmake --build . -j$(nproc)

# Summary
echo ""
echo "========================================="
echo " Build Complete!"
echo "========================================="
echo ""

if [ -f "black-hole-raytracer" ]; then
    echo "✓ CPU Ray Tracer: ./build/black-hole-raytracer"
fi

if [ -f "black-hole-gpu" ]; then
    echo "✓ GPU Real-Time Renderer: ./build/black-hole-gpu"
fi

echo ""
echo "Run examples:"
echo "  ./build/black-hole-raytracer output/test.ppm 400 300"
if [ -f "black-hole-gpu" ]; then
    echo "  ./build/black-hole-gpu"
fi

echo ""
