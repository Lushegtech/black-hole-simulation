#ifndef VEC4_H
#define VEC4_H

#include <cmath>
#include <array>
#include <iostream>

/**
 * 4-component vector class for spacetime coordinates and velocities
 * Components: (t, x, y, z) or (t, r, theta, phi) depending on coordinate system
 */
class Vec4 {
public:
    std::array<double, 4> data;

    // Constructors
    Vec4() : data{0.0, 0.0, 0.0, 0.0} {}
    Vec4(double t, double x, double y, double z) : data{t, x, y, z} {}
    explicit Vec4(const std::array<double, 4>& arr) : data(arr) {}

    // Accessors
    double& operator[](size_t i) { return data[i]; }
    const double& operator[](size_t i) const { return data[i]; }

    double& t() { return data[0]; }
    double& x() { return data[1]; }
    double& y() { return data[2]; }
    double& z() { return data[3]; }

    const double& t() const { return data[0]; }
    const double& x() const { return data[1]; }
    const double& y() const { return data[2]; }
    const double& z() const { return data[3]; }

    // Vector operations
    Vec4 operator+(const Vec4& other) const {
        return Vec4(data[0] + other.data[0], data[1] + other.data[1],
                   data[2] + other.data[2], data[3] + other.data[3]);
    }

    Vec4 operator-(const Vec4& other) const {
        return Vec4(data[0] - other.data[0], data[1] - other.data[1],
                   data[2] - other.data[2], data[3] - other.data[3]);
    }

    Vec4 operator*(double scalar) const {
        return Vec4(data[0] * scalar, data[1] * scalar,
                   data[2] * scalar, data[3] * scalar);
    }

    Vec4 operator/(double scalar) const {
        return Vec4(data[0] / scalar, data[1] / scalar,
                   data[2] / scalar, data[3] / scalar);
    }

    Vec4& operator+=(const Vec4& other) {
        for (size_t i = 0; i < 4; ++i) data[i] += other.data[i];
        return *this;
    }

    // Euclidean norm (for debugging/testing, not physical in spacetime)
    double norm() const {
        return std::sqrt(data[0]*data[0] + data[1]*data[1] +
                        data[2]*data[2] + data[3]*data[3]);
    }

    // Print for debugging
    friend std::ostream& operator<<(std::ostream& os, const Vec4& v) {
        os << "(" << v.data[0] << ", " << v.data[1] << ", "
           << v.data[2] << ", " << v.data[3] << ")";
        return os;
    }
};

// Scalar multiplication from left
inline Vec4 operator*(double scalar, const Vec4& v) {
    return v * scalar;
}

#endif // VEC4_H
