#ifndef MAT4_H
#define MAT4_H

#include "vec4.h"
#include <array>

/**
 * 4x4 matrix class for metric tensor g_μν and its inverse g^μν
 */
class Mat4 {
public:
    std::array<std::array<double, 4>, 4> data;

    // Constructors
    Mat4() {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                data[i][j] = (i == j) ? 1.0 : 0.0;  // Identity by default
    }

    // Accessors
    std::array<double, 4>& operator[](size_t i) { return data[i]; }
    const std::array<double, 4>& operator[](size_t i) const { return data[i]; }

    double& operator()(size_t i, size_t j) { return data[i][j]; }
    const double& operator()(size_t i, size_t j) const { return data[i][j]; }

    // Matrix-vector multiplication: g_μν * v^ν (lowers index)
    Vec4 operator*(const Vec4& v) const {
        Vec4 result;
        for (int mu = 0; mu < 4; ++mu) {
            result[mu] = 0.0;
            for (int nu = 0; nu < 4; ++nu) {
                result[mu] += data[mu][nu] * v[nu];
            }
        }
        return result;
    }

    // Compute determinant (for metric tensor)
    double determinant() const {
        // Full 4x4 determinant calculation
        double det = 0.0;

        // Using cofactor expansion along first row
        for (int j = 0; j < 4; ++j) {
            det += (j % 2 == 0 ? 1 : -1) * data[0][j] * minor(0, j);
        }

        return det;
    }

    // Compute 3x3 minor for determinant calculation
    double minor(int row, int col) const {
        double m[3][3];
        int mi = 0;

        for (int i = 0; i < 4; ++i) {
            if (i == row) continue;
            int mj = 0;
            for (int j = 0; j < 4; ++j) {
                if (j == col) continue;
                m[mi][mj] = data[i][j];
                ++mj;
            }
            ++mi;
        }

        // 3x3 determinant
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // Compute inverse (using Gauss-Jordan elimination)
    Mat4 inverse() const {
        Mat4 result;
        Mat4 temp = *this;

        // Create augmented matrix [A|I]
        double aug[4][8];
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                aug[i][j] = temp.data[i][j];
                aug[i][j+4] = (i == j) ? 1.0 : 0.0;
            }
        }

        // Forward elimination
        for (int i = 0; i < 4; ++i) {
            // Find pivot
            int maxRow = i;
            for (int k = i + 1; k < 4; ++k) {
                if (std::abs(aug[k][i]) > std::abs(aug[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Swap rows
            if (maxRow != i) {
                for (int j = 0; j < 8; ++j) {
                    std::swap(aug[i][j], aug[maxRow][j]);
                }
            }

            // Scale pivot row
            double pivot = aug[i][i];
            if (std::abs(pivot) < 1e-12) continue;  // Singular matrix

            for (int j = 0; j < 8; ++j) {
                aug[i][j] /= pivot;
            }

            // Eliminate column
            for (int k = 0; k < 4; ++k) {
                if (k != i) {
                    double factor = aug[k][i];
                    for (int j = 0; j < 8; ++j) {
                        aug[k][j] -= factor * aug[i][j];
                    }
                }
            }
        }

        // Extract inverse from right half
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                result.data[i][j] = aug[i][j+4];
            }
        }

        return result;
    }

    // Print for debugging
    friend std::ostream& operator<<(std::ostream& os, const Mat4& m) {
        for (int i = 0; i < 4; ++i) {
            os << "[";
            for (int j = 0; j < 4; ++j) {
                os << m.data[i][j];
                if (j < 3) os << ", ";
            }
            os << "]";
            if (i < 3) os << "\n";
        }
        return os;
    }
};

#endif // MAT4_H
