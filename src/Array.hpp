#pragma once
#include<array>

using Array2 = std::array<double, 2>;
using Array2x2 = std::array<std::array<double, 2>, 2>;

inline Array2& operator+=(Array2& a, const Array2& b) {
    a[0] += b[0];
    a[1] += b[1];
    return a;
}
inline Array2& operator-=(Array2& a, const Array2& b) {
    a[0] -= b[0];
    a[1] -= b[1];
    return a;
}
inline Array2 operator+(Array2 a, const Array2& b) {
    a += b;
    return a;
}
inline Array2 operator-(Array2 a, const Array2& b) {
    a -= b;
    return a;
}
inline Array2x2& operator+=(Array2x2& a, const Array2x2& b) {
    a[0][0] += b[0][0];
    a[0][1] += b[0][1];
    a[1][0] += b[1][0];
    a[1][1] += b[1][1];
    return a;
}
inline Array2x2 operator+(Array2x2 a, const Array2x2& b) {
    a += b;
    return a;
}

inline double dot(Array2 a, Array2 b){
    return a[0] * b[0] + a[1] * b[1];
}

inline double norm_sq(Array2 a){
    return a[0] * a[0] + a[1] * a[1];
}

inline double cross(Array2 a, Array2 b) {
    return a[0] * b[1] - a[1] * b[0];
}

inline Array2 matmul(const Array2x2& A, const Array2& x){
    return
    {A[0][0] * x[0] + A[0][1] * x[1],
     A[1][0] * x[0] + A[1][1] * x[1]};
}

inline Array2& operator*=(Array2& a, double c) {
    a[0] *= c;
    a[1] *= c;
    return a;
}

inline Array2& operator/=(Array2& a, double c) {
    a[0] /= c;
    a[1] /= c;
    return a;
}

inline constexpr double sq(double x){return x * x;}