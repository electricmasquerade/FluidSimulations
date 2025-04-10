
#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

//Define a vector class for handling 3D vectors. For 2D simulation, only the x and y components will be used (z=0).

class Vec3 {
public:
    float e[3];

    Vec3() = default;
    Vec3(float e0, float e1, float e2) : e{e0, e1, e2} {}

    [[nodiscard]] float x() const {return e[0]; }
    [[nodiscard]] float y() const {return e[1]; }
    [[nodiscard]] float z() const {return e[2]; }

    Vec3 operator-() const { return Vec3{-e[0], -e[1], -e[2]}; }
    float operator[](int i) const { return e[i]; }
    float &operator[](int i) { return e[i]; }

    void setX(float x){ e[0] = x; }
    void setY(float y){ e[1] = y; }
    void setZ(float z){ e[2] = z; }


    Vec3& operator+=(const Vec3 &v){
        e[0] += v.e[0];
        e[1] += v.e[1];
        e[2] += v.e[2];
        return *this;
    }

    Vec3& operator*=(float t) {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }

    Vec3& operator/=(float t){
        return *this *= 1/t;
    }

    [[nodiscard]] float length() const {
        return std::sqrt(lengthSquared());
    }

    [[nodiscard]] float lengthSquared() const {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
    }
};

inline std::ostream& operator<<(std::ostream &out, const Vec3 &v){
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline Vec3 operator+(const Vec3 &u, const Vec3 &v){
    return Vec3{u.e[0] + v.e[0], u.e[1] + v.e[1], u.e[2] + v.e[2]};
}

inline Vec3 operator-(const Vec3 &u, const Vec3 &v) {
    return {u.e[0] - v.e[0], u.e[1] - v.e[1], u.e[2] - v.e[2]};
}

inline Vec3 operator*(const Vec3 &u, const Vec3 &v) {
    return {u.e[0] * v.e[0], u.e[1] * v.e[1], u.e[2] * v.e[2]};
}

inline Vec3 operator*(float t, const Vec3 &v) {
    return {t*v.e[0], t*v.e[1], t*v.e[2]};
}

inline Vec3 operator*(const Vec3 &v, float t) {
    return t * v;
}

inline Vec3 operator/(Vec3 v, float t) {
    return (1/t) * v;
}

using Position3 = Vec3;

float dot(const Vec3 &u, const Vec3 &v);

Vec3 cross(const Vec3 &u, const Vec3 &v);

Vec3 unitVector(const Vec3 &v);

Vec3 unitNormal(const Vec3 &u, const Vec3 &v);

Vec3 tangentVector(const Vec3 &u, const Vec3 &v);

#endif //VEC3_H
