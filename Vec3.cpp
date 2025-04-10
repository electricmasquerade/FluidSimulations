#include "Vec3.h"

float dot(const Vec3 &u, const Vec3 &v) {
    return u.e[0]*v.e[0] + u.e[1]*v.e[1] + u.e[2]*v.e[2];
}

Vec3 cross(const Vec3 &u, const Vec3 &v){
    return {u.e[1] * v.e[2] - u.e[2] * v.e[1],
            u.e[2] * v.e[0] - u.e[0] * v.e[2],
            u.e[0] * v.e[1] - u.e[1] * v.e[0]};
}

Vec3 unitVector(const Vec3 &v){
    if (v.length() < 1e-8){
        return {0,0,0};
    }
    return v / v.length();

}

Vec3 unitNormal(const Vec3 &u, const Vec3 &v){
    return unitVector(v - u);
}