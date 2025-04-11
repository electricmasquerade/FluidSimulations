
#include "HelperFunctions.h"
#include <cmath>

std::vector<int> HelperFunctions::HSVtoRGB(int hue, int sat, int val) {
    float c = val * sat;
    float x = c * (1 - std::abs(std::fmod(hue / 60.0f, 2) - 1));
    float m = val - c;
    float r, g, b;

    if (hue < 60) {
        r = c;
        g = x;
        b = 0;
    }
    else if (hue < 120) {
        r = x;
        g = c;
        b = 0;
    }
    else if (hue < 180) {
        r = 0;
        g = c;
        b = x;
    }
    else if (hue < 240) {
        r = 0;
        g = x;
        b = c;
    }
    else if (hue < 300) {
        r = x;
        g = 0;
        b = c;
    }
    else {
        r = c;
        g = 0;
        b = x;
    }
    return {
        static_cast<int>(255 * (r+m)),
        static_cast<int>(255 * (g+m)),
        static_cast<int>(255 * (b+m))
    };
}
