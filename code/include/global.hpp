#ifndef SIMGLOBAL
#define SIMGLOBAL

#include <cmath>

#undef M_PI
#define M_PI 3.141592653589793f

// deltaT can be changed
static float deltaT = 4e-3;

const float BSPLINE_EPSILON = 1e-4;

static const float EPSILON = 0.001;
static const float kInfinity = std::numeric_limits<float>::max();
static const float FRAMERATE = 60.;
// Adaptive timestep adjustment
static const float CFL = .04;
static const float MAX_TIMESTEP = 5.e-5;
// Percentage that FLIP TAKES, PIC is 1 - that
static const float FLIP_PERCENT = .95;
// Percentage that should be implicit vs explicit
static const float IMPLICIT_RATIO = 0;
static const int MAX_IMPLICIT_ITERS = 30;
static const float MAX_IMPLICIT_ERR = 1.e4;
static const float MIN_IMPLICIT_ERR = 1.e-4;
static const float GRAVITY = -9.8;

inline float clamp(const float &lo, const float &hi, const float &v)
{
    return std::max(lo, std::min(hi, v));
}

inline bool solveQuadratic(const float &a, const float &b, const float &c,
                           float &x0, float &x1)
{
    float discr = b * b - 4 * a * c;
    if (discr < 0)
        return false;
    else if (discr == 0)
        x0 = x1 = -0.5 * b / a;
    else
    {
        float q = (b > 0) ? -0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
        x0 = q / a;
        x1 = c / q;
    }
    if (x0 > x1) std::swap(x0, x1);
    return true;
}

// Cubic B-spline shape/basis/interpolation function
// A smooth curve from (0,1) to (1,0)
inline float bspline(float x)
{
    x = fabs(x);
    float w;
    if (x < 1)
        w = x * x * (x / 2 - 1) + 2 / 3.0;
    else if (x < 2)
        w = x * (x * (-x / 6 + 1) - 2) + 4 / 3.0;
    else
        return 0;
    // Clamp between 0 and 1... if needed
    if (w < BSPLINE_EPSILON) return 0;
    return w;
}
// Slope of interpolation function
inline float bsplineSlope(float x)
{
    float abs_x = fabs(x);
    if (abs_x < 1)
        return 1.5 * x * abs_x - 2 * x;
    else if (x < 2)
        return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
    else
        return 0;
    // Clamp between -2/3 and 2/3... if needed
}

// some compile options
#define ENABLE_TESTONTHERUN true
#define ENABLE_IMPLICIT false

// #define ENABLE_EFFECTIVE

#endif