#ifndef SIMGLOBAL
#define SIMGLOBAL

#include <eigen3/Eigen/Eigen>

#undef M_PI
#define M_PI 3.141592653589793f

// deltaT can be changed
static float deltaT = 1.e-4;
const float BSPLINE_EPSILON = 1e-4;
const int BSPLINE_RADIUS = 2;

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

// Cubic B-spline shape/basis/interpolation function
// A smooth curve from (0,1) to (1,0)
static float bspline(float x)
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
    if (w < BSPLINE_EPSILON)
        return 0;
    return w;
}
// Slope of interpolation function
static float bsplineSlope(float x)
{
    float abs_x = fabs(x), w;
    if (abs_x < 1)
        return 1.5 * x * abs_x - 2 * x;
    else if (x < 2)
        return -x * abs_x / 2 + 2 * x - 2 * x / abs_x;
    else
        return 0;
    // Clamp between -2/3 and 2/3... if needed
}