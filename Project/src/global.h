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