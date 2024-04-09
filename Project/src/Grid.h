#include <vector>

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using namespace CGL;

struct GridNode
{
    float mass;
    // bool active;
    Vector3D old_velocity, new_velocity;
    Vector3D old_position, new_position;
    Vector3D deform_gradient;
    Vector3D force;
    Vector3D velocity;
    Vector3D velocity_star;
};

class Grid
{
public:
    /* store all the particles in the grid */
    Vector3D origin;
    size_t x_length, y_length, z_length; // size of the grid
    Vector3D node_size;
    Vector3D cell_num;
    std::vector<Particle *> particles;
    std::vector<GridNode *> gridnodes;
    float max_velocity;

    // Grid be at least one cell; there must be one layer of cells surrounding all particles
    Grid(Vector3D pos, Vector3D dims, Vector3D cells, PointCloud *obj);
    Grid(const Grid &orig);
    virtual ~Grid();

    // Map particles to grid
    void Rasterize_Particles_to_Grid();
    // void initializeMass();
    void initializeVelocities();
    // Map grid volumes back to particles (first timestep only)
    void calculateVolumes() const;
    // Compute grid velocities
    void explicitVelocities(const Vector3D &gravity);
#if ENABLE_IMPLICIT
    void implicitVelocities();
    void recomputeImplicitForces();
#endif
    // Map grid velocities back to particles
    void updateVelocities() const;

    // Collision detection
    void collisionGrid();
    void collisionParticles() const;

private:
    GridNode *get_GridNode(size_t x, size_t y, size_t z) const
    {
        return gridnodes[x * y_length * z_length + y * z_length + z];
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
    ic:
        void applyPlasticity()
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
};