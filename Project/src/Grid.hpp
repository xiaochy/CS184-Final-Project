#ifndef GRID_H
#define GRID_H

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include <eigen3/Eigen/src/Core/Matrix.h>
#include "Particle.hpp"
#include "global.hpp"
#include "bounds3.hpp"
// #include "bbox.h"
#include <vector>

// using Eigen::Matrix3f;
// using Eigen::Vector3f, Eigen::Vector3i;
using namespace Eigen;

class GridNode
{
public:
    float mass;
    Vector3f old_v, new_v;
    /* never change, is the index on the grid; need to initialize in main */
    Vector3i index;
    Vector3f force;
    Vector3f v_star;
    bool active;

    GridNode(float mass, Vector3i index);
    ~GridNode() = default;
    void update_velocity_star();
    void explicit_velocity_update();
};

class Grid
{
public:
    /* store all the particles in the grid */
    size_t x_length, y_length, z_length; // size of the grid
    size_t num_nodes;
    Vector3f node_size;
    //std::vector<Particle *> particles;
    SnowParticleSet* Global_Set;
    std::vector<GridNode *> gridnodes;
    float max_velocity;
    /* need to add a class Bound. global_bbox is the wall */
    Bounds3 global_bbox;

    Grid();
    virtual ~Grid();

    // Map particles to grid: update weight. weight_gradient, mass, velocity, force
    //void Rasterize_Particles_to_Grid();
    void Simulate_Once();
    // Map grid volumes back to particles (first timestep only)
    void calculateVolumes() const;
    void collision_grid_node();
    void collision_grid_particle();
    void collision_object_node();
    void collision_object_particle();

private:
    GridNode *get_GridNode(size_t x, size_t y, size_t z) const
    {
        return gridnodes[x * y_length * z_length + y * z_length + z];
    }
};

#endif