#ifndef GRID_H
#define GRID_H

#include "CGL/matrix3x3.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "Particle.h"
#include "global.h"
#include <vector>

using namespace CGL;

class GridNode
{
public:
    float mass;
    Vector3f old_v, new_v;
    Vector3f old_pos, new_po;
    Vector3f force;
    Vector3f v_star;

    GridNode();
    ~GridNode();
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
    std::vector<Particle *> particles;
    std::vector<GridNode *> gridnodes;
    float max_velocity;

    Grid();
    virtual ~Grid();

    // Map particles to grid: update weight. weight_gradient, mass, velocity, force
    void Rasterize_Particles_to_Grid();
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