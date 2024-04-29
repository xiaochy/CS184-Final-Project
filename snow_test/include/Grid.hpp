#ifndef SNOWSIM_GRID
#define SNOWSIM_GRID

#include "SnowParticle.hpp"
#include "global.hpp"
#include <cstring>
#include <math.h>
#include <stdio.h>
using namespace Eigen;


class SnowParticleMaterial;
class SnowParticle;
class SnowParticleSet;

struct GridNode
{
   public:
    GridNode();
    ~GridNode();
    float mass; 
    bool active; // the gridnode is active if its being used in the current timestep (it has some particles nearby)
    Vector3f old_v, new_v; // the velocity at last timestep, at current timstep
    Vector3f v_star; // the explicit velocity
    Vector3i index; // the index of the gridnode in the gridMesh
    Vector3f force; // the force applied on the node
    void resetNode();
};

class GridMesh
{
   private:
   public:
    Bounds3 bbox; // bounding box of the grid
    size_t x_length, y_length, z_length; // the number of nodes in x, y, z direction
    Vector3f node_size; // the size of each node in x, y, z direction
    SnowParticleSet* SPS; // a set of all snow particles in this grid
    float eachNodeVolume; // volume of a node
    size_t num_nodes; // total number of nodes in the grid
    std::vector<GridNode*> gridnodes; // all nodes in the grid
    int totalEffectiveNodeNum; // number effective nodes at the current time step
    std::vector<GridNode*> effectiveNodes; // effective nodes at the current time step

    void collision_object_node();
    void collision_object_particle();


    // helper function: index->gridnode
    GridNode *get_GridNode(size_t x, size_t y, size_t z) const
    {
        return gridnodes[x * y_length * z_length + y * z_length + z];
    }
    

    GridMesh(const Bounds3& bbox, const Vector3f& nodeSize,
             SnowParticleSet* SPS);
    ~GridMesh();

    // Map particles to grid
    // Initialize the gridnodes' velocity and mass based on each node's nearby particles
    // Only call this once at the beginning
    void initialize_grid_mass_velocity();

    // Recalculate the gridnodes' mass and velocity at each time step
    void rasterize_particles_to_grid();

    // Estimate the particle's denisty and volume based on its surrounding gridnodees
    void calculate_particle_volume() const;

    // Update each node's velocity star
    void update_node_velocity_star();
    
    // Handle the collision between gridnodes and the environment
    void collision_grid_node();

    // Map the nodes' values back to particles, update the the particle velocity after collision
    void update_particle_velocity();

    // Handle the collision between particles and the environment
    void collision_grid_particle();

    // Update the particle position at the end of the timestep
    void update_particle_position();
};

#endif  // SNOWSIM_GRID
