#include "simDomain.hpp"
#include <algorithm>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

SimDomain::SimDomain(SnowParticleSet* SPS, GridMesh* gridMesh)
{
    currentTime = 0;
    this->SPS = SPS;
    this->gridMesh = gridMesh;
}

SimDomain::~SimDomain()
{
}

void SimDomain::initializeSimulator()
{
    gridMesh->initialize_grid_mass_velocity(); // initialized step 1
    gridMesh->calculate_particle_volume(); // step 2
}

void SimDomain::oneTimeSimulate()
{
    gridMesh->rasterize_particles_to_grid(); // step 1 and step 3
    
    deltaT = correctedDeltaT();
    
    gridMesh->update_node_velocity_star(); // step 4

    gridMesh->collision_grid_node(); // step 5

    gridMesh->update_particle_velocity(); // step 7 and 8

    gridMesh->collision_grid_particle(); // step 9

    gridMesh->update_particle_position(); // step 10
    currentTime += deltaT;
}

float SimDomain::correctedDeltaT()
{
    float prevMaxVelocity = SPS->max_velocity;
    float dt;
    if (prevMaxVelocity > 1.e-8)
    {
        float minNodeSize =
            std::min(std::min(gridMesh->node_size[0], gridMesh->node_size[1]),
                     gridMesh->node_size[2]);
        dt = CFL * minNodeSize / prevMaxVelocity;
        dt = std::min(dt, 1.f / FRAMERATE);
    }
    else
    {
        dt = 1. / FRAMERATE;
    }
    return dt > MAX_TIMESTEP ? MAX_TIMESTEP : dt;
}