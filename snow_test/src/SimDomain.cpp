#include "simDomain.hpp"
#include <iostream>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
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
    gridMesh->initializeGridMeshActiveMassAndMomentum(); // initialized step 1
    gridMesh->calculateParticleVolume(); // step 2
    // gridMesh->calculateVolumes();
}

void SimDomain::oneTimeSimulate()
{
    gridMesh->rasterize_particles_to_grid(); // step 1 and step 3
    // correct delta T after initialize grid mass
    // this is because the first one has to be serial
    // so why not just get the maxV here
    
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
    float prevMaxVelocity = SPS->maxVelocity;
    float dt;
    if (prevMaxVelocity > 1.e-8)
    {
        float minCellSize =
            std::min(std::min(gridMesh->node_size[0], gridMesh->node_size[1]),
                     gridMesh->node_size[2]);
        float dt = CFL * minCellSize / prevMaxVelocity;
        dt = dt > 1. / FRAMERATE ? 1. / FRAMERATE : dt;
    }
    else
    {
        dt = 1. / FRAMERATE;
    }
    return dt > MAX_TIMESTEP ? MAX_TIMESTEP : dt;
}