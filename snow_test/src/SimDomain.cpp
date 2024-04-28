#include "simDomain.hpp"
#include <algorithm>
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

    auto f = [&](){ return SPS->particles[random() % SPS->particles.size()]; };
    std::vector<SnowParticle*> v;
    for (int i = 0; i < 5; i++) {
        v.push_back(f());
        std::cout << v[i]->position << "\n\n";
    }
    std::cout << 4 << std::endl;
    
    deltaT = correctedDeltaT();
    gridMesh->update_node_velocity_star(); // step 4

    for (int i = 0; i < 5; i++) {
        std::cout << v[i]->velocity << "\n\n";
        // std::cout << "vstar:" << v[i]->v_star << "\n\n";
    }
    std::cout << 5 << std::endl;

    // gridMesh->collision_grid_node(); // step 5

    for (int i = 0; i < 5; i++) {
        std::cout << v[i]->velocity << "\n\n";
    }
    std::cout << 6 << std::endl;

    gridMesh->update_particle_velocity(); // step 7 and 8

    for (int i = 0; i < 5; i++) {
        //std::cout << v[i]->v_FLIP << std::endl;

        std::cout << v[i]->velocity << "\n\n";
    }
    std::cout << 78 << std::endl;

    // gridMesh->collision_grid_particle(); // step 9

    for (int i = 0; i < 5; i++) {
        std::cout << v[i]->velocity << "\n\n";
    }
    std::cout << 9 << std::endl;

    gridMesh->update_particle_position(); // step 10
    currentTime += deltaT;

    for (int i = 0; i < 5; i++) {
        std::cout << v[i]->position << ' ';
    }
    std::cout << 10 << std::endl;
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