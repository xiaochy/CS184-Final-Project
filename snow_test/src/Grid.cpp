#include "Grid.hpp"
#include <iostream>
#include "global.hpp"

GridNode::GridNode(){
    this->active = false;
    this->mass = 0;
    this->old_v.setZero();
    this->new_v.setZero();
    this->v_star.setZero();
    this->force.setZero();
    this->index.setZero();
}

GridNode::~GridNode()
{
}


void GridNode::resetNode()
{
    this->active = false;
    this->mass = 0;
    this->old_v.setZero();
    this->new_v.setZero();
    this->v_star.setZero();
    this->force.setZero();
}

GridMesh::GridMesh(const Bounds3& bbox, const Vector3f& nodeSize,
                   SnowParticleSet* SPS)
{
    this->SPS = SPS;
    this->node_size = nodeSize;

    // estimate the number of nodes based on the size of the bbox
    int num_x, num_y, num_z;
    Vector3f diagonal = bbox.Diagonal();
    Vector3f estimated_num_nodes = diagonal.cwiseProduct(node_size.cwiseInverse());
    num_x = std::max((int)estimated_num_nodes.x() + 1, 3);
    num_y = std::max((int)estimated_num_nodes.y() + 1, 3);
    num_z = std::max((int)estimated_num_nodes.z() + 1, 3);
    this->num_nodes = num_x * num_y * num_z;
    this->x_length = num_x;
    this->y_length = num_y;
    this->z_length = num_z;
    this->eachNodeVolume = node_size.x() * node_size.y() * node_size.z();
    

    // update the bbox
    Vector3f updated_diagonal = Vector3f((float)num_x * node_size.x(), (float)num_y * node_size.y(), (float)num_z * node_size.z());
    std::cout << "x:" << (float)num_x * node_size.x() << std::endl;
    std::cout << "y:" << (float)num_y * node_size.y() << std::endl;
    std::cout << "z:" << (float)num_y * node_size.y() << std::endl;
    Vector3f diff = updated_diagonal - diagonal;
    Vector3f pmin = bbox.pMin - 0.5 * diff;
    Vector3f pmax = bbox.pMax + 0.5 * diff;
    this->bbox = Bounds3(pmin, pmax);

    // contruct the nodes
    this->gridnodes.resize(this->num_nodes);
    for (int i = 0; i < num_x; i++) {
        for (int j = 0; j < num_y; j++) {
            for (int k = 0; k < num_z; k++) {
                GridNode* node = new GridNode();
                node->index = Vector3i(i, j, k);
                this->gridnodes[num_y * num_z * i + num_z * j + k] = node;
            }
        }
    }
}

GridMesh::~GridMesh()
{
    // TODO check if this destroy way is okay
    for (auto& oneNode : gridnodes)
    {
        delete (oneNode);
    }
}

// Maps mass and velocity to the grid
void GridMesh::initialize_grid_mass_velocity()
{
    //SPSbbox = Bounds3(SPS->particles[0]->position);
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetNode();
    }
    //for (auto p : SPS->particles)
    for (int i = 0; i < SPS->particles.size(); i++)
    {
        SnowParticle *p = SPS->particles[i];
        Vector3f particle_pos = p->position;
        Vector3f particle_velocity = p->velocity;
        //SPSbbox = Union(SPSbbox, particle_pos);
        int x_begin = floor((particle_pos.x() - bbox.pMin.x()) / node_size.x());
        int y_begin = floor((particle_pos.y() - bbox.pMin.y()) / node_size.y());
        int z_begin = floor((particle_pos.z() - bbox.pMin.z()) / node_size.z());

        for (int i = x_begin - 1; i <= x_begin + 2; i++)
        {
            for (int j = y_begin - 1; j <= y_begin + 2; j++)
            {
                for (int k = z_begin - 1; k <= z_begin + 2; k++)
                {
                    if (i >= 0 && j >= 0 && k >= 0 && i < x_length && j < y_length && k < z_length)
                    {
                        // // call helper function: from (i,j,k)->node of gridnode
                        GridNode *node = get_GridNode(i, j, k);
                        float offset_x = (particle_pos.x()-bbox.pMin.x()) / node_size.x() - (float)i;
                        float offset_y = (particle_pos.y()-bbox.pMin.y()) / node_size.y() - (float)j;
                        float offset_z = (particle_pos.z()-bbox.pMin.z()) / node_size.z() - (float)k;
                        float wx = bspline(offset_x);
                        float wy = bspline(offset_y);
                        float wz = bspline(offset_z);
                        float weight = wx * wy * wz;
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        p->weights[p_i * 16 + p_j * 4 + p_k] = weight;
                        Vector3f wGrad(
                            wy * wz * bsplineSlope(offset_x) / node_size.x(),
                            wz * wx * bsplineSlope(offset_y) / node_size.y(),
                            wx * wy * bsplineSlope(offset_z) / node_size.z());
                        p->weight_gradient[p_i * 16 + p_j * 4 + p_k] = wGrad;
                        node->mass += weight * p->mass;
                        node->old_v += weight * p->mass * p->velocity;
#ifdef ENABLE_EFFECTIVE                        
                        if (!node->active) {
                            effectiveNodes.push_back(node);
                            node->active = true;
                        }
#else
                        node->active = true;
#endif                        
                   }
                }
            }
        }
    }
}

// after the first time initialize, we will call it to calculate node's mass, velocity and force
void GridMesh::rasterize_particles_to_grid()
{
    //SPSbbox = Bounds3(SPS->particles[0]->position);
#ifndef ENABLE_EFFECTIVE
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetNode();
    }
#else
    for (auto node: effectiveNodes){
        node->resetNode();
    }
    effectiveNodes.clear();
#endif

    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->position;
        Vector3f particle_velocity = p->velocity;
        //SPSbbox = Union(SPSbbox, particle_pos);
        int x_begin = floor((particle_pos.x() - bbox.pMin.x()) / node_size.x());
        int y_begin = floor((particle_pos.y() - bbox.pMin.y()) / node_size.y());
        int z_begin = floor((particle_pos.z() - bbox.pMin.z()) / node_size.z());

        for (int i = x_begin - 1; i <= x_begin + 2; i++)
        {
            for (int j = y_begin - 1; j <= y_begin + 2; j++)
            {
                for (int k = z_begin - 1; k <= z_begin + 2; k++)
                {
                    if (i >= 0 && j >= 0 && k >= 0 && i < x_length && j < y_length && k < z_length)
                    {
                        // call helper function: from (i,j,k)->node of gridnode
                        GridNode *node = get_GridNode(i, j, k);
                        float offset_x = (particle_pos.x()-bbox.pMin.x()) / node_size.x() - (float)i;
                        float offset_y = (particle_pos.y()-bbox.pMin.y()) / node_size.y() - (float)j;
                        float offset_z = (particle_pos.z()-bbox.pMin.z()) / node_size.z() - (float)k;
                        float wx = bspline(offset_x);
                        float wy = bspline(offset_y);
                        float wz = bspline(offset_z);
                        float weight = wx * wy * wz;
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        int idx = p_i * 16 + p_j * 4 + p_k;
                        assert(idx >= 0 && idx < 64);
                        p->weights[idx] = weight;
                        Vector3f wGrad(
                            wy * wz * bsplineSlope(offset_x) / node_size.x(),
                            wz * wx * bsplineSlope(offset_y) / node_size.y(),
                            wx * wy * bsplineSlope(offset_z) / node_size.z());
                        p->weight_gradient[idx] = wGrad;
                        // TODO: need to write class particles: mass
                        node->mass += weight * p->mass;
                        //assert(!p->velocity.hasNaN());
                        node->old_v += weight * p->mass * p->velocity;
                        node->force -= p->energy_derivative() * wGrad;
#ifdef ENABLE_EFFECTIVE                        
                        if (!node->active) {
                            effectiveNodes.push_back(node);
                            node->active = true;
                        }
#else
                        node->active = true;
#endif                                                
                    }
                }
            }
        }
    }
}

void GridMesh::calculate_particle_volume() const
{
    float node_volume = this->eachNodeVolume;
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->position;
        Vector3f particle_velocity = p->velocity;
        p->density = 0;
        int x_begin = floor((particle_pos.x() - bbox.pMin.x()) / node_size.x());
        int y_begin = floor((particle_pos.y() - bbox.pMin.y()) / node_size.y());
        int z_begin = floor((particle_pos.z() - bbox.pMin.z()) / node_size.z());
        
        for (int i = x_begin - 1; i <= x_begin + 2; i++)
        {
            for (int j = y_begin - 1; j <= y_begin + 2; j++)
            {
                for (int k = z_begin - 1; k <= z_begin + 2; k++)
                {
                    if (i >= 0 && j >= 0 && k >= 0 && i < x_length && j < y_length && k < z_length)
                    {
                        // call helper function: from (i,j,k)->node of gridnode
                        GridNode *node = get_GridNode(i, j, k);
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        float weight = p->weights[p_i * 16 + p_j * 4 + p_k];
                        p->density += weight * node->mass;
                    }
                }
            }
        }
        p->density /= node_volume;
        p->volume = p->mass / p->density;
    }
}

void GridMesh::update_node_velocity_star() {
    Vector3f gravity(0, GRAVITY,0);
#ifdef ENABLE_EFFECTIVEx
    for (auto node: effectiveNodes)
#else
    for (auto node: gridnodes)
#endif
    {
        Vector3i index = node->index;
        //std::cout << index.y() << std::endl;
        if (index.y() == 0){
            gravity = Vector3f(0,0,0);
        }
        if (node->mass) {
            node->old_v /= node->mass; // step 1: update node->old_v (we don't divide it before)
            node->v_star = node->old_v + deltaT * (gravity - node->force / node->mass);
        } else {
            node->old_v = Vector3f::Zero();
            node->v_star = Vector3f::Zero();
        }
        assert(!node->old_v.hasNaN());
        assert(!node->v_star.hasNaN());
    }
}

void GridMesh::collision_grid_node()
{
#ifdef ENABLE_EFFECTIVE
    for (auto node: effectiveNodes)
#else
    for (auto &node: gridnodes)
#endif
    {
        Vector3i node_idx = node->index;
        Vector3f node_pos = node->index.cast<float>().cwiseProduct(node_size);
        Vector3f node_tmp_pos = node_pos + deltaT * node->v_star;

        assert(!node->v_star.hasNaN());

        Vector3f Vco, Vrel, Vn, Vt;

        auto f = [&]() {
            float vn = Vn.norm();
            float vt = Vt.norm();
            float mu = SPS->particles[0]->m->sticky;
            assert(!Vt.hasNaN());
            if (vt <= mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = (1. - mu * vn / vt) * Vt;
            }
            assert(!Vrel.hasNaN());
            node->v_star = Vrel + Vco;
        };
        
        // if collision with the x-plane
        if (node_tmp_pos.x() > bbox.pMax.x() - bbox.pMin.x() || node_tmp_pos.x() < 0)
        {
            Vco = Vector3f::Zero();
            Vrel = node->v_star - Vco;
            Vn = Vector3f(Vrel.x(), 0, 0);
            Vt = Vector3f(0, Vrel.y(), Vrel.z());
            f();
        }
        // if collision with y-plane
        if (node_tmp_pos.y() > bbox.pMax.y() - bbox.pMin.y() || node_tmp_pos.y() < 0)
        {
            Vco = Vector3f::Zero();
            Vrel = node->v_star - Vco;
            Vn = Vector3f(0, Vrel.y(), 0);
            Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            f();
        }
        // if collision with z-plane
        if (node_tmp_pos.z() > bbox.pMax.z() - bbox.pMin.z() || node_tmp_pos.z() < 0)
        {
            Vco = Vector3f::Zero();
            Vrel = node->v_star - Vco;
            Vn = Vector3f(0, 0, Vrel.z());
            Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            f();
        }
    }
}

void GridMesh::update_particle_velocity(){
#ifdef ENABLE_EFFECTIVE
    for (auto node: effectiveNodes)
#else
    for (auto node: gridnodes)
#endif
    {
        node->new_v = node->v_star; 
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->position;
        p->v_PIC = Vector3f(0,0,0);
        p->v_FLIP = p->velocity;
        int x_begin = floor((particle_pos.x() - bbox.pMin.x()) / node_size.x());
        int y_begin = floor((particle_pos.y() - bbox.pMin.y()) / node_size.y());
        int z_begin = floor((particle_pos.z() - bbox.pMin.z()) / node_size.z());

        for (int i = x_begin - 1; i <= x_begin + 2; i++)
        {
            for (int j = y_begin - 1; j <= y_begin + 2; j++)
            {
                for (int k = z_begin - 1; k <= z_begin + 2; k++)
                {
                    if (i >= 0 && j >= 0 && k >= 0 && i < x_length && j < y_length && k < z_length)
                    {
                        // call helper function: from (i,j,k)->node of gridnode
                        GridNode *node = get_GridNode(i, j, k);
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        assert(!node->new_v.hasNaN());
                        assert(!(node->new_v - node->old_v).hasNaN());

                        p->v_PIC += node->new_v*p->weights[p_i * 16 + p_j * 4 + p_k];
                        p->v_FLIP += (node->new_v - node->old_v)*p->weights[p_i * 16 + p_j * 4 + p_k];
                    }
                }
            }
        }
        p->update_deform_gradient();
        p->update_velocity();
    }
}

void GridMesh::collision_grid_particle()
{
    // Traverse all the particles and test whether it will collide with the wall
    for (auto p : SPS->particles)
    {
        assert(!p->position.hasNaN());
        assert(!p->velocity.hasNaN());
        Vector3f pos = p->position;
        Vector3f tmp_pos = pos + deltaT * p->velocity;

        Vector3f Vco, Vrel, Vn, Vt;

        auto f = [&]() {
            float vn = Vn.norm();
            float vt = Vt.norm();
            float mu = SPS->particles[0]->m->sticky;
            assert(!Vt.hasNaN());
            if (vt <= mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = (1.f - mu * vn / vt) * Vt;
            }
            try {
                if (Vrel.hasNaN()) {
                    throw std::runtime_error("1234567989");
                }
            } catch (...) {
                std::cout << Vrel << std::endl;
                std::cout << Vco << std::endl;
                std::cout << Vn << std::endl;
                std::cout << Vt << std::endl;
            }
            p->velocity = Vrel + Vco;
        };
        // if collision with the x-plane
        if (tmp_pos.x() > bbox.pMax.x() || tmp_pos.x() < bbox.pMin.x())
        {
            // since the wall is static
            Vco = Vector3f::Zero();
            Vrel = p->velocity - Vco;
            Vn = Vector3f(Vrel.x(), 0, 0);
            Vt = Vector3f(0, Vrel.y(), Vrel.z());
            f();
        }
        // if collision with y-plane
        if (tmp_pos.y() > bbox.pMax.y() || tmp_pos.y() < bbox.pMin.y())
        {
            Vco = Vector3f::Zero();
            Vrel = p->velocity - Vco;
            Vn = Vector3f(0, Vrel.y(), 0);
            Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            f();
        }
        // if collision with z-plane
        if (tmp_pos.z() > bbox.pMax.z() || tmp_pos.z() < bbox.pMin.z())
        {
            Vco = Vector3f::Zero();
            Vrel = p->velocity - Vco;
            Vn = Vector3f(0, 0, Vrel.z());
            Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            f();
        }
    }
}

void GridMesh::update_particle_position(){
    for (auto p : SPS->particles)
    {
        p->update_pos();
    }
}

