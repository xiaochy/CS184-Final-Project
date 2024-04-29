#include "Grid.hpp"

GridCell::GridCell()
{
    // active = false;
    // // the mass will be mapped from adj pars each time of simulation
    // mass = 0;
    this->velocity.setZero();
    //velocityStar.setZero();
    // force.setZero();
    // affectedParticles = {};
    // affectedParticlesWeightID = {};

    // start
    this->active = false;
    this->mass = 0;
    //this->old_v.setZero();
    //this->new_v.setZero();
    this->v_star.setZero();
    this->index.setZero();
    this->force.setZero();
    // end

    error.setZero();
    residual.setZero();
    p.setZero();
    Ep.setZero();
    Er.setZero();
    rEr = 0;
}


GridCell::~GridCell()
{
}


void GridCell::resetCell()
{
    this->velocity.setZero();
    //v_star.setZero();
    // start
    this->active = false;
    this->mass = 0;
    this->old_v.setZero();
    this->new_v.setZero();
    this->v_star.setZero();
    //this->index.setZero();
    this->force.setZero();
    // end

    error.setZero();
    residual.setZero();
    p.setZero();
    Ep.setZero();
    Er.setZero();
    rEr = 0;
}

GridMesh::GridMesh(const Bounds3& bbox, const Vector3f& nodeSize,
                   SnowParticleSet* SPS)
{
    // our implementation start
    this->SPS = SPS;
    this->node_size = nodeSize;

    // estimate the number of nodes based on the size of the bbox
    int num_x, num_y, num_z;
    Vector3f diagonal = bbox.Diagonal();
    Vector3f estimated_num_nodes = diagonal.cwiseProduct(node_size.cwiseInverse());
    num_x = std::max((int)estimated_num_nodes.x() + 1, 3);
    num_y = std::max((int)estimated_num_nodes.y() + 1, 3);
    num_z = std::max((int)estimated_num_nodes.z() + 1, 3);
    // std::cout << "x:" << num_x << std::endl;
    // std::cout << "y:" << num_y << std::endl;
    // std::cout << "z:" << num_z << std::endl;
    this->num_nodes = num_x * num_y * num_z;
    this->x_length = num_x;
    this->y_length = num_y;
    this->z_length = num_z;
    this->eachCellVolume = node_size.x() * node_size.y() * node_size.z();
    

    // update the bbox
    Vector3f updated_diagonal = Vector3f((float)num_x * node_size.x(), (float)num_y * node_size.y(), (float)num_z * node_size.z());
    //Vector3f updated_diagonal = Vector3f(num_x, num_y, num_z);
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
                GridCell* node = new GridCell();
                node->index = Vector3i(i, j, k);
                this->gridnodes[num_y * num_z * i + num_z * j + k] = node;
            }
        }
    }
}

GridMesh::~GridMesh()
{
    // TODO check if this destroy way is okay
    for (auto& oneCell : gridnodes)
    {
        delete (oneCell);
    }
}

// Maps mass and velocity to the grid
void GridMesh::initializeGridMeshActiveMassAndMomentum()
{
    SPSbbox = Bounds3(SPS->particles[0]->position);
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetCell();
    }
    //for (auto p : SPS->particles)
    for (int i = 0; i < SPS->particles.size(); i++)
    {
        //std::cout << "entered initialize for loop" << std::endl;
        SnowParticle *p = SPS->particles[i];
        Vector3f particle_pos = p->position;
        //Vector3f particle_velocity = p->old_v;
        Vector3f particle_velocity = p->velocity;
        // Vector3f particle_pos = p->position;
        // Vector3f particle_velocity = p->velocity;
        SPSbbox = Union(SPSbbox, particle_pos);
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
                        GridCell *node = get_GridNode(i, j, k);
                        // float offset_x = particle_pos.x() - (float)i * node_size.x();
                        // float offset_y = particle_pos.y() - (float)j * node_size.y();
                        // float offset_z = particle_pos.z() - (float)k * node_size.z();
                        float offset_x = (particle_pos.x()-bbox.pMin.x()) / node_size.x() - (float)i;
                        float offset_y = (particle_pos.y()-bbox.pMin.y()) / node_size.y() - (float)j;
                        float offset_z = (particle_pos.z()-bbox.pMin.z()) / node_size.z() - (float)k;
                        float wx = bspline(offset_x);
                        float wy = bspline(offset_y);
                        float wz = bspline(offset_z);
                        // std::cout << "offset.x:" << offset_x << std::endl;
                        // std::cout << "y:" << offset_y << std::endl;
                        // std::cout << "z:" << offset_z << std::endl;
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
                        // TODO: need to write class particles: mass
                        node->mass += weight * p->mass;
                        // diff: node->old_v & p->old_v
                        //node->old_v += weight * p->old_v * p->mass;
                        node->old_v += weight * p->mass * p->velocity;
                        // node->velocity += weight*p->velocity*p->mass;
                        node->active = true;
                        // node->force -= p->volume_cauchy_stress() * wGrad;
                   }
                }
            }
        }
    }
    //std::cout << "passed initialization" << std::endl;
    // assert(SPS->particles[0]);
    //std::cout << "enter volume" << std::endl;
}

// after the first time initialize, we will call it to calculate node's mass, velocity and force
void GridMesh::rasterize_particles_to_grid()
{
    SPSbbox = Bounds3(SPS->particles[0]->position);
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetCell();
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->position;
        //Vector3f particle_velocity = p->old_v;
        // Vector3f particle_pos = p->position;
        Vector3f particle_velocity = p->velocity;
        SPSbbox = Union(SPSbbox, particle_pos);
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
                        GridCell *node = get_GridNode(i, j, k);
                        // float offset_x = (particle_pos.x() - bbox.pMin()) - i * node_size.x();
                        // float offset_y = particle_pos.y() - j * node_size.y();
                        // float offset_z = particle_pos.z() - k * node_size.z();
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
                        // diff: node->old_v & p->old_v
                        //node->old_v += weight * p->old_v * p->mass;
                        assert(!p->velocity.hasNaN());
                        node->old_v += weight * p->mass * p->velocity;
                        // node->velocity += weight*p->velocity*p->mass;
                        node->active = true;
                        //node->force -= p->volume_cauchy_stress() * wGrad;
                        node->force -= p->energyDerivative() * wGrad;
                    }
                }
            }
        }
    }
}

void GridMesh::calculateParticleVolume() const
{
    //std::cout << "enter volume" << std::endl;
    float node_volume = node_size.x() * node_size.y() * node_size.z();
    //assert(SPS->particles[0]);
    //assert(SPS->particles[1]);
    for (auto p : SPS->particles)
    {
        //std::cout << "enter for loop" << std::endl;
        //std::out << length(sp)
        Vector3f particle_pos = p->position;
        //Vector3f particle_pos = p->old_pos;
        //std::cout << particle_pos.x() << std::endl;
        //Vector3f particle_velocity = p->old_v;
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
                        // std::cout << "inside" << std:: endl;
                        GridCell *node = get_GridNode(i, j, k);
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        float weight = p->weights[p_i * 16 + p_j * 4 + p_k];
                        // std::cout << weight << ' ' << node->mass << std::endl;
                        p->density += weight * node->mass;
                        //std :: cout << weight << std:: endl;
                    }
                }
            }
        }
        p->density /= (node_size.x() * node_size.y() * node_size.z());
        p->volume = p->mass / p->density;
        //std :: cout << p->density << std:: endl;
    }
}

void GridMesh::update_node_velocity_star() {
    Vector3f gravity(0, GRAVITY,0);
    for (int i = 0; i < num_nodes; i++)
    {
        if (gridnodes[i]->mass) {
            gridnodes[i]->old_v /= gridnodes[i]->mass; // step 1: update node->old_v (we don't divide it before)
            gridnodes[i]->v_star = gridnodes[i]->old_v + deltaT * (gravity - gridnodes[i]->force / gridnodes[i]->mass);
        } else {
            gridnodes[i]->old_v = Vector3f::Zero();
            gridnodes[i]->v_star = Vector3f::Zero();
        }
        assert(!gridnodes[i]->old_v.hasNaN());
        assert(!gridnodes[i]->v_star.hasNaN());
        // //gridnodes[i]->new_v = gridnodes[i]->v_star; // step 6 explicit (without implicit)
    }
}

void GridMesh::collision_grid_node()
{
    for (auto &node: gridnodes)
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
            // Vector3f n(0, 1, 0); // normal of the y-plane
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
        
        // node->new_v = node->v_star; // step 6 : explicit update
    }
}

void GridMesh::update_particle_velocity(){
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->new_v = gridnodes[i]->v_star; 
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->position;
        p->v_PIC = Vector3f(0,0,0);
        //p->v_FLIP = Vector3f(0,0,0);
        p->v_FLIP = p->velocity;
        // int x_begin = floor(particle_pos.x() / node_size.x());
        // int y_begin = floor(particle_pos.y() / node_size.y());
        // int z_begin = floor(particle_pos.z() / node_size.z());
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
                        GridCell *node = get_GridNode(i, j, k);
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
        //p->v_FLIP += p->velocity;
        //std::cout << p->v_FLIP << std::endl;
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
            assert(!Vrel.hasNaN());
            p->velocity = Vrel + Vco;
        };
        // if collision with the x-plane
        if (tmp_pos.x() > bbox.pMax.x() || tmp_pos.x() < bbox.pMin.x())
        {
            //std::cout << "collide x-plane" << std::endl;
            // since the wall is static
            Vector3f Vco(0, 0, 0);
            //Vector3f Vrel = p->new_v - Vco;
            Vector3f Vrel = p->velocity - Vco;
            Vector3f Vn = Vector3f(Vrel.x(), 0, 0);
            Vector3f Vt = Vector3f(0, Vrel.y(), Vrel.z());
            f();
            // float vn = Vn.norm();
            // float vt = Vt.norm();
            // float mu = SPS->particles[0]->m->sticky;
            // if (vt < -mu * vn)
            // {
            //     Vrel.setZero();
            // }
            // else
            // {
            //     Vrel = (1. + mu * vn / vt) * Vt;
            // }
            // p->velocity = Vrel + Vco;
            //p->old_v = p->new_v;
        }
        // if collision with y-plane
        if (tmp_pos.y() > bbox.pMax.y() || tmp_pos.y() < bbox.pMin.y())
        {
            //std::cout << "collide y-plane" << std::endl;
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->velocity - Vco;
            Vector3f Vn = Vector3f(0, Vrel.y(), 0);
            Vector3f Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            f();
            // float vn = Vn.norm();
            // float vt = Vt.norm();
            // // TODO this mu calculation should be cleverer
            // float mu = SPS->particles[0]->m->sticky;
            // if (vt < -mu * vn)
            // {
            //     Vrel.setZero();
            // }
            // else
            // {
            //     Vrel = (1. + mu * vn / vt) * Vt;
            // }
            // p->velocity = Vrel + Vco;
            // //p->old_v = p->new_v;
        }
        // if collision with z-plane
        if (tmp_pos.z() > bbox.pMax.z() || tmp_pos.z() < bbox.pMin.z())
        {
            //std::cout << "collide z-plane" << std::endl;
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->velocity - Vco;
            Vector3f Vn = Vector3f(0, 0, Vrel.z());
            Vector3f Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            f();
            // float vn = Vn.norm();
            // float vt = Vt.norm();
            // float mu = SPS->particles[0]->m->sticky;
            // if (vt < -mu * vn)
            // {
            //     Vrel.setZero();
            // }
            // else
            // {
            //     Vrel = (1. + mu * vn / vt) * Vt;
            // }
            // p->velocity = Vrel + Vco;
            // //p->old_v = p->new_v;
        }
    }
}

void GridMesh::update_particle_position(){
    for (auto p : SPS->particles)
    {
        p->update_pos();
    }
}


// #if ENABLE_IMPLICIT
// Solve linear system for implicit velocities
void GridMesh::implicitUpdateGridVelocity()
{
    // TODO
}

void GridMesh::recomputeImplicitForces()
{
    // TODO
}
// #endif
