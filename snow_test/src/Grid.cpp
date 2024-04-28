#include "Grid.hpp"

GridCell::GridCell()
{
    // active = false;
    // // the mass will be mapped from adj pars each time of simulation
    // mass = 0;
    velocity.setZero();
    //velocityStar.setZero();
    // force.setZero();
    // affectedParticles = {};
    // affectedParticlesWeightID = {};

    // start
    this->active = false;
    this->mass = 0;
    this->old_v.setZero();
    this->new_v.setZero();
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
    velocity.setZero();
    //v_star.setZero();
    // start
    this->active = false;
    this->mass = 0;
    this->old_v.setZero();
    this->new_v.setZero();
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
    this->num_nodes = num_x * num_y * num_z;
    this->x_length = num_x;
    this->y_length = num_y;
    this->z_length = num_z;
    this->eachCellVolume = node_size.x() * node_size.y() * node_size.z();

    // update the bbox
    Vector3f updated_diagonal = Vector3f(this->x_length, this->y_length, this->z_length);
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
    SPSbbox = Bounds3(SPS->particles[0]->old_pos);
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetCell();
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->old_pos;
        Vector3f particle_velocity = p->old_v;
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
                        // call helper function: from (i,j,k)->node of gridnode
                        GridCell *node = get_GridNode(i, j, k);
                        node->old_v = Vector3f();
                        float offset_x = particle_pos.x() - i * node_size.x();
                        float offset_y = particle_pos.y() - j * node_size.y();
                        float offset_z = particle_pos.z() - k * node_size.z();
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
                        p->weight_gradient[i * 16 + j * 4 + k] = wGrad;
                        // TODO: need to write class particles: mass
                        node->mass += weight * p->mass;
                        // diff: node->old_v & p->old_v
                        node->old_v += weight * p->old_v * p->mass;
                        // node->velocity += weight*p->velocity*p->mass;
                        node->active = true;
                        // node->force -= p->volume_cauchy_stress() * wGrad;
                    }
                }
            }
        }
    }
}

// after the first time initialize, we will call it to calculate node's mass, velocity and force
void GridMesh::rasterize_particles_to_grid()
{
    SPSbbox = Bounds3(SPS->particles[0]->old_pos);
    for (int i = 0; i < num_nodes; i++){
        gridnodes[i]->resetCell();
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->old_pos;
        Vector3f particle_velocity = p->old_v;
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
                        // call helper function: from (i,j,k)->node of gridnode
                        GridCell *node = get_GridNode(i, j, k);
                        node->old_v = Vector3f();
                        float offset_x = particle_pos.x() - i * node_size.x();
                        float offset_y = particle_pos.y() - j * node_size.y();
                        float offset_z = particle_pos.z() - k * node_size.z();
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
                        p->weight_gradient[i * 16 + j * 4 + k] = wGrad;
                        // TODO: need to write class particles: mass
                        node->mass += weight * p->mass;
                        // diff: node->old_v & p->old_v
                        node->old_v += weight * p->old_v * p->mass;
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
    float node_volume = node_size.x() * node_size.y() * node_size.z();
    for (auto p : SPS->particles)
    {
        // Vector3f particle_pos = p->position;
        Vector3f particle_pos = p->old_pos;
        Vector3f particle_velocity = p->old_v;
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
                        GridCell *node = get_GridNode(i, j, k);
                        int p_i = i - (x_begin - 1);
                        int p_j = j - (y_begin - 1);
                        int p_k = k - (z_begin - 1);
                        float weight = p->weights[p_i * 16 + p_j * 4 + p_k];
                        p->density += weight * node->mass;
                    }
                }
            }
        }
        p->density /= (node_size.x() * node_size.y() * node_size.z());
        p->volume = p->mass / p->density;
    }
}

void GridMesh::update_node_velocity_star() {
    Vector3f gravity(0, GRAVITY,0);
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->old_v /= gridnodes[i]->mass; // step 1: update node->old_v (we don't divide it before)
        gridnodes[i]->v_star = gridnodes[i]->old_v + deltaT * (gravity - gridnodes[i]->force / gridnodes[i]->mass);
        gridnodes[i]->new_v = gridnodes[i]->v_star; // step 6 explicit (without implicit)
    }
}

void GridMesh::collision_grid_node()
{
    for (int i = 0; i < num_nodes; i++)
    {
        GridCell *node = gridnodes[i];
        Vector3i node_idx = node->index;
        Vector3i node_pos = node->index.cast<float>().cwiseProduct(node_size).cast<int>();
        Vector3f node_tmp_pos = node_pos.cast<float>() + deltaT * node->v_star;
        // if collision with the x-plane
        if (node_tmp_pos.x() > bbox.pMax.x() || node_tmp_pos.x() < bbox.pMin.x())
        {
            // since the wall is static
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = node->v_star - Vco;
            Vector3f Vn = Vector3f(Vrel.x(), 0, 0);
            Vector3f Vt = Vector3f(0, Vrel.y(), Vrel.z());
            float vn = Vn.norm();
            float vt = Vt.norm();
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            node->v_star = Vrel + Vco;
        }
        // if collision with y-plane
        if (node_tmp_pos.y() > bbox.pMax.y() || node_tmp_pos.y() < bbox.pMin.y())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = node->v_star - Vco;
            Vector3f Vn = Vector3f(0, Vrel.y(), 0);
            Vector3f Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            node->v_star = Vrel + Vco;
        }
        // if collision with z-plane
        if (node_tmp_pos.z() > bbox.pMax.z() || node_tmp_pos.z() < bbox.pMin.z())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = node->v_star - Vco;
            Vector3f Vn = Vector3f(0, 0, Vrel.z());
            Vector3f Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            node->v_star = Vrel + Vco;
        }
    }
}

void GridMesh::update_particle_velocity(){
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->explicit_velocity_update();
    }
    for (auto p : SPS->particles)
    {
        Vector3f particle_pos = p->old_pos;
        p->v_PIC = Vector3f(0,0,0);
        p->v_FLIP = Vector3f(0,0,0);
        int x_begin = floor(particle_pos.x() / node_size.x());
        int y_begin = floor(particle_pos.y() / node_size.y());
        int z_begin = floor(particle_pos.z() / node_size.z());

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
                        p->v_PIC += node->new_v*p->weights[p_i * 16 + p_j * 4 + p_k];
                        p->v_FLIP += (node->new_v - node->old_v)*p->weights[p_i * 16 + p_j * 4 + p_k];
                    }
                }
            }
        }
        p->v_FLIP += p->old_v;
        p->update_deform_gradient();
        p->update_velocity();
    }
}

void GridMesh::collision_grid_particle()
{
    // Traverse all the particles and test whether it will collide with the wall
    for (auto p : SPS->particles)
    {
        Vector3f pos = p->old_pos;
        Vector3f tmp_pos = pos + deltaT * p->new_v;
        // if collision with the x-plane
        if (tmp_pos.x() > bbox.pMax.x() || tmp_pos.x() < bbox.pMin.x())
        {
            // since the wall is static
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->new_v - Vco;
            Vector3f Vn = Vector3f(Vrel.x(), 0, 0);
            Vector3f Vt = Vector3f(0, Vrel.y(), Vrel.z());
            float vn = Vn.norm();
            float vt = Vt.norm();
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            p->new_v = Vrel + Vco;
            p->old_v = p->new_v;
        }
        // if collision with y-plane
        if (tmp_pos.y() > bbox.pMax.y() || tmp_pos.y() < bbox.pMin.y())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->new_v - Vco;
            Vector3f Vn = Vector3f(0, Vrel.y(), 0);
            Vector3f Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            p->new_v = Vrel + Vco;
            p->old_v = p->new_v;
        }
        // if collision with z-plane
        if (tmp_pos.z() > bbox.pMax.z() || tmp_pos.z() < bbox.pMin.z())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->new_v - Vco;
            Vector3f Vn = Vector3f(0, 0, Vrel.z());
            Vector3f Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = SPS->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                Vrel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            p->new_v = Vrel + Vco;
            p->old_v = p->new_v;
        }
    }
}

void GridMesh::update_particle_position(){
    for (auto p : SPS->particles)
    {
        p->update_pos();
    }
}

// void GridMesh::updateVelocityInGrids(const Vector3f& gravity)
// {
// #pragma omp parallel for
//     for (int iCell = 0; iCell < totalEffectiveCellNum; iCell++)
//     {
//         // get vel
//         {
//             GridCell* cell = effectiveCells[iCell];
//             cell->velocity /= cell->mass;
//         }
//         // get vel star
//         {
//             // First, compute the forces
//             GridCell* cell = effectiveCells[iCell];
//             Vector3f force = Vector3f::Zero();
//             for (int i = 0; i < cell->affectedParticles.size(); i++)
//             {
//                 SnowParticle* p = cell->affectedParticles[i];
//                 int countForAdjGrid = cell->affectedParticlesWeightID[i];
//                 float w = p->weights[countForAdjGrid];
//                 Vector3f wGrad = p->weight_gradient[countForAdjGrid];
//                 Matrix3f energyD = p->energyDerivative();
//                 if (w > BSPLINE_EPSILON)
//                 {
//                     force += energyD * wGrad;
//                 }
//             }
//             cell->v_star =
//                 cell->velocity + deltaT * (gravity - force / cell->mass);
//         }
//         // todo consider implicit

//         // collision test
//         {
//             GridCell* cell = gridnodes[iCell];
//             int i = cell->index.x();
//             int j = cell->index.y();
//             int k = cell->index.z();
//             // collision test
//             Vector3f oldCellPos((i + 0.5) * node_size.x(),
//                                 (j + 0.5) * node_size.y(),
//                                 (k + 0.5) * node_size.z());
//             Vector3f tempPos = oldCellPos + deltaT * cell->v_star;
//             // TODO consider wrap this fixing vel as a function
//             // x dir
//             if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = cell->v_star - VSolid;
//                 Vector3f Vn = Vector3f(VRel.x(), 0, 0);
//                 Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 // TODO this mu calculation should be cleverer
//                 float mu = SPS->particles[0]->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 cell->v_star = VRel + VSolid;
//                 // std::cout << " cell velocity star is " << cell->velocityStar
//                 //           << std::endl;
//             }
//             // y dir
//             if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = cell->v_star - VSolid;
//                 Vector3f Vn = Vector3f(0, VRel.y(), 0);
//                 Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 // TODO this mu calculation should be cleverer
//                 float mu = SPS->particles[0]->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 cell->v_star = VRel + VSolid;
//                 // std::cout << " cell velocity star is " << cell->velocityStar
//                 //           << std::endl;
//             }
//             // z dir
//             if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = cell->v_star - VSolid;
//                 Vector3f Vn = Vector3f(0, 0, VRel.z());
//                 Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 // TODO this mu calculation should be cleverer
//                 float mu = SPS->particles[0]->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 cell->v_star = VRel + VSolid;
//                 // std::cout << " cell velocity star is " << cell->velocityStar
//                 //           << std::endl;
//             }
//         }
//     }
// }


//     // TODO the above loop shall be combined with collision Grid when Parfor
//     gridCollisionTest();
// }

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

// Map grid velocities back to particles
// void GridMesh::mapVelocityToSPS() const
// {
// #pragma omp parallel for
//     for (int i = 0; i < SPS->particles.size(); i++)
//     {
//         // map back to particles
//         {
//             SnowParticle* p = SPS->particles[i];
//             // position
//             Vector3f pos = p->position - bbox.pMin;
//             // We calculate Vpic and VFLIP velocities separately
//             Vector3f Vpic = Vector3f::Zero();
//             Vector3f Vflip = p->velocity;
//             // Also keep track of velocity gradient
//             Matrix3f& grad = p->v_grad;
//             grad.setZero();
//             // VISUALIZATION PURPOSES ONLY:
//             // Recompute density
//             p->density = 0;

//             for (int countForAdjGrid = 0, px = (int)(pos.x() / x_length),
//                      ipx = px - 1, pxEnd = px + 3;
//                  ipx < pxEnd; ipx++)
//             {
//                 if (ipx >= 0 && ipx < x_length)
//                 {
//                     for (int py = (int)(pos.y() / y_length), ipy = py - 1,
//                              pyEnd = py + 3;
//                          ipy < pyEnd; ipy++)
//                     {
//                         if (ipy >= 0 && ipy < y_length)
//                         {
//                             for (int pz = (int)(pos.z() / z_length),
//                                      ipz = pz - 1, pzEnd = pz + 3;
//                                  ipz < pzEnd; ipz++)
//                             {
//                                 if (ipz >= 0 && ipz < z_length)
//                                 {
//                                     float w = p->weights[countForAdjGrid];
//                                     Vector3f wGrad =
//                                         p->weight_gradient[countForAdjGrid];
//                                     if (w > BSPLINE_EPSILON)
//                                     {
//                                         GridCell* cell = gridnodes[(
//                                             int)(ipx * y_length *
//                                                      z_length +
//                                                  ipy * z_length + ipz)];
//                                         // Particle in cell
//                                         Vpic += cell->v_star * w;
//                                         // Fluid implicit particle
//                                         Vflip += (cell->v_star -
//                                                   cell->velocity) *
//                                                  w;
//                                         // Velocity gradient
//                                         grad += cell->v_star *
//                                                 wGrad.transpose();
//                                         // VISUALIZATION ONLY: Update density
//                                         p->density += w * cell->mass;
//                                     }
//                                 }
//                                 // do it or not, it has goto next cell
//                                 countForAdjGrid++;
//                             }
//                         }
//                         else
//                         {
//                             // do nothing and skip this z edge
//                             countForAdjGrid += 4;
//                         }
//                     }
//                 }
//                 else
//                 {
//                     // do nothing and skip this yz face
//                     countForAdjGrid += 16;
//                 }
//             }
//             // Final velocity is a linear combination of Vpic and VFLIP
//             // components
//             p->velocity = Vflip * FLIP_PERCENT + Vpic * (1. - FLIP_PERCENT);
//             // VISUALIZATION: Update density
//             p->density /= eachCellVolume;
//         }
//         // collision test
//         {
//             // collision test
//             SnowParticle* p = SPS->particles[i];
//             Vector3f tempPos = p->position + deltaT * p->velocity;
//             // x dir
//             if (tempPos.x() > bbox.pMax.x() || tempPos.x() < bbox.pMin.x())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = p->velocity - VSolid;
//                 Vector3f Vn = Vector3f(VRel.x(), 0, 0);
//                 Vector3f Vt = Vector3f(0, VRel.y(), VRel.z());
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 float mu = p->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 p->velocity = VRel + VSolid;
//             }
//             // y dir
//             if (tempPos.y() > bbox.pMax.y() || tempPos.y() < bbox.pMin.y())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = p->velocity - VSolid;
//                 Vector3f Vn = Vector3f(0, VRel.y(), 0);
//                 Vector3f Vt = Vector3f(VRel.x(), 0, VRel.z());
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 float mu = p->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 p->velocity = VRel + VSolid;
//             }
//             // z dir
//             if (tempPos.z() > bbox.pMax.z() || tempPos.z() < bbox.pMin.z())
//             {
//                 Vector3f VSolid(0, 0, 0);
//                 Vector3f VRel = p->velocity - VSolid;
//                 Vector3f Vn = Vector3f(0, 0, VRel.z());
//                 Vector3f Vt = Vector3f(VRel.x(), VRel.y(), 0);
//                 float vn = Vn.norm();
//                 float vt = Vt.norm();
//                 float mu = p->m->sticky;
//                 if (vt < mu * vn || vt < 1.e-12)
//                 {
//                     VRel.setZero();
//                 }
//                 else
//                 {
//                     VRel = (1. - mu * vn / vt) * Vt;
//                 }
//                 p->velocity = VRel + VSolid;
//             }
//             // TODO there are more collisions
//             // Finish collision test
//             // Now push Particles Postion
//             p->update_pos();
//             // update pure elastic deformation
//             p->updatePureElasticGradient();
//             // push excessive part to plastic deformation
//             p->updateCombinedPElasticGradient();
//         }
//     }
// }