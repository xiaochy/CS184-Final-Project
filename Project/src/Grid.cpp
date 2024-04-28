#include "Grid.hpp"
#include "Bounds3.hpp"
#include "global.hpp"

using namespace Eigen;

GridNode::GridNode(float mass, Vector3i index) : 
    mass(mass), index(index), force(Vector3f{}) {
        this->active = false;
        this->mass = mass;
        this->old_v.setZero();
        this->new_v.setZero();
        this->v_star.setZero();
        this->index.setZero();
        this->force.setZero();
}

Grid::Grid(Bounds3& bbox, const Vector3f node_size, SnowParticleSet* sps){
    this->Global_Set = sps;
    this->node_size = node_size;

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
    this->global_bbox = Bounds3(pmin, pmax);

    // contruct the nodes
    this->gridnodes.resize(this->num_nodes);
    for (int i = 0; i < num_x; i++) {
        for (int j = 0; j < num_y; j++) {
            for (int k = 0; k < num_z; k++) {
                GridNode* node = new GridNode();
                node->index = Vector3f(i, j, k);
                this->gridnodes[num_y * num_z * i + num_z * j + k] = node;
            }
        }
    }
}

Grid::~Grid() {
    for (auto &node: gridnodes) {
        delete node;
    }
}

void GridNode::update_velocity_star()
{
    Vector3f gravity(0, GRAVITY,0);
    v_star = old_v + deltaT * (gravity-force / mass);
}

void GridNode::explicit_velocity_update()
{
    new_v = v_star;
}

// Maps volume from the grid to particles
// This should only be called once, at the beginning of the simulation
// after that, each iteration should call Rasterize_Particles_to_Grid()
// begin from iteration 1
void Grid::calculateVolumes() const
{
    float node_volume = node_size.x() * node_size.y() * node_size.z();
    for (auto p : Global_Set->particles)
    {
        Vector3f particle_pos = p->old_pos;
        Vector3f particle_velocity = p->old_v;
        p->density = 0;
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
        p->density /= (node_size.x() * node_size.y() * node_size.z());
        p->volume = p->mass / p->density;
    }
}
// n: old  n+1: new
// At the end of an iteration, we will set old to be the new, because in the next iteration,
// we should use old to calculate new
// we should set paticle's old to be the new
// Maps mass and velocity to the grid
void Grid::Simulate_Once()
{
    for (auto p : Global_Set->particles)
    {
        Vector3f particle_pos = p->old_pos;
        Vector3f particle_velocity = p->old_v;
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
                        GridNode *node = get_GridNode(i, j, k);
                        node->mass = 0;
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
                        node->old_v += weight * p->old_v * p->mass;
                        node->force -= p->volume_cauchy_stress() * wGrad;
                    }
                }
            }
        }
    }
    int num_nodes = x_length * y_length * z_length;
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->old_v /= gridnodes[i]->mass;
    }
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->update_velocity_star();
    }
    // collision between nodes and objects
    collision_grid_node();
    // TODO : collision_object_node();
    // maybe a little bit slow
    for (int i = 0; i < num_nodes; i++)
    {
        gridnodes[i]->explicit_velocity_update();
    }
    for (auto p : Global_Set->particles)
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
                        GridNode *node = get_GridNode(i, j, k);
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
    collision_grid_particle();
    // TODO: collision_object_particle();
    for (auto p : Global_Set->particles)
    {
        p->update_pos();
    }
}

void Grid::collision_grid_node()
{
    for (int i = 0; i < num_nodes; i++)
    {
        GridNode *node = gridnodes[i];
        Vector3i node_idx = node->index;
        Vector3i node_pos = node->index.cwiseProduct(node_size);
        Vector3f node_tmp_pos = node_pos + deltaT * node->v_star;
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
            float mu = Global_Set->particles[0]->m->sticky;
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
            float mu = Global_Set->particles[0]->m->sticky;
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
            float mu = Global_Set->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                VRel.setZero();
            }
            else
            {
                Vrel = Vt + mu * vn * Vt.normalized();
            }
            node->v_star = Vrel + Vco;
        }
    }
}

void Grid::collision_grid_particle()
{
    // Traverse all the particles and test whether it will collide with the wall
    for (auto p : Global_Set->particles)
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
            float mu = Global_Set->particles[0]->m->sticky;
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
        if (node_tmp_pos.y() > bbox.pMax.y() || node_tmp_pos.y() < bbox.pMin.y())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->new_v - Vco;
            Vector3f Vn = Vector3f(0, Vrel.y(), 0);
            Vector3f Vt = Vector3f(Vrel.x(), 0, Vrel.z());
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = Global_Set->particles[0]->m->sticky;
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
        if (node_tmp_pos.z() > bbox.pMax.z() || node_tmp_pos.z() < bbox.pMin.z())
        {
            Vector3f Vco(0, 0, 0);
            Vector3f Vrel = p->new_v - Vco;
            Vector3f Vn = Vector3f(0, 0, Vrel.z());
            Vector3f Vt = Vector3f(Vrel.x(), Vrel.y(), 0);
            float vn = Vn.norm();
            float vt = Vt.norm();
            // TODO this mu calculation should be cleverer
            float mu = Global_Set->particles[0]->m->sticky;
            if (vt < -mu * vn)
            {
                VRel.setZero();
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

void Grid::collision_object_node()
{
}

void Grid::collision_object_particle()
{
}
