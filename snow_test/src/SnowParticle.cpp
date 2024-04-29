#include "SnowParticle.hpp"
#include <iostream>

SnowParticleMaterial::SnowParticleMaterial()
{
    mu = youngsModule / (2. + 2. * PoissonsRatio);
    lambda = youngsModule * PoissonsRatio /
             ((1. + PoissonsRatio) * (1. - 2. * PoissonsRatio));
}

SnowParticleMaterial::~SnowParticleMaterial()
{
}

SnowParticle::SnowParticle() : m(nullptr)
{
    // volume mass and density will be calculated initially
    volume = 0;
    mass = 0;
    density = 0;
    position = Vector3f(0, 0, 0);
    velocity = Vector3f(0, 0, 0);
    v_grad = Matrix3f::Zero();
    deform_elastic_grad = Matrix3f::Identity();
    deform_plastic_grad = Matrix3f::Identity();
    svd_u = Matrix3f::Identity();
    svd_v = Matrix3f::Identity();
    svd_s = Vector3f(1, 1, 1);
    v_PIC = Vector3f(0, 0, 0);
    v_FLIP = Vector3f(0, 0, 0);
}

SnowParticle::SnowParticle(const Vector3f& pos, const Vector3f& vel,
                           const float Mass, SnowParticleMaterial* material)
    : mass(Mass), position(pos), velocity(vel), m(material)
{
    // volume mass and density will be calculated initially
    volume = 0;
    density = 0;
    v_grad = Matrix3f::Zero();
    deform_elastic_grad = Matrix3f::Identity();
    deform_plastic_grad = Matrix3f::Identity();
    svd_u = Matrix3f::Identity();
    svd_v = Matrix3f::Identity();
    svd_s = Vector3f(1, 1, 1);
    v_PIC = Vector3f(0, 0, 0);
    v_FLIP = Vector3f(0, 0, 0);
}

SnowParticle::~SnowParticle()
{
}

void SnowParticle::update_pos()
{
    // Simple euler integration
    position += deltaT * velocity;
}

void SnowParticle::update_velocity() // step 8
{
    velocity = (1 - m->alpha) * v_PIC + m->alpha * v_FLIP;
}


void SnowParticle::update_deform_gradient()
{
    deform_elastic_grad = (Matrix3f::Identity() + deltaT * v_grad) * deform_elastic_grad;
    deform_grad = deform_elastic_grad * deform_plastic_grad;
    JacobiSVD svd(deform_elastic_grad, ComputeFullU | ComputeFullV);
    svd_s = svd.singularValues();
    svd_u = svd.matrixU();
    svd_v = svd.matrixV();
    // clamp the singular value within the limit
    for (int i = 0; i < 3; i++)
    {
        svd_s[i] = clamp(1 - m->criticalStress, m->criticalStretch, svd_s[i]);
    }
    // compute the elastic and plastic gradient
    for (int i = 0; i < 3; i++)
    {
        svd_s[i] = 1. / svd_s[i];
    }
    deform_elastic_grad = svd_u * svd_s.asDiagonal() * svd_v.transpose();
    deform_plastic_grad = svd_v * svd_s.asDiagonal() * svd_u.transpose() * deform_grad;
}

const Matrix3f SnowParticle::energy_derivative()
{
    // Adjust lame parameters to account for m->hardening
    float Jp = deform_plastic_grad.determinant();
    float harden = exp(m->hardening * (1. - Jp));
    float Je = svd_s(0) * svd_s(1) * svd_s(2);
    // Compute co-rotational term
    Matrix3f deform_elastic_polar_r = svd_u * svd_v.transpose();
    Matrix3f temp = 2. * m->mu *
                    (deform_elastic_grad - deform_elastic_polar_r) *
                    deform_elastic_grad.transpose();
    // Add in the primary contour term
    temp += m->lambda * Je * (Je - 1.) * Matrix3f::Identity();
    
    return volume * harden * temp;
}

// Matrix3f SnowParticle::volume_cauchy_stress()
// {
//     float Jp = deform_plastic_grad.determinant();
//     float Je = deform_elastic_grad.determinant();
//     float J = deform_grad.determinant();
//     float harden = std::exp(m->hardening * (1 - Jp));
//     float mu_grad = m->mu * harden;
//     float lambda_grad = m->lambda * harden;
//     Matrix3f deform_elastic_polar_r = svd_u * svd_v.transpose();
//     return 2.0 * volume * mu_grad * (deform_elastic_grad - deform_elastic_polar_r) * deform_elastic_grad.transpose() + MatrixXf::Constant(3, 3, (lambda_grad * (Je - 1.0) * Je * volume));
// }

/* SPS */

SnowParticleSet::SnowParticleSet() : particles()
{
}

SnowParticleSet::~SnowParticleSet()
{
    for (auto& oneParticle : particles)
    {
        delete (oneParticle);
    }
}

void SnowParticleSet::addParticle(SnowParticle* sp)
{
    particles.push_back(sp);
}

void SnowParticleSet::addParticle(const Vector3f& pos, const Vector3f& vel,
                                  const float Mass, SnowParticleMaterial* m)
{
    // std::cout << " sp mass is " << Mass << std::endl;
    SnowParticle* sp = new SnowParticle(pos, vel, Mass, m);
    assert(sp);
    particles.push_back(sp);
}

void SnowParticleSet::addParticlesInShape(Shape* s, const Vector3f& vel,
                                           SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    // std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        float total_mass = s->getVolume() * m->initialDensity;
        float mass_per_particle = total_mass / (float)temp;

        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, vel, mass_per_particle, m);
        }
    }
}

void SnowParticleSet::appendSet(SnowParticleSet& anotherSet)
{
    for (SnowParticle* p : anotherSet.particles)
    {
        particles.push_back(p);
    }
    anotherSet.particles.clear();
    // the other set will be cleared in the end
    // this is to avoid deleting particle* error when the latter is destroyed
}

void SnowParticleSet::CreateMirror(const SnowParticleSet& anotherSet, float a,
                                   float b, float c, float d, const Vector3f p)
{
    for (SnowParticle* sp : anotherSet.particles)
    {
        Vector3f pos = sp->position;
        Vector3f vel = sp->velocity;
        Vector3f relPos = p - pos;
        Vector3f relPosNormalToPlane =
            relPos.dot(Vector3f(a, b, c)) * Vector3f(a, b, c);
        Vector3f newPos = 2. * relPosNormalToPlane + pos;

        Vector3f relVelTOPlane = vel.dot(Vector3f(a, b, c)) * Vector3f(a, b, c);
        Vector3f relVelTangent = vel - relVelTOPlane;
        Vector3f newVel = relVelTangent - relVelTOPlane;

        addParticle(newPos, newVel, sp->mass, sp->m);
    }
}

void SnowParticleSet::addParticlesInShape(Shape* s, SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    Vector3f vel(0, 0, 0);
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        float total_mass = s->getVolume() * m->initialDensity;
        float mass_per_particle = total_mass / (float)temp;

        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, vel, mass_per_particle, m);
        }
    }
}

void SnowParticleSet::update()
{
    max_velocity = 0;
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i]->update_pos();
        //particles[i]->updatePureElasticGradient();
        particles[i]->update_deform_gradient();
        // Update max velocity, if needed
        float vel = particles[i]->velocity.norm();
        if (vel > max_velocity) max_velocity = vel;
    }
}

