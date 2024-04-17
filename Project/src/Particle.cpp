#include <algorithm>
#include "Particle.hpp"
#include "global.hpp"
#include <eigen3/Eigen/Eigen>

using namespace Eigen;

Particle::Particle() : m(nullptr)
{
    // volume mass and density will be calculated initially
    // so no need to assign now
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
    // polarR = Matrix3f::Identity();
    // polarTheta = Matrix3f::Identity();
    // polarPhi = Matrix3f::Identity();
}

Particle::Particle(const Vector3f& pos, const Vector3f& vel,
                           const float Mass, SnowParticleMaterial* material)
    : position(pos), velocity(vel), mass(Mass), m(material)
{
    // volume mass and density will be calculated initially
    // so no need to assign now
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
    // polarR = Matrix3f::Identity();
    // polarTheta = Matrix3f::Identity();
    // polarPhi = Matrix3f::Identity();

}

Particle::~Particle()
{
}

SnowParticleMaterial::SnowParticleMaterial()
{
    mu = youngsModule / (2. + 2. * PoissonsRatio);
    lambda = youngsModule * PoissonsRatio /
             ((1. + PoissonsRatio) * (1. - 2. * PoissonsRatio));
}

SnowParticleMaterial::~SnowParticleMaterial()
{
}

void Particle::update_pos()
{
    new_pos = old_pos + deltaT * new_v;
    old_pos = new_pos;
}

void Particle::update_velocity()
{
    new_v = (1 - m->alpha) * v_PIC + m->alpha * v_FLIP;
    old_v = new_v;
}

void Particle::update_deform_gradient()
{
    deform_elastic_grad = (Matrix3f::Identity() + deltaT * v_grad) * deform_elastic_grad;
    deform_grad = deform_elastic_grad * deform_plastic_grad;
    JacobiSVD svd(deform_elastic_grad, ComputeFullU | ComputeFullV);
    MatrixXf s = svd.singularValues();
    svd_u = svd.matrixU();
    svd_v = svd.matrixV();
    // clamp the singular value within the limit
    for (int i = 0; i < 3; i++)
    {
        s[i] = clamp(1 - m->critical_compress, m->criticalStretch, s[i]);
    }
    // compute the elastic and plastic gradient
    svd_s = s.asDiagonal();
    for (int i = 0; i < 3; i++)
    {
        s[i] = 1. / s[i];
    }
    deform_elastic_grad = svd_u * svd_s * svd_v.transpose();
    deform_plastic_grad = svd_v * s.asDiagonal() * svd_u.transpose() * deform_grad;
}

Matrix3f Particle::volume_cauchy_stress()
{
    float Jp = deform_plastic_grad.determinant();
    float Je = deform_elastic_grad.determinant();
    float J = deform_grad.determinant();
    float harden = std::exp(m->hardening * (1 - Jp));
    float mu_grad = m->mu * harden;
    float lambda_grad = m->lambda * harden;
    Matrix3f deform_elastic_polar_r = svd_u * svd_v.transpose();
    return 2.0 * volume * mu_grad * (deform_elastic_grad - deform_elastic_polar_r) * deform_elastic_grad.transpose() + MatrixXf::Constant(3, 3, (lambda_grad * (Je - 1.0) * Je * volume));
}

void SnowParticleSet::addParticle(Particle* sp)
{
    particles.push_back(sp);
}

void SnowParticleSet::addParticle(const Vector3f& pos, const Vector3f& vel,
                                  const float Mass, SnowParticleMaterial* m)
{
    Particle* sp = new Particle(pos, vel, Mass, m);
    particles.push_back(sp);
}

void SnowParticleSet::addParticlesInAShape(Shape* s, const Vector3f& vel,
                                           SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    if (temp > 0)
    {
        float totMass = s->getVolume() * m->initialDensity;
        float massPerP = totMass / (float)temp;

        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, vel, massPerP, m);
        }
    }
}

void SnowParticleSet::addParticlesInAShape(Shape* s, SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    Vector3f vel(0, 0, 0);
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    if (temp > 0)
    {
        float totMass = s->getVolume() * m->initialDensity;
        float massPerP = totMass / (float)temp;

        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, vel, massPerP, m);
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

void SnowParticleSet::update()
{
    maxVelocity = 0;
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i]->updatePos();
        particles[i]->updatePureElasticGradient();
        particles[i]->updateCombinedPElasticGradient();
        // Update max velocity, if needed
        float vel = particles[i]->velocity.norm();
        if (vel > maxVelocity) maxVelocity = vel;
    }
}
