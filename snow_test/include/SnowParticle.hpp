#ifndef SNOWSIM_SNOWPARTICLES
#define SNOWSIM_SNOWPARTICLES

#include "Shape.hpp"
#include <iostream>
using namespace Eigen;

class SnowParticleMaterial
{
   private:
   public:
    float initialDensity = 400.;
    // line number density when generating particles
    // initial value 138 = 0.0072 diameter
    int lNumDensity = 138;
    //float criticalStress = 1. - 2.5e-2;
    float criticalStress = 1. - 2.5e-1;
    //float criticalStretch = 1. + 7.5e-3;
    float criticalStretch = 1. + 7.5e-2;
    float hardening = 10.;
    float youngsModule = 1.4e5;
    float PoissonsRatio = .2;
    // Collision stickiness (lower = stickier)
    float sticky = .9;
    // Lame parameters
    float mu;
    float lambda;
    // change begin
    float alpha = .95;
    float beta;
    // change end

    // etc

    SnowParticleMaterial();
    ~SnowParticleMaterial();
};

class SnowParticle
{
   private:
   public:
    // change begin
    float volume, mass, density;
    //Vector3f old_pos, new_pos;
    //Vector3f old_v, new_v;
    //Vector3f velocity;
    Vector3f position;
    Vector3f velocity;
    Vector3f v_PIC;
    Vector3f v_FLIP;
    Matrix3f v_grad;
    Matrix3f svd_u, svd_v;
    Vector3f svd_s;
    Matrix3f deform_grad;
    Matrix3f deform_elastic_grad;
    Matrix3f deform_plastic_grad;
    Vector3f weight_gradient[64];
    float weights[64];
    SnowParticleMaterial *m;
    // change end

    SnowParticle();
    SnowParticle(const Vector3f& pos, const Vector3f& vel, const float mass,
                 SnowParticleMaterial* material);
    ~SnowParticle();

    // change begin
    void update_pos();
    void update_velocity();
    void update_deform_gradient();
    Matrix3f volume_cauchy_stress();
    // change end

    // void updatePos();
    // // Update deformation gradient
    void updatePureElasticGradient();
    void updateCombinedPElasticGradient();
    // // Compute stress tensor
    const Matrix3f energyDerivative();

    // Computes stress force delta, for implicit velocity update
    const Vector3f deltaForce(const Vector2f& u, const Vector2f& weight_grad);
};

class SnowParticleSet
{
   private:
   public:
    std::vector<SnowParticle*> particles;
    float maxVelocity;

    SnowParticleSet();
    ~SnowParticleSet();
    void addParticle(SnowParticle* sp);
    void addParticle(const Vector3f& pos, const Vector3f& vel, const float Mass,
                     SnowParticleMaterial* m);
    // TODO can consider initial rotation in a snow shape
    void addParticlesInAShape(Shape* s, SnowParticleMaterial* m);
    void addParticlesInAShape(Shape* s, const Vector3f& vel,
                              SnowParticleMaterial* m);
    void appendSet(SnowParticleSet& anotherSet);
    void CreateMirror(const SnowParticleSet& anotherSet, float a, float b,
                      float c, float d, const Vector3f p);
    // inline SnowParticleSet unionSet(const SnowParticleSet& set1,
    //                                 const SnowParticleSet& set2);
    // inline SnowParticleSet unionSet(const std::vector<SnowParticleSet>&
    // sets);
    void update();
};

#endif  // SNOWSIM_SNOWPARTICLES
