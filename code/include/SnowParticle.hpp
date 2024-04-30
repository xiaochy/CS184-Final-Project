#ifndef SNOWSIM_SNOWPARTICLES
#define SNOWSIM_SNOWPARTICLES

#include "Shape.hpp"
#include "Sphere.hpp"

using namespace Eigen;

class SnowParticleMaterial
{
   private:
   public:
    float initialDensity = 400.;
    // line number density when generating particles
    // initial value 138 = 0.0072 diameter
    int lNumDensity = 138;
    float criticalStress = 1. - 2.5e-2;
    float criticalStretch = 1. + 7.5e-3;
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
    float volume, mass, density;
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

    SnowParticle();
    SnowParticle(const Vector3f& pos, const Vector3f& vel, const float mass,
                 SnowParticleMaterial* material);
    ~SnowParticle();
    void update_pos();
    void update_velocity();
    void update_deform_gradient();
    // Update deformation gradient
    void updatePureElasticGradient();
    void updateCombinedPElasticGradient();
    // Compute stress tensor
    const Matrix3f energy_derivative();
    Matrix3f volume_cauchy_stress();

};

class SnowParticleSet
{
   private:
   public:
    std::vector<SnowParticle*> particles;
    float max_velocity;

    SnowParticleSet();
    ~SnowParticleSet();
    void addParticle(SnowParticle* sp);
    void addParticle(const Vector3f& pos, const Vector3f& vel, const float Mass,
                     SnowParticleMaterial* m);

    void addParticlesInShape(Shape* s, SnowParticleMaterial* m);
    void addParticlesInShape(Shape* s, const Vector3f& vel,
                              SnowParticleMaterial* m);
    void addParticlesInSphere(Sphere*s, const Vector3f& vel, SnowParticleMaterial*m);
    void appendSet(SnowParticleSet& anotherSet);
    void CreateMirror(const SnowParticleSet& anotherSet, float a, float b,
                      float c, float d, const Vector3f p);
    void update();
};

#endif  // SNOWSIM_SNOWPARTICLES
