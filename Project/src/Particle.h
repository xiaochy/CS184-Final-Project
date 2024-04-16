#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"

// using namespace CGL;

class SnowParticleMaterial
{
private:
public:
    float initialDensity = 400.;
    // line number density when generating particles
    // initial value 138 = 0.0072 diameter
    int lNumDensity = 138;
    float critical_compress = 2.5e-2;
    float criticalStretch = 7.5e-3;
    float hardening = 10.;
    float youngsModule = 1.4e5;
    float PoissonsRatio = .2;
    // Collision stickiness (lower = stickier)
    float sticky = .9;
    // Lame parameters
    float mu;
    float lambda;
    float alpha;
    float beta;

    SnowParticleMaterial();
    ~SnowParticleMaterial();
};

class Particle
{
public:
    float volume, mass, density;
    Vector3f old_pos, new_pos;
    Vector3f old_v, new_v;
    Vector3f v_PIC;
    Vector3f v_FLIP;
    Matrix3f v_grad;
    Matrix3f svd_u, svd_s, svd_v;
    Matrix3f deform_grad;
    Matrix3f deform_elastic_grad;
    Matrix3f deform_plastic_grad;
    Vector3f weight_gradient[64];
    float weights[64];
    SnowParticleMaterial *m;

    Particle();
    virtual ~Particle();
    void update_pos();
    void update_velocity();
    void update_deform_gradient();
    Matrix3f volume_cauchy_stress();
};

#endif
