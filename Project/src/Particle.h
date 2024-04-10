#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"

using namespace CGL;

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

    // etc

    SnowParticleMaterial();
    ~SnowParticleMaterial();
};

class Particle
{
public:
    float volume, mass, density;
    Vector3D position, velocity;
    Vector3D velocity_star;
    Matrix3x3 velocity_gradeint;
    Matrix3x3 deform_gradient_elastic, deform_gradient_plastic;
    // Cached SVD's for elastic deformation gradient
    Matrix3x3 svd_w, svd_v;
    Vector2D svd_e;
    // Cached polar decomposition
    Matrix3x3 polar_r, polar_s;
    // Grid interpolation weights
    Vector3D grid_position;
    Vector3D weight_gradient[64];
    float weights[64];
    SnowParticleMaterial *m;
    Matrix3x3 deformationGradientElastic;
    Matrix3D deformationGradientPlastic;

    Particle();
    Particle(const Vector3D &pos, const Vector3D &vel, float mass, float lambda, float mu);
    virtual ~Particle();

    // Update position, based on velocity
    void updatePos();
    // Update deformation gradient
    void updateGradient();
    void applyPlasticity();
    // Compute stress tensor
    const Matrix3x3 energyDerivative();

    // Computes stress force delta, for implicit velocity update
    const Vector3D deltaForce(const Vector3D &u, const Vector3D &weight_grad);

private:
    void computeWeights();
};

#endif
