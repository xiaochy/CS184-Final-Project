#ifndef PARTICLE_H
#define PARTICLE_H

#include <cmath>
#include "CGL/Vector3D.h"
#include "Matrix2f.h"
#include "SimConstants.h"
#include "CGL/matrix3x3.h"

class Particle
{
public:
    float volume, mass, density;
    Vector3D position, velocity;
    Vector3D velocity_star;
    Matrix3D velocity_gradeint;
    // // Lame parameters (_s denotes starting configuration)
    // float lambda,
    //     mu;
    // Deformation gradients
    Matrix3x3 deform_gradient_elastic, deform_gradient_plastic;
    // Cached SVD's for elastic deformation gradient
    Matrix2f svd_w,
        svd_v;
    Vector2f svd_e;
    // Cached polar decomposition
    Matrix2f polar_r, polar_s;
    // Grid interpolation weights
    Vector3D grid_position;
    Vector3D weight_gradient[64];
    float weights[64];

    Particle();
    Particle(const Vector2f &pos, const Vector2f &vel, float mass, float lambda, float mu);
    virtual ~Particle();

    // Update position, based on velocity
    void updatePos();
    // Update deformation gradient
    void updateGradient();
    void applyPlasticity();
    // Compute stress tensor
    const Matrix2f energyDerivative();

    // Computes stress force delta, for implicit velocity update
    const Vector2f deltaForce(const Vector2f &u, const Vector2f &weight_grad);

private:
    void computeWeights();
};

#endif
