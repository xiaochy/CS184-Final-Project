#include <algorithm>
#include "Particle.h"
#include "global.h“
#include <eigen3/Eigen/Eigen>

using namespace Eigen;

void Particle::update_pos()
{
    new_pos = old_pos + deltaT * new_v;
}

void Particle::update_velocity()
{
    new_v = (1 - m->alpha) * v_PIC + m->alpha * v_FLIP;
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

Matrix3f Particle::cauchy_stress()
{
    float Jp = deform_plastic_grad.determinant();
    float Je = deform_elastic_grad.determinant();
    float J = deform_grad.determinant();
    float harden = std::exp(m->hardening * (1 - Jp));
    float mu_grad = m->mu * harden;
    float lambda_grad = m->lambda * harden;
    Matrix3f deform_elastic_polar_r = svd_u * svd_v.transpose();
    return 2.0 * mu_grad / J * (deform_elastic_grad - deform_elastic_polar_r) * deform_elastic_grad.transpose() + MatrixXf::Constant(3, 3, (lambda_grad * (Je - 1.0) * Je / J));
}

// // Collision stickiness (lower = stickier)
// float sticky = .9;

// void Particle::computeWeights()
// {
//     float x_offset = position.x - floor(position.x);

//     float wx[4] = {
//         bspline(x_offset - 1),
//         bspline(x_offset),
//         bspline(x_offset + 1),
//         bspline(x_offset + 2)};

//     float y_offset = position.y - floor(position.y);

//     float wy[4] = {
//         bspline(y_offset - 1),
//         bspline(y_offset),
//         bspline(y_offset + 1),
//         bspline(y_offset + 2)};

//     float z_offset = position.z - floor(position.z);

//     float wz[4] = {
//         bspline(z_offset - 1),
//         bspline(z_offset),
//         bspline(z_offset + 1),
//         bspline(z_offset + 2)};

//     for (int i = 0; i < 4; i++)
//     {
//         for (int j = 0; j < 4; j++)
//         {
//             for (int k = 0; k < 4; k++)
//             {
//                 weights[i * 16 + j * 4 + k] = wx[i] * wy[j] * wz[k];
//             }
//         }
//     }
// }

// const Matrix3D Particle::energyDerivative()
// {
//     // Adjust lame parameters to account for m->hardening
//     float harden =
//         exp(m->hardening * (1. - deformationGradientPlastic.determinant()));
//     float Je = SVDS.x() * SVDS.y() * SVDS.z();
//     // This is the co-rotational term
//     Matrix3D temp = 2. * m->mu *
//                     (deformationGradientElastic - SVDU * SVDV.transpose()) *
//                     deformationGradientElastic.transpose();
//     // Add in the primary contour term
//     temp += m->lambda * Je * (Je - 1.) * Matrix3D::Identity();

//     return volume * harden * temp;
// }
