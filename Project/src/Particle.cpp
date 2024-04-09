#include <algorithm>
#include "Particle.h"

void Particle::computeWeights()
{
    float x_offset = position.x - floor(position.x);

    float wx[4] = {
        bspline(x_offset - 1),
        bspline(x_offset),
        bspline(x_offset + 1),
        bspline(x_offset + 2)};

    float y_offset = position.y - floor(position.y);

    float wy[4] = {
        bspline(y_offset - 1),
        bspline(y_offset),
        bspline(y_offset + 1),
        bspline(y_offset + 2)};

    float z_offset = position.z - floor(position.z);

    float wz[4] = {
        bspline(z_offset - 1),
        bspline(z_offset),
        bspline(z_offset + 1),
        bspline(z_offset + 2)};

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                weights[i * 16 + j * 4 + k] = wx[i] * wy[j] * wz[k];
            }
        }
    }
}

const Matrix3D Particle::energyDerivative()
{
    // Adjust lame parameters to account for m->hardening
    float harden =
        exp(m->hardening * (1. - deformationGradientPlastic.determinant()));
    float Je = SVDS.x() * SVDS.y() * SVDS.z();
    // This is the co-rotational term
    Matrix3D temp = 2. * m->mu *
                    (deformationGradientElastic - SVDU * SVDV.transpose()) *
                    deformationGradientElastic.transpose();
    // Add in the primary contour term
    temp += m->lambda * Je * (Je - 1.) * Matrix3D::Identity();

    return volume * harden * temp;
}