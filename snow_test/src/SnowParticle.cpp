#include "SnowParticle.hpp"

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
    // so no need to assign now
    volume = 0;
    mass = 0;
    density = 0;
    position = Vector3f(0, 0, 0);
    velocity = Vector3f(0, 0, 0);
    // velocityGradient = Matrix3f::Zero();
    v_grad = Matrix3f::Zero();
    // deformationGradientElastic = Matrix3f::Identity();
    // deformationGradientPlastic = Matrix3f::Identity();
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

SnowParticle::SnowParticle(const Vector3f& pos, const Vector3f& vel,
                           const float Mass, SnowParticleMaterial* material)
    : position(pos), velocity(vel), mass(Mass), m(material)
{
    // volume mass and density will be calculated initially
    // so no need to assign now
    volume = 0;
    // mass = 0;
    density = 0;
    // velocityGradient = Matrix3f::Zero();
    // deformationGradientElastic = Matrix3f::Identity();
    // deformationGradientPlastic = Matrix3f::Identity();
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

    // std::cout << " sp mass is " << mass << std::endl;
}

SnowParticle::~SnowParticle()
{
}

void SnowParticle::update_pos()
{
    // Simple euler integration
    position += deltaT * velocity;
    new_pos = old_pos + deltaT * new_v;
    old_pos = new_pos;
}

void SnowParticle::update_velocity() // step 8
{
    new_v = (1 - m->alpha) * v_PIC + m->alpha * v_FLIP;
    old_v = new_v;
    //velocity = new_v;
}


void SnowParticle::update_deform_gradient()
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
        s[i] = clamp(1 - m->criticalStress, m->criticalStretch, s[i]);
    }
    // compute the elastic and plastic gradient
    svd_s = s.asDiagonal();
    for (int i = 0; i < 3; i++)
    {
        s[i] = 1. / s[i];
    }
    deform_elastic_grad = svd_u * svd_s * svd_v.transpose();
    deform_plastic_grad = svd_v * s.asDiagonal() * svd_u.transpose() * deform_grad;
    // Matrix3f VDivideSVDS(svd_v), UTimeSVDS(svd_u);
    // VDivideSVDS.col(0) /= svd_s[0];
    // VDivideSVDS.col(1) /= svd_s[1];
    // VDivideSVDS.col(2) /= svd_s[2];
    // UTimeSVDS.col(0) *= svd_s[0];
    // UTimeSVDS.col(1) *= svd_s[1];
    // UTimeSVDS.col(2) *= svd_s[2];
    // deform_plastic_grad =
    //     VDivideSVDS * svd_u.transpose() * deform_grad;
    // deform_elastic_grad = UTimeSVDS * svd_v.transpose();
}

/* energyDerivative() */
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

const Matrix3f SnowParticle::energyDerivative()
{
    // Adjust lame parameters to account for m->hardening
    float harden =
        exp(m->hardening * (1. - deform_plastic_grad.determinant()));
    float Je = svd_s.x() * svd_s.y() * svd_s.z();
    // This is the co-rotational term
    Matrix3f temp = 2. * m->mu *
                    (deform_elastic_grad - svd_u * svd_v.transpose()) *
                    deform_elastic_grad.transpose();
    // Add in the primary contour term
    temp += m->lambda * Je * (Je - 1.) * Matrix3f::Identity();
    
    return volume * harden * temp;
}


// #if ENABLE_IMPLICITF
const Vector3f SnowParticle::deltaForce(const Vector2f& u,
                                        const Vector2f& weight_grad)
{
    return Vector3f(1., 1., 1.);
    // TODO

    // // For detailed explanation, check out the implicit math pdf for details
    // // Before we do the force calculation, we need deltaF, deltaR, and
    // // delta(JF^-T)

    // // Finds delta(Fe), where Fe is the elastic deformation gradient
    // // Probably can optimize this expression with parentheses...
    // Matrix2f del_elastic =
    //     deltaT * u.outer_product(weight_grad) * deformationGradientElastic;

    // // Check to make sure we should do these calculations?
    // if (del_elastic[0][0] < MATRIX_EPSILON &&
    //     del_elastic[0][1] < MATRIX_EPSILON &&
    //     del_elastic[1][0] < MATRIX_EPSILON &&
    //     del_elastic[1][1] < MATRIX_EPSILON)
    //     return Vector2f(0);

    // // Compute R^T*dF - dF^TR
    // // It is skew symmetric, so we only need to compute one value (three for
    // 3D) float y =
    //     (polar_r[0][0] * del_elastic[1][0] +
    //      polar_r[1][0] * del_elastic[1][1]) -
    //     (polar_r[0][1] * del_elastic[0][0] + polar_r[1][1] *
    //     del_elastic[0][1]);
    // // Next we need to compute MS + SM, where S is the hermitian matrix
    // // (symmetric for real valued matrices) of the polar decomposition and M
    // is
    // // (R^T*dR); This is equal to the matrix we just found (R^T*dF ...), so
    // we
    // // set them equal to eachother Since everything is so symmetric, we get a
    // // nice system of linear equations once we m->multiply everything out.
    // (see
    // // pdf for details) In the case of 2D, we only need to solve for one
    // // variable (three for 3D)
    // float x = y / (polar_s[0][0] + polar_s[1][1]);
    // // Final computation is deltaR = R*(R^T*dR)
    // Matrix2f del_rotate = Matrix2f(-polar_r[1][0] * x, polar_r[0][0] * x,
    //                                -polar_r[1][1] * x, polar_r[0][1] * x);

    // // We need the cofactor matrix of F, JF^-T
    // Matrix2f cofactor = deformationGradientElastic.cofactor();

    // // The last matrix we need is delta(JF^-T)
    // // Instead of doing the complicated matrix derivative described in the
    // paper
    // // we can just take the derivative of each individual entry in JF^-T;
    // JF^-T
    // // is the cofactor matrix of F, so we'll just hardcode the whole thing
    // For
    // // example, entry [0][0] for a 3x3 matrix is 	cofactor = e*i - f*h
    // // derivative
    // //= (e*Di + De*i) - (f*Dh + Df*h) 	where the derivatives (capital D)
    // come
    // // from our precomputed delta(F) In the case of 2D, this turns out to be
    // // just
    // // the cofactor of delta(F) For 3D, it will not be so simple
    // Matrix2f del_cofactor = del_elastic.cofactor();

    // // Calculate "A" as given by the paper
    // // Co-rotational term
    // Matrix2f Ap = del_elastic - del_rotate;
    // Ap *= 2 * m->mu;
    // // Primary contour term
    // cofactor *= cofactor.frobeniusInnerProduct(del_elastic);
    // del_cofactor *= (deformationGradientElastic.determinant() - 1);
    // cofactor += del_cofactor;
    // cofactor *= m->lambda;
    // Ap += cofactor;

    // // Put it all together
    // // Parentheses are important; M*M*V is slower than M*(M*V)
    // return volume *
    //        (Ap * (deformationGradientElastic.transpose() * weight_grad));
}



// #endif

SnowParticleSet::SnowParticleSet() : particles()
{
}

SnowParticleSet::~SnowParticleSet()
{
    // TODO check if this destroy way is okay
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
    particles.push_back(sp);
}

void SnowParticleSet::addParticlesInAShape(Shape* s, const Vector3f& vel,
                                           SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    // std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        float totMass = s->getVolume() * m->initialDensity;
        float massPerP = totMass / (float)temp;
        // std::cout << " totMass is " << totMass << std::endl;
        // std::cout << " massPerP is " << massPerP << std::endl;

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

void SnowParticleSet::addParticlesInAShape(Shape* s, SnowParticleMaterial* m)
{
    std::vector<Vector3f> tempPos;
    Vector3f vel(0, 0, 0);
    int temp = s->generateParticlesInside(m->lNumDensity, tempPos);
    // std::cout << " size is " << temp << std::endl;
    if (temp > 0)
    {
        float totMass = s->getVolume() * m->initialDensity;
        float massPerP = totMass / (float)temp;
        // std::cout << " totMass is " << totMass << std::endl;
        // std::cout << " massPerP is " << massPerP << std::endl;

        for (const auto& onePos : tempPos)
        {
            addParticle(onePos, vel, massPerP, m);
        }
    }
}

void SnowParticleSet::update()
{
    maxVelocity = 0;
    for (int i = 0; i < particles.size(); i++)
    {
        particles[i]->update_pos();
        //particles[i]->updatePureElasticGradient();
        particles[i]->update_deform_gradient();
        // Update max velocity, if needed
        float vel = particles[i]->velocity.norm();
        if (vel > maxVelocity) maxVelocity = vel;
    }
}

