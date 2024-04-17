#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "shape.hpp"

class Triangle : public Shape
{
   public:
    Vector3f v0, v1, v2;  // vertices A, B, C, counter-clockwise order
    Vector3f e1, e2;      // 2 edges v1-v0, v2-v0;
    Vector3f normal;
    float area;

    Triangle(Vector3f _v0, Vector3f _v1, Vector3f _v2)
        : v0(_v0), v1(_v1), v2(_v2)
    {
        e1 = v1 - v0;
        e2 = v2 - v0;
        normal = e1.cross(e2).normalized();
        area = e1.cross(e2).norm() * 0.5f;
    }
    Bounds3 getBounds()
    {
        return Union(Bounds3(v0, v1), v2);
    }
    float getArea()
    {
        return area;
    }
    float getVolume()
    {
        return 0.;
    }
    Vector3f centroid()
    {
        return (v0 + v1 + v2) / 3.;
    }
};