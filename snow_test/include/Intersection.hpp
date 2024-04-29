//
// Created by LEI XU on 5/16/19.
//

#ifndef RAYTRACING_INTERSECTION_H
#define RAYTRACING_INTERSECTION_H
// #include "Shape.hpp"
#include <eigen3/Eigen/Eigen>
#include <limits>

using namespace Eigen;

class Shape;

class Intersection
{
public:
    Intersection()
    {
        happened = false;
        coords = Vector3f();
        normal = Vector3f();
        distance = std::numeric_limits<double>::max();
        shape = nullptr;
    }
    bool happened;
    Vector3f coords;
    Vector3f normal;
    double distance;
    Shape* shape;
};
#endif  // RAYTRACING_INTERSECTION_H
