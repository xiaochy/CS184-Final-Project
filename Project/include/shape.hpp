#ifndef SNOWSIM_SHAPE
#define SNOWSIM_SHAPE
#include "Bounds3.hpp"

using namespace Eigen;

class Shape
{
public:
    std::vector<Vector3f> points;

    Shape();
    Shape(const Shape &orig);
    virtual ~Shape();

    void addPoint(float x, float y, float z);

    // bool contains(float x, float y, float z);

    float area();

    float volume();

    Bounds3 bounds();

    void draw();
};
#endif