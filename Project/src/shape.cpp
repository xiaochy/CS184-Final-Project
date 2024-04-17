#include "shape.h"
#include "bounds3.h"

void Shape::addPoint(float x, float y, float z)
{
    Vector3f new_pt = Vector3f(x, y, z);
    points.push_back(new_pt);
}

// void Shape::contains(float x, float y, float z)
// {
// }

// https://www.mathopenref.com/coordpolygonarea2.html
float Shape::area()
{
    int num_points = points.size();
    float area = 0;
    // j is the previous vertex to i
    int j = num_points - 1;

    for (int i = 0; i < num_points; i++)
    {
        Vector2f pt = points[i];
        Vector2f pt_prev = points[j];
        area += (pt_prev.x() + pt.x()) * (pt_prev.y() - pt.y());
        j = i;
    }
    return fabs(area / 2); // return the absolute value
}

float Shape::volume()
{
    float len = sqrt(area());
    return len * len * len;
}

Bounds3 Shape::bounds()
{
    float bounds[6];
    // X-bound
    bounds[0] = points[0][0];
    bounds[1] = bounds[0];
    // Y-bound
    bounds[2] = points[0][1];
    bounds[3] = bounds[2];
    // Z-bound
    bounds[4] = points[0][2];
    bounds[5] = bounds[4];
    for (int i = 0, len = points.size(); i < len; i++)
    {
        Vector3f &p = points[i];
        // updata X-bound
        if (p[0] < bounds[0])
        {
            bounds[0] = p[0];
        }
        else if (p[0] > bounds[1])
            bounds[1] = p[0];
        // update Y-bound
        if (p[1] < bounds[2])
            bounds[2] = p[1];
        else if (p[1] > bounds[3])
            bounds[3] = p[1];
        // update Z-bound
        if (p[2] < bounds[4])
            bounds[4] = p[2];
        else if (p[2] > bounds[5])
            bounds[5] = p[2];
    }
    Vector3f pmin = Vector3f(bounds[0], bounds[2], bounds[4]);
    Vector3f pmax = Vector3f(bounds[1], bounds[3], bounds[5]);
    return Bounds3(pmin, pmax);
}
