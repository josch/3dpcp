#include "slam6d/point.h"
#include <vector>
#include <iostream>

#ifndef __POINTNEIGHBOR_H__
#define __POINTNEIGHBOR_H__

class PointNeighbor : Point {
  public:
    /// x coordinate in 3D space
    double x;
    /// y coordinate in 3D space
    double y;
    /// z coordinate in 3D space
    double z;
    /// additional information about the point, e.g., semantic
    ///  also used in veloscan for distiuguish moving or static
    int type;

    /////////////////////////for veloslam/////////////////////////////
    double rad;
    ///    tang in  cylindrical coordinates for veloscan
    double tan_theta;
    // point id in points for veloscan , you can use it find point.
    long point_id;
    /////////////////////////for veloslam/////////////////////////////

    // color information of the point between 0 and 255
    // rgb
    unsigned char rgb[3];

    float reflectance;
    float temperature;
    float amplitude;
    float deviation;

    std::vector<Point> neighbors;

    PointNeighbor() : Point() {
        neighbors.clear();
    }    
    
    PointNeighbor(const Point& p, const std::vector<Point>& neighbors) : Point(p) {
        this->neighbors = neighbors;
    }

    PointNeighbor(const double _x, const double _y, const double _z) : Point(_x, _y, _z) {
        neighbors.clear();
    }

    PointNeighbor(const double *p) : Point(p) {
        neighbors.clear();
    }
};

#endif
