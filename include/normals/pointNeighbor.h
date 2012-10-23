#include "slam6d/point.h"
#include <vector>
#include <iostream>

#ifndef __POINTNEIGHBOR_H__
#define __POINTNEIGHBOR_H__

class PointNeighbor {
  public:
    Point point;
    std::vector<Point> neighbors;

    float range_neighbors[3][3];
    int range_image_row;
    int range_image_col;

    PointNeighbor() {
        point.x = 0.; 
        point.y = 0.;
        point.z = 0.;
        neighbors.clear();
    }    
    
    PointNeighbor(const Point& point, const std::vector<Point>& neighbors) {
        this->point = point;  
        this->neighbors = neighbors;
    }

    PointNeighbor(const double _x, const double _y, const double _z) {
        point.x = _x;
        point.y = _y;
        point.z = _z;
        neighbors.clear();
    }

    PointNeighbor(const double *p) {
        point.x = p[0];
        point.y = p[1];
        point.z = p[2];
        neighbors.clear();
    }
};

#endif
