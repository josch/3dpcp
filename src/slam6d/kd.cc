/*
 * kd implementation
 *
 * Copyright (C) Andreas Nuechter, Kai Lingemann, Thomas Escher
 *
 * Released under the GPL version 3.
 *
 */

/** @file 
 *  @brief An optimized k-d tree implementation
 *  @author Andreas Nuechter. Institute of Computer Science, University of Osnabrueck, Germany.
 *  @author Kai Lingemann. Institute of Computer Science, University of Osnabrueck, Germany.
 *  @author Thomas Escher Institute of Computer Science, University of Osnabrueck, Germany.
 */

#ifdef _MSC_VER
#define  _USE_MATH_DEFINES
#endif

#include "slam6d/kd.h"
#include "slam6d/globals.icc"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
#include <algorithm>
using std::swap;
#include <cmath>
#include <cstring>

// KDtree class static variables
template<class PointData, class AccessorData, class AccessorFunc>
KDParams KDTreeImpl<PointData, AccessorData, AccessorFunc>::params[MAX_OPENMP_NUM_THREADS];

/**
 * struct used for sorting the neighbor arrays
 * comparing by distance to the origin point
 */

struct ComparePoints {
    double *origin;
    bool operator() (double* i, double* j) { 
        return (Dist2(i, origin) < Dist2(j, origin));
    }
} distance_comparison;


/**
 * Constructor
 *
 * Create a KD tree from the points pointed to by the array pts
 *
 * @param pts 3D array of points
 * @param n number of points
 */
KDtree::KDtree(double **pts, int n)
{
    create(Void(), pts, n);
}

KDtree::~KDtree()
{
}

/**
 * Finds the closest point within the tree,
 * wrt. the point given as first parameter.
 * @param _p point
 * @param maxdist2 maximal search distance.
 * @param threadNum Thread number, for parallelization
 * @return Pointer to the closest point
 */
double *KDtree::FindClosest(double *_p, double maxdist2, int threadNum) const
{
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  _FindClosest(Void(), threadNum);
  return params[threadNum].closest;
}

void KDtree::FindClosestKNNRange(double *_p, double maxdist2, vector<double*>& closest_list, size_t knn, int threadNum) const
{
  closest_list.clear();
  params[threadNum].closest_list.clear();

  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  _FindClosestKNNRange(Void(), threadNum);
  if (knn < 1) {
    /// range search
    closest_list = params[threadNum].closest_list;
  } else {
    /// knn range search
    // only sort if the list size is greater than knn
    if (params[threadNum].closest_list.size() > knn) {
        /// sort neighbor vector
        distance_comparison.origin = _p;
        sort(params[threadNum].closest_list.begin(), params[threadNum].closest_list.end(), distance_comparison);
    }
    for (size_t i = 0; i < params[threadNum].closest_list.size(); ++i)
      if (Dist2(_p, params[threadNum].closest_list[i]) <= maxdist2 && i < knn)
        closest_list.push_back(params[threadNum].closest_list[i]);
  }

}
