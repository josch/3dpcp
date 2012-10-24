#ifndef SRI_H
#define SRI_H

#include "normals/point.h"
#include <vector>

class SRI {
 public: 
  SRI(int _method = 0, int _factor = 10);
  ~SRI();

  void addPoint(double _x, double _y, double _z);
  
  int input(const char* name);
  void writePGM(const char* name);
  
  void calculateNormals();
  
   
  int factor, method;
  int width, height, max;
  double** image;
  std::vector<PointN*> points;
  
  void sobel(int i, int j, double& dRdTheta, double& dRdPhi);
  void smooth();

  void getNormalSRI(PointN *p, double* n);
  void getNormalLS(PointN *p, double* n);


};




#endif // SRI_H
