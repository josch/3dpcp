#ifndef POINT_H
#define POINT_H

#include <cmath>

enum construction {
  Cartesian, 
  Spherical, 
  Image
};




class Normal {
 public:
  Normal(double _x, double _y, double _z);
  ~Normal();
  int* toColor(int max = 1000);

 private:
  void normalize();
  double x, y, z;
  int* color;

};




class PointN {
 public:
  PointN(construction c, double _a, double _b, double _c, int factor); 
   
  void getCartesian(double &_x, double &_y, double &_z);
  void getSpherical(double &_rho, double &_phi, double &_theta);
  void getImage(int &_x, int &_y, double &_r);
  void setNormal(Normal _n) { n = _n; }
  Normal getNormal() const { return n; }

  static double sqr(double x) { return x*x; }
    

 private:
  double x, y, z;
  double rho, phi, theta;
  int ix, iy;
  Normal n;

  
};




#endif // POINT_H
