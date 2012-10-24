#include "normals/point.h"
#include <cassert>
#include <cstdlib>


static const double rad2deg = 180.0 / 3.1415;
static const double deg2rad = 3.1415 / 180.0;

// constructor
PointN::PointN(construction c, double _x, double _y, double _z, int factor) :
  n(Normal(0,0,0))
{
  switch (c) {

  case Cartesian:
    x = _x; 
    y = _y; 
    z = _z;
    
    rho = sqrt(sqr(x) + sqr(y) + sqr(z));
    phi = acos(z/rho)- 30.0 * deg2rad;
    theta = atan2(y, x) + 180.0 * deg2rad;
    
    ix =  (int)floor(theta*rad2deg*factor + 0.5);
    iy =  (int)floor(phi*rad2deg*factor + 0.5);
    
    break;

  case Spherical:
    rho = _x; theta = _y; phi = _z;
    x = rho * sin((phi+30)*deg2rad) * cos(theta*deg2rad);
    y = rho * sin((phi+30)*deg2rad) * sin(theta*deg2rad);
    z = rho * cos((phi+30)*deg2rad);
    ix =  (int)floor(theta*factor + 0.5);
    iy =  (int)floor(phi*factor + 0.5);
    break;

  case Image:
    ix = (int)_x; 
    iy = (int)_y; 

    rho = _z+1;    
    phi = ((double)(iy) / factor) * deg2rad;
    theta = ((double)(ix) / factor) * deg2rad;
    
    x = rho * sin(phi + 30.0*deg2rad) * cos(theta - 180*deg2rad);
    y = rho * sin(phi + 30.0*deg2rad) * sin(theta - 180*deg2rad);
    z = rho * cos(phi + 30.0*deg2rad);
    
    break;
  }

}


void PointN::getCartesian(double &_x, double &_y, double &_z) {
  _x = x; _y = y; _z = z;
}


void PointN::getSpherical(double &_rho, double &_phi, double &_theta) {
  _rho = rho; _phi = phi; _theta = theta;
}


void PointN::getImage(int &_x, int &_y, double &_r ) {
  _x = ix; _y = iy; _r = rho;
}







Normal::Normal(double _x, double _y, double _z) : x(_x), y(_y), z(_z) 
{  
  normalize();
  color = NULL;
}


Normal::~Normal() {
  delete [] color;
}


void Normal::normalize() {
  double m = sqrt(x*x + y*y + z*z);
  x /= m; y /= m; z /= m;
  assert(abs(x) <= 1.0);
  assert(abs(y) <= 1.0);
  assert(abs(z) <= 1.0);
}


  
int* Normal::toColor(int max) {
  
  if (color == NULL) {
    color = new int[3];
    
    color[0] = floor(x*max/2.0 + 0.5) + max/2.0;
    color[1] = floor(y*max/2.0 + 0.5) + max/2.0;
    color[2] = floor(z*max/2.0 + 0.5) + max/2.0;
  }

  return color;
}
