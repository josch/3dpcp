#include <newmat/newmatap.h>
#include "normals/SRI.h"
#include <slam6d/globals.icc>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace NEWMAT;

static const double deg2rad = 3.1415 / 180.0;


SRI::SRI(int _method, int _factor) {
  method = _method;
  factor = _factor;
  width = height = 0;
  image = NULL;
} 


SRI::~SRI() {
  if (image != NULL) {
    for (int i = 0; i < height; i++)
      delete [] image[i];
    delete [] image;
  }
}



// add a point with its cartesian coordinates x, y, z
void SRI::addPoint(double x, double y, double z) {

  PointN *p = new PointN(Cartesian, x, z, y, factor);
    
  if (image == NULL) {
    height = (int)factor*100;
    width = (int)factor*360;
    image = new double*[height];
    for (int i = 0; i < height; i++)
      image[i] = new double[width];
    for (int i = 0; i < height; i++)
      for (int j = 0; j < width; j++)
	image[i][j] = 0;
    max = 0;
  }


  int ix, iy;
  double r;
  p->getImage(ix, iy, r);

 
  // save the point in the image array if applicable
  if (iy >= 0 && ix >= 0 && iy < height && ix < width && 
      (image[iy][ix] > r || image[iy][ix] == 0) /*&& r < 1000*/) {
    image[iy][ix] = r;
    points.push_back(p);
    if (r > max)
      max = r;
  }
  else
    delete p;
  
}


void SRI::calculateNormals() {
  ofstream show_out("normals.nrm", ios::out);
  ofstream out;

  cerr << "Calculating normals using ";
  void (SRI::*calcMethod)(PointN*, double*) = NULL;
  if (method == 0) {
    cerr << "SRI method... ";
    calcMethod = &SRI::getNormalSRI;
    out.open("normals_sri.pgm", ios::out);
  }
  else {
    cerr << "LS method... ";
    calcMethod = &SRI::getNormalLS;
    out.open("normals_ls.pgm", ios::out);
  }
  

  vector<double*> calculated_normals;
  
  unsigned long ms = GetCurrentTimeInMilliSec();
  
  for (unsigned int i = 0 ; i < points.size(); i++) { 
    double* n = new double[3];
    (this->*calcMethod)(points[i], n);   
    calculated_normals.push_back(n);
  }

  ms = GetCurrentTimeInMilliSec() - ms;

  cerr << "done in " << ms << " ms. Total normals: " << calculated_normals.size() << endl;
  
  out << "P3\n" << width << " " << height-20 << endl << 1000 << endl;

  for (unsigned int i = 0 ; i < calculated_normals.size(); i++) {
    double x,y,z;
    points[i]->getCartesian(x,z,y);
    show_out << x << " " << y << " " << z  << " " << 
      calculated_normals[i][0] << " " << calculated_normals[i][1] << " " << 
      calculated_normals[i][2] << endl;
    

    if (i > width*(height-20)) continue;

    int rgbN[3];
    for (int j = 0; j < 3; j++) 
      //rgbN[j] = floor(abs(calculated_normals[i][j])*1000 + 0.5);
      rgbN[j] = floor(calculated_normals[i][j]*500.0 + 0.5)+500; 

    out << rgbN[0] << " " << rgbN[1] << " " << rgbN[2] << " ";
  }


  show_out.close();  
  out.close();  

  for (unsigned i = 0; i < calculated_normals.size(); i++)
    delete [] calculated_normals[i];
  calculated_normals.clear();

}




// read a grayscale SRI and get the points from there
int SRI::input(const char* name) {
 
  ifstream in(name, ios::in);
  string buffer;
  
  if (image != NULL) {
    for (int i = 0; i < height; i++)
      delete [] image[i];
    delete [] image;
  }

  for (unsigned int i = 0; i < points.size(); i++)
    delete points[i];
  points.clear();

  in >> buffer;
  in >> width >> height >> max;
  
  image = new double*[height];
  for (int i = 0; i < height; i++)
    image[i] = new double[width];

  for (int i = 0; i < height; i++)
    for (int j = 0; j < width; j++) 
      in >> image[i][j];

  in.close();

  // fix any "holes" 
  smooth();

  // save the read points
  for (int i = 0; i < height; i++) 
    for (int j = 0; j < width; j++) 
      points.push_back(new PointN(Image, j, i, image[i][j], factor));

  return 0;
}




// sobel operator estimating the parital derivatives at a point
void SRI::sobel(int i, int j, double& dRdTheta, double& dRdPhi) {

  dRdTheta = dRdPhi = 0.0;

  if (i == 0 || i == height-1 || j == 0 || j == width-1)
    return; 

  dRdPhi += 10*image[i-1][j];
  dRdPhi += 3 *image[i-1][j-1];
  dRdPhi += 3 *image[i-1][j+1];
  dRdPhi -= 10*image[i+1][j];
  dRdPhi -= 3 *image[i+1][j-1];
  dRdPhi -= 3 *image[i+1][j+1];
    
  dRdTheta += 10*image[i][j-1];
  dRdTheta += 3 *image[i-1][j-1];
  dRdTheta += 3 *image[i+1][j-1];
  dRdTheta -= 10*image[i][j+1];
  dRdTheta -= 3 *image[i-1][j+1];
  dRdTheta -= 3 *image[i+1][j+1];
}  



// LS method 
void SRI::getNormalLS(PointN *p, double* n) {
  
  int i, j;
  double rho;
  double phi, theta;
  p->getImage(j,i,rho);
  double _r;
  p->getSpherical(_r, phi, theta);

  if (_r < 10) {
    for (int i = 0; i < 3; i++) 
      n[i] = 0;
    return;
  }


  phi += (30.0*deg2rad);

  Matrix b(3, 1); b = 0;
  Matrix M(3, 3); M = 0; 


  if (i>0) {
    Matrix v(3,1);
    double _th, _ph; double _rh;
    _rh = image[i-1][j]+1;
    _ph = phi - (1.0 / (double)factor) * deg2rad;
    _th = theta;
    v << cos(_th)*sin(_ph) << sin(_th)*sin(_ph) << cos(_ph);
    M += v * v.t();
    v /= _rh;
    b += v;
  }
  

  if (j>0) {
    Matrix v(3,1);
    double _th, _ph; double _rh;
    _rh = image[i][j-1]+1;
    _ph = phi;
    _th = theta - (1.0 / (double)factor) * deg2rad;
    v << cos(_th)*sin(_ph) << sin(_th)*sin(_ph) << cos(_ph);
    M += v * v.t();
    v /= _rh;
    b += v;
  }


  if (i<height-1) {
    Matrix v(3,1);
    double _th, _ph; double _rh;
    _rh = image[i+1][j]+1;
    _ph = phi + (1.0 / (double)factor) * deg2rad;
    _th = theta;
    v << cos(_th)*sin(_ph) << sin(_th)*sin(_ph) << cos(_ph);
    M += v * v.t();
    v /= _rh;
    b += v;
  }


  if (j<width-1) {
    Matrix v(3,1);
    double _th, _ph; double _rh;
    _rh = image[i][j+1]+1;
    _ph = phi;
    _th = theta + (1.0 / (double)factor) * deg2rad;
    v << cos(_th)*sin(_ph) << sin(_th)*sin(_ph) << cos(_ph);
    M += v * v.t();
    v /= _rh;
    b += v;
  }
    
  Matrix N = M.i() * b;
  
  n[0] = N(1,1);
  n[1] = N(2,1);
  n[2] = -N(3,1);

  double m = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] /= m; n[1] /= m; n[2] /= m;

}



// SRI method
void SRI::getNormalSRI(PointN *p, double* n) {

  double theta, phi, rho;
  int i, j;
  double r;
  p->getSpherical(rho, phi, theta);
  p->getImage(j, i, r); 

  // if no point return 0 normal
  if (r < 10) {
    for (int i = 0; i < 3; i++) 
      n[i] = 0;
    return;
  }

  // partial derivative values
  double dRdTheta, dRdPhi;

  sobel(i, j, dRdTheta, dRdPhi);

  // add more weight to derivatives
  dRdTheta *= factor*2;
  dRdPhi *= factor*2;

  phi += (30.0*deg2rad);  
  
  n[0] = cos(theta) * sin(phi) - sin(theta) * dRdTheta / rho / sin(phi) + 
    cos(theta) * cos(phi) * dRdPhi / rho;

  n[1] = sin(theta) * sin(phi) + cos(theta) * dRdTheta / rho / sin(phi) + 
    sin(theta) * cos(phi) * dRdPhi / rho;
  
  n[2] =  cos(phi) - sin(phi) * dRdPhi / rho;

  n[2] = -n[2];

  double m = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
  n[0] /= m; n[1] /= m; n[2] /= m;
  
}




void SRI::writePGM(const char* name) {

  ofstream pgm (name, ios::out);
  ofstream human_pgm("human_sri.pgm", ios::out);

  smooth();

  pgm << "P2\n" << width << " " << height << "\n" << max << "\n";
  human_pgm << "P2\n" << width << " " << height << "\n" << max << "\n";
  
  for (int i = 0; i < height; i++) {
    for (int j = width-1; j >= 0; j--) {
      human_pgm << (int)(image[i][j]) << " ";
      pgm << (image[i][j]) << " ";
    }
    human_pgm << endl;
    pgm << endl;
  }
  pgm.close();
  human_pgm.close();
}




void SRI::smooth() {

  double ** copy = new double*[height];
  for (int i = 0; i < height; i++)
    copy[i] = new double[width];
  for (int i = 1; i < height-1; i++)
    for (int j = 1; j < width-1; j++) 
      copy[i][j] = image[i][j];

  for (int i = 2; i < height-2; i++)
    for (int j = 2; j < width-2; j++) 
      if (image[i][j] == 0 || image[i][j] > max)
	copy[i][j] = (image[i+1][j] + image[i-1][j] + image[i+1][j-1] + image[i-1][j-1] 
		      + image[i+1][j+1] + image[i-1][j+1] + image[i][j-1] + image[i][j+1]) / 8;

  for (int i = 1; i < height-1; i++)
    for (int j = 1; j < width-1; j++) 
      image[i][j] = copy[i][j];
  
    for (int i = 0; i < height; i++)
      delete [] copy[i];
    delete [] copy;

}
