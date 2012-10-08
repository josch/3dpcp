#include "segment/SRI.h"
#include <iostream>
#include <complex>

using namespace std;

inline double sqr(double x) { return x*x; }

SRI::point::point(double _x, double _y, double _z, double _r) {
	  static const double rad2deg = 180.0 / 3.1415;
	  
	  r     = sqrt(sqr(_x) + sqr(_y) + sqr(_z));
	  theta = acos(_y/ r) * rad2deg;
	  phi   = atan2(_z, _x) * rad2deg;

	  x =  floor(phi  *10 + 0.5) + 1800;
	  y =  floor(theta*10 + 0.5);

//	  orig_x = _x; orig_y = _y; orig_z = _z;
	  reflectance = _r;
}

SRI::SRI() {
	minR = 10000; maxR = 0; avgR = 0;
	minX = 10000; maxX = 0;
	minY = 10000; maxY = 0;
	minRef = 10000; maxRef = -10000;
	width = -1; height = -1;
}

SRI::SRI(int width, int height, int point_size)
{
	minR = 10000; maxR = 0; avgR = 0;
	minX = 10000; maxX = 0;
	minY = 10000; maxY = 0;
	minRef = 10000; maxRef = -10000;
	this->width = width;
	this->height = height;
	ps = point_size;
}

void SRI::addPoint(double x, double y, double z, double reflectance, rgb color)
{
	point p(x,y,z, reflectance);
	p.color = color;
	points.push_back(p);
	minR = min(minR, p.r); maxR = max(maxR, p.r);
	minX = min(minX, p.x); maxX = max(maxX, p.x);
	minY = min(minY, p.y); maxY = max(maxY, p.y);
	minRef = min(minRef, p.reflectance); maxRef = max(maxRef, p.reflectance);
}

inline rgb modulate(rgb c, float alpha)
{

	float r = c.r / 255.0;
	float g = c.g / 255.0;
	float b = c.b / 255.0;

	float h, s, v;

	// convert RGB to HSV
	float _min = min(r, min(g,b));
	float _max = max(r, max(g,b));
	float _d = _max - _min;
	s = _max == 0 ? 0 : _d/_max;

	if ( s==0 )
		return c;

	if ( _max == _min )
		h = 0;
	else
	{
		if ( _max==r )
			h = (g-b) / _d + (g<b ? 6 : 0);
		else if (_max == g)
			h = (b-r) / _d + 2.0;
		else
			h = (r-g) / _d + 4.0;
	}
	v = _max;

	// adjust brightness
	v = alpha;


	// convert HSV back to RGB
	int i = h;
	float f = h - i;
	float p = v * (1-s);
	float q = v * (1-f*s);
	float t = v * (1-(1-f)*s);

	switch ( i%6 ) {
    case 0: r = v, g = t, b = p; break;
    case 1: r = q, g = v, b = p; break;
    case 2: r = p, g = v, b = t; break;
    case 3: r = p, g = q, b = v; break;
    case 4: r = t, g = p, b = v; break;
    case 5: r = v, g = p, b = q; break;
	}

	rgb ret;
	ret.r = r*255;
	ret.g = g*255;
	ret.b = b*255;

	return ret;
}

image< rgb >* SRI::getImage(bool useColor, bool useReflectance)
{
	if ( width == -1 )
		width = maxX-minX+1;
	double coeffX = 1.0 * width / (maxX-minX+1);

	if ( height == -1 )
		height = maxY - minY + 1;
	double coeffY = 1.0 * height / (maxY-minY+1);

	image<rgb>* rez = new image<rgb>(width, height, true);
	for (size_t i=0; i<points.size(); ++i)
	{
		uchar r = (points[i].r - minR) / (maxR - minR) * 250;
		int X = (points[i].x - minX) * coeffX;
		int Y = (points[i].y - minY) * coeffY;
		float ref = (points[i].reflectance-minRef) / (maxRef - minRef);

		for (int x = max(0, X-ps+1); x<=min(width-1, X+ps-1); ++x)
			for (int y=max(0, Y-ps+1); y<=min(height-1, Y+ps-1); ++y)
			{
				if ( useColor )
				{
					if ( !useReflectance )
						rez->access[y][x] = points[i].color;
					else
					{
						rgb c = modulate(points[i].color, ref);
						rez->access[y][x] = c;
					}
				}
				else
				{
					if ( !useReflectance )
						rez->access[y][x].r = rez->access[y][x].g = rez->access[y][x].b = r;
					else
						rez->access[y][x].r = rez->access[y][x].g = r, rez->access[y][x].b = ref;
				}
			}
	}
	return rez;
}

void SRI::reserve(int qty)
{
	points.reserve(qty);
}
