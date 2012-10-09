#include <vector>

#include "image.h"
#include "misc.h"

class SRI {
public:
	class point
	{
	public:
		double theta; // azimuth
		double phi; // elevation
		double r; // range
		int x, y; // image coordinates
		double reflectance;

//		double orig_x, orig_y, orig_z;
		rgb color;
		point(double _x, double _y, double _z, double reflectance=0);
	};

	SRI();

	SRI(int width, int height, int point_size);

	void addPoint(double x, double y, double z, double reflectance, rgb color);

	image<rgb>* getImage(bool useColor = true, bool useReflectance = false);

	void reserve(int qty);

private:
	std::vector<point> points;
	image<int>* references;
	double minR, maxR, avgR;
	int minX, maxX, minY, maxY;
	double minRef, maxRef;
	int width, height, ps;
};
