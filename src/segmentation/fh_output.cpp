/*
 * fh_output.cpp
 *
 *  Created on: May 8, 2012
 *      Author: msima
 */

#include <unistd.h>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#include "slam6d/scan_io_xyzr.h"
#include "slam6d/scan_io_rxp.h"
#include "slam6d/scan_io_velodyne.h"
#include "slam6d/scan_io.h"
#include "slam6d/scan.h"

#include "segmentation/segment-image.h"
#include "segmentation/SRI.h"
#include "segmentation/pnmfile.h"

template<typename T> T parse(string X)
{
	stringstream ss(X);
	T result;
	ss >> result;
	return result;
}

int main(int argc, char** argv) {
	reader_type type;
	int start = 0;
	int end = -1;
	int height = -1;
	int width = -1;
	string output = "output.ppm";
	bool with_reflectance = false;
	int point_size = 1;
	int reserve = -1;
	int reserve2 = -1;
    bool oneColor = false;


    const char* options = "f:s:e:h:w:o:rp:R:T:1";
	char c;
	while ( (c=getopt(argc, argv, options)) != -1 )
	{
		switch ( c )
		{
			case 'f':
				if ( ! Scan::toType(optarg, type) )
				{
					cerr << "Invalid type " << optarg << endl;
				}
				break;
			case 's':
				start = parse<int>(optarg);
				break;
			case 'e':
				end = parse<int>(optarg);
				break;
			case 'h':
				height = parse<int>(optarg);
				break;
			case 'w':
				width = parse<int>(optarg);
				break;
			case 'o':
				output = optarg;
				break;
			case 'r':
				with_reflectance = true;
				break;
			case 'p':
				point_size = parse<int>(optarg);
				break;
			case 'R':
				reserve = parse<int>(optarg);
				break;
			case 'T':
				reserve2 = parse<int>(optarg);
				break;
            case '1':
                oneColor = true;
                break;
		}
	}

	SRI sri(width, height, point_size);
	if ( reserve>0 )
		sri.reserve(reserve);

	for (; optind < argc; ++optind )
	{
		string file = argv[optind];
		cerr << file << endl;
		if ( file[file.length()-1]!= '/' )
			file += '/';

		Scan s;
		s.readScans(type, start, end, file, -1, -1, false, reserve2);
		cerr << s.allScans.size() << endl;
		for (int i=s.allScans.size()-1; i>=0; --i) {
			const vector<Point>* v = s.allScans[i]->get_points();
			rgb c;
			if ( oneColor )
				c.r = 250, c.g = c.b = 255;
			else
				c = random_rgb();
			for (size_t j=0; j<v->size(); ++j)
			{
				Point p = (*v)[j];
				sri.addPoint(p.x, p.y, p.z, p.reflectance, c);
			}
			s.allScans[i]->clearPoints();
		}
	}
	savePPM(sri.getImage(true, with_reflectance), output.c_str());

	cout << "DONE!" << endl;
	return 0;
}


