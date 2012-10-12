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

#include "slam6d/globals.icc"
#include "slam6d/io_utils.h"
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
	IOType type = RXP;
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
  bool scanserver = false;
  string dir;

  WriteOnce<IOType> w_type(type);
  WriteOnce<int> w_start(start), w_end(end);

	char c;
	while ( (c=getopt(argc, argv, "f:s:e:h:w:o:rp:R:T:1")) != -1 )
		switch ( c )
		{
			 case 'f': 
        try {
          w_type = formatname_to_io_type(optarg);
        } catch (...) { // runtime_error
          cerr << "Format " << optarg << " unknown." << endl;
          abort();
        }
        break;
			case 's':
				w_start = parse<int>(optarg);
				break;
			case 'e':
				w_end = parse<int>(optarg);
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
      default:
	     abort ();
		}

  if (optind != argc-1) {
    cerr << "\n*** Directory missing ***" << endl;
    exit(1);
  }
  dir = argv[optind];

#ifndef _MSC_VER
  if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
  if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif
  if(start > end) {
    cerr << "\n*** scan start number has to be smaller than scan end number ***" << endl;
    exit(1);
  } 

  parseFormatFile(dir, w_type, w_start, w_end);

	SRI sri(width, height, point_size);
	if ( reserve>0 )
		sri.reserve(reserve);

  Scan::openDirectory(scanserver, dir, type, start, end);

  if (Scan::allScans.size() == 0) {
    cerr << "No scans found. Did you use the correct format?" << endl;
    exit(-1);
  }

  for (std::vector<Scan*>::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
    Scan* scan = *it;
    DataXYZ xyz = scan->get("xyz");
    DataReflectance xyz_reflectance = scan->get("reflectance");
    unsigned int nPoints = xyz.size();
    
    rgb c;
    if ( oneColor) 
      c.r = 250, c.g = c.b = 255;
    else
      c = random_rgb();

    for(unsigned int i = 0; i < nPoints; i++) {
      float x, y, z, reflectance;
      x = xyz[i][0];
      y = xyz[i][1];
      z = xyz[i][2];
      reflectance = xyz_reflectance[i];

      //normalize the reflectance
      reflectance += 32;
      reflectance /= 64;
      reflectance -= 0.2;
      reflectance /= 0.3;
      if (reflectance < 0) reflectance = 0;
      if (reflectance > 1) reflectance = 1;

      sri.addPoint(x, y, z, reflectance, c);
    }
  }

	savePPM(sri.getImage(true, with_reflectance), output.c_str());

  Scan::closeDirectory();

	cout << "DONE!" << endl;
	return 0;
}


