/*
 * scan_red implementation
 *
 * Copyright (C) Dorit Borrmann
 *
 * Released under the GPL version 3.
 *
 */


/**
 * @file
 * @brief Main program for reducing 3D scans.
 * 
 * Program to reduce scans for use with slam6d 
 * Usage: bin/scan_red -r <NR> 'dir',
 * Use -r for octree based reduction  (voxel size=<NR>)
 * and 'dir' the directory of a set of scans
 * Reduced scans will be written to 'dir/reduced'
 *
 * @author Dorit Borrmann. Automation Group, Jacobs University Bremen gGmbH, Germany. 
 */
#ifdef _MSC_VER
#if !defined _OPENMP && defined OPENMP 
#define _OPENMP
#endif
#endif

#define WANT_STREAM ///< define the WANT stream :)
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <errno.h>

#include "slam6d/scan.h"
#ifdef WITH_SCANSERVER
#include "scanserver/clientInterface.h"
#endif //WITH_SCANSERVER

#include "slam6d/globals.icc"

#ifdef _OPENMP
#include <omp.h>
#endif


#ifndef _MSC_VER
#include <getopt.h>
#else
#include "XGetopt.h"
#endif

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#include <windows.h>
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <strings.h>
#include <dlfcn.h>
#endif

#include <opencv2/opencv.hpp>

enum red_method {
OCTREE,
EQUIRECTANGULAR,
CYLINDRICAL,
MERCATOR,
CONIC
};

/**
 * Explains the usage of this program's command line parameters
 */
void usage(char* prog)
{
#ifndef _MSC_VER
  const string bold("\033[1m");
  const string normal("\033[m");
#else
  const string bold("");
  const string normal("");
#endif
  cout << endl
	  << bold << "USAGE " << normal << endl
	  << "   " << prog << " [options] -r <NR> directory" << endl << endl;
  cout << bold << "OPTIONS" << normal << endl

	  << bold << "  -e" << normal << " NR, " << bold << "--end=" << normal << "NR" << endl
	  << "         end after scan NR" << endl
	  << endl
	  << bold << "  -f" << normal << " F, " << bold << "--format=" << normal << "F" << endl
	  << "         using shared library F for input" << endl
	  << "         (chose F from {uos, uos_map, uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, riegl_txt, riegl_rgb, riegl_bin, zahn, ply})" << endl
	  << endl
	  << bold << "  -m" << normal << " NR, " << bold << "--max=" << normal << "NR" << endl
	  << "         neglegt all data points with a distance larger than NR 'units'" << endl
	  << endl
	  << bold << "  -M" << normal << " NR, " << bold << "--min=" << normal << "NR" << endl
	  << "         neglegt all data points with a distance smaller than NR 'units'" << endl
	  << endl
	  << bold << "  -r" << normal << " NR, " << bold << "--reduce=" << normal << "NR" << endl
	  << "         turns on octree based point reduction (voxel size=<NR>)" << endl
	  << endl
	  << bold << "  -s" << normal << " NR, " << bold << "--start=" << normal << "NR" << endl
	  << "         start at scan NR (i.e., neglects the first NR scans)" << endl
	  << "         [ATTENTION: counting naturally starts with 0]" << endl
	  << endl
         << bold << "  -S, --scanserver" << normal << endl
         << "         Use the scanserver as an input method and handling of scan data" << endl
	     << endl
         << bold << "  -R, --reductionmethod" << normal << endl
         << "         Method for reduction [OCTREE|EQUIRECTANGULAR|CYLINDRICAL|MERCATOR|CONIC] default: OCTREE" << endl
	     << endl
         << bold << "  -z, --resize" << normal << endl
         << "         factor by which to resize the range image" << endl
	     << endl
         << bold << "  -w, --width" << normal << endl
         << "         width of range image (default: 800)" << endl
	     << endl
         << bold << "  -h, --height" << normal << endl
         << "         height of range image (default: 600)" << endl
    	  << endl << endl;
  
  cout << bold << "EXAMPLES " << normal << endl
	  << "   " << prog << " -m 500 -r 5 dat" << endl
	  << "   " << prog << " --max=5000 -r 10.2 dat" << endl
	  << "   " << prog << " -s 2 -e 10 -r dat" << endl << endl;
  exit(1);
}

red_method string_to_reduction_method(string method)
{
    if(strcasecmp(method.c_str(), "OCTREE") == 0) return OCTREE;
    else if(strcasecmp(method.c_str(), "EQUIRECTANGULAR") == 0) return EQUIRECTANGULAR;
    else if(strcasecmp(method.c_str(), "CYLINDRICAL") == 0) return CYLINDRICAL;
    else if(strcasecmp(method.c_str(), "MERCATOR") == 0) return MERCATOR;
    else if(strcasecmp(method.c_str(), "CONIC") == 0) return CONIC;
    else throw std::runtime_error(std::string("reduction method ")+method+std::string(" is unknown"));
}

/** A function that parses the command-line arguments and sets the respective flags.
 * @param argc the number of arguments
 * @param argv the arguments
 * @param dir the directory
 * @param red using point reduction?
 * @param rand use randomized point reduction?
 * @param start starting at scan number 'start'
 * @param end stopping at scan number 'end'
 * @param maxDist - maximal distance of points being loaded
 * @param minDist - minimal distance of points being loaded
 * @param quiet switches on/off the quiet mode
 * @param veryQuiet switches on/off the 'very quiet' mode
 * @return 0, if the parsing was successful. 1 otherwise
 */
int parseArgs(int argc, char **argv, string &dir, double &red, 
		    int &start, int &end, int &maxDist, int &minDist, int &octree, 
		    IOType &type, bool &scanserver, red_method &rmethod, double &resize,
            unsigned int &width, unsigned int &height)
{
  bool reduced = false;
  bool resized = false;
  bool setwidth = false;
  bool setheight = false;
  int  c;
  // from unistd.h:
  extern char *optarg;
  extern int optind;

  /* options descriptor */
  // 0: no arguments, 1: required argument, 2: optional argument
  static struct option longopts[] = {
    { "format",          required_argument,   0,  'f' },  
    { "max",             required_argument,   0,  'm' },
    { "min",             required_argument,   0,  'M' },
    { "start",           required_argument,   0,  's' },
    { "end",             required_argument,   0,  'e' },
    { "reduce",          required_argument,   0,  'r' },
    { "octree",          optional_argument,   0,  'O' },
    { "scanserver",      no_argument,         0,  'S' },
    { "reductionmethod", required_argument,   0,  'R' },
    { "resize",          required_argument,   0,  'z' },
    { "width",           required_argument,   0,  'w' },
    { "height",          required_argument,   0,  'h' },
    { 0,           0,   0,   0}                    // needed, cf. getopt.h
  };

  cout << endl;
  while ((c = getopt_long(argc, argv, "f:r:s:e:m:M:O:SR:z:w:h:", longopts, NULL)) != -1)
    switch (c)
	 {
	 case 'r':
	   red = atof(optarg);
     reduced = true;
	   break;
	 case 's':
	   start = atoi(optarg);
	   if (start < 0) { cerr << "Error: Cannot start at a negative scan number.\n"; exit(1); }
	   break;
	 case 'e':
	   end = atoi(optarg);
	   if (end < 0)     { cerr << "Error: Cannot end at a negative scan number.\n"; exit(1); }
	   if (end < start) { cerr << "Error: <end> cannot be smaller than <start>.\n"; exit(1); }
	   break;
	 case 'f': 
    try {
      type = formatname_to_io_type(optarg);
    } catch (...) { // runtime_error
      cerr << "Format " << optarg << " unknown." << endl;
      abort();
    }
    break;
   case 'm':
	   maxDist = atoi(optarg);
	   break;
	 case 'O':
     if (optarg) {
       octree = atoi(optarg);
     } else {
       octree = 1;
     }
	   break;
	 case 'M':
	   minDist = atoi(optarg);
	   break;
    case 'S':
       scanserver = true;
       break;
    case 'R':
       rmethod = string_to_reduction_method(optarg);
       break;
    case 'z':
       resize = atof(optarg);
       resized = true;
       break;
    case 'w':
       width = atoi(optarg);
       setwidth = true;
       break;
    case 'h':
       height = atoi(optarg);
       setheight = true;
       break;
   case '?':
	   usage(argv[0]);
	   return 1;
      default:
	   abort ();
      }

  // some arguments only hold when octree is chosen
  switch (rmethod) {
      case OCTREE:
          if(!reduced) {
            cerr << "\n*** Reduction method missed ***" << endl;
            usage(argv[0]);
          }
          if(resized) {
              cerr << "\n*** No resizing allowed in octree mode ***" << endl;
              usage(argv[0]);
          }
          if(setwidth) {
              cerr << "\n*** You cannot set width in octree mode ***" << endl;
              usage(argv[0]);
          }
          if(setheight) {
              cerr << "\n*** You cannot set height in octree mode ***" << endl;
              usage(argv[0]);
          }
          break;
      case EQUIRECTANGULAR:
      case CYLINDRICAL:
      case MERCATOR:
      case CONIC:
          if(octree) {
              cerr << "\n*** -O can only be specified when choosing octree reduction ***" << endl;
              usage(argv[0]);
          }
          if(reduced) {
              cerr << "\n*** -r can only be specified when choosing octree reduction ***" << endl;
              usage(argv[0]);
          }
  }

  if (optind != argc-1) {
    cerr << "\n*** Directory missing ***" << endl;
    usage(argv[0]);
  }
  dir = argv[optind];

#ifndef _MSC_VER
  if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
  if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif

  return 0;
}

void octree_reduction(Scan *scan, double red, int octree, std::vector<cv::Vec3d> &result)
{
    scan->setReductionParameter(red, octree);
    scan->toGlobal();
    DataXYZ xyz_r(scan->get("xyz reduced"));
    for(unsigned int j = 0; j < xyz_r.size(); j++) {
        result.push_back(cv::Vec3d(xyz_r[j][0], xyz_r[j][1], xyz_r[j][2]));
    }
}

//Vertical angle of view of scanner
#define MAX_ANGLE 60.0
#define MIN_ANGLE -40.0

void equirectangular_reduction(Scan *scan, double resize, unsigned int width,
       unsigned int height, std::vector<cv::Vec3d> &result)
{
    cv::Mat iMap(width, height, CV_32FC(1), cv::Scalar::all(0));

    double kart[3], polar[3], phi, theta, range;
    unsigned int x, y;
    double xFactor = (double)width / 2 / M_PI;
    double yFactor = (double)height / ((MAX_ANGLE - MIN_ANGLE) / 180 * M_PI);
    //shift all the values to positive points on image
    cv::MatIterator_<cv::Vec4f> it, end;
    DataXYZ xyz(scan->get("xyz"));
    for(unsigned int i = 0; i < xyz.size(); ++i) {
        kart[0] = xyz[i][2]/100;
        kart[1] = xyz[i][0]/-100;
        kart[2] = xyz[i][1]/100;
        toPolar(kart, polar);
        //horizontal angle of view of [0:360] and vertical of [-40:60]
        theta = 0.5*M_PI - polar[0] - MIN_ANGLE/180.0 * M_PI;
        phi = 2.0*M_PI - polar[1];
        range = polar[2];
        x = (int) (xFactor * phi);
        if (x < 0) x = 0;
        if (x > width - 1) x = width - 1;
        y = (int) (yFactor * theta);
        if (y < 0) y = 0;
        if (y > height - 1) y = height - 1;
        // only map the nearest
        if (iMap.at<float>(x, y) > range || iMap.at<float>(x, y) == 0)
            iMap.at<float>(x, y) = range;
    }

    cv::Mat destMat;
    if (resize != 1.0) {
        destMat.create((int)(width*resize), (int)(height*resize), CV_32FC(1));
        cv::resize(iMap, destMat, destMat.size(), 0, 0);
    } else {
        destMat.create(width, height, CV_32FC(1));
        destMat = iMap;
    }

    for (int i = 0; i < destMat.size().width; ++i) {
        for (int j = 0; j < destMat.size().height; ++j) {
            range = destMat.at<float>(j, i);
            if (range == 0.0)
                continue;
            polar[2] = range;
            theta = i/yFactor/resize;
            phi = j/xFactor/resize;
            polar[0] = 0.5*M_PI - theta - MIN_ANGLE/180.0 * M_PI;
            polar[1] = 2.0*M_PI - phi;
            toCartesian(polar, kart);
            result.push_back(cv::Vec3d(kart[1]*(-100), kart[2]*100, kart[0]*100));
        }
    }
}

void cylindrical_reduction(Scan *scan, double resize, unsigned int width,
       unsigned int height, std::vector<cv::Vec3d> &result)
{
    cv::Mat iMap(width, height, CV_32FC(1), cv::Scalar::all(0));

    double kart[3], polar[3], phi, theta, range;
    unsigned int x, y;
    double xFactor = (double)width / 2 / M_PI;
    double yFactor = (double)height / (tan(MAX_ANGLE / 180.0 * M_PI) - tan(MIN_ANGLE / 180.0 * M_PI));
    double heightLow = tan((MIN_ANGLE)/180.0*M_PI);
    //shift all the values to positive points on image
    cv::MatIterator_<cv::Vec4f> it, end;
    DataXYZ xyz(scan->get("xyz"));
    for(unsigned int i = 0; i < xyz.size(); ++i) {
        kart[0] = xyz[i][2]/100;
        kart[1] = xyz[i][0]/-100;
        kart[2] = xyz[i][1]/100;
        toPolar(kart, polar);
        //horizontal angle of view of [0:360] and vertical of [-40:60]
        phi = 2.0* M_PI - polar[1];
        theta = tan(0.5*M_PI - polar[0]) - heightLow;
        range = polar[2];
        x = (int) (xFactor * phi);
        if (x < 0) x = 0;
        if (x > width - 1) x = width - 1;
        y = (int) (yFactor * theta);
        if (y < 0) y = 0;
        if (y > height - 1) y = height - 1;
        // only map the nearest
        if (iMap.at<float>(x, y) > range || iMap.at<float>(x, y) == 0)
            iMap.at<float>(x, y) = range;
    }

    cv::Mat destMat;
    if (resize != 1.0) {
        destMat.create((int)(width*resize), (int)(height*resize), CV_32FC(1));
        cv::resize(iMap, destMat, destMat.size(), 0, 0);
    } else {
        destMat.create(width, height, CV_32FC(1));
        destMat = iMap;
    }

    for (int i = 0; i < destMat.size().width; ++i) {
        for (int j = 0; j < destMat.size().height; ++j) {
            range = destMat.at<float>(j, i);
            if (range == 0.0)
                continue;
            polar[2] = range;
            theta = i/yFactor/resize;
            phi = j/xFactor/resize;
            polar[0] = 0.5*M_PI - atan(theta+heightLow);
            polar[1] = 2.0*M_PI - phi;
            toCartesian(polar, kart);
            result.push_back(cv::Vec3d(kart[1]*(-100), kart[2]*100, kart[0]*100));
        }
    }
}


/* FIXME: DOESNT WORK YET!! */
void mercator_reduction(Scan *scan, double resize, unsigned int width,
       unsigned int height, std::vector<cv::Vec3d> &result)
{
    cv::Mat iMap(width, height, CV_32FC(1), cv::Scalar::all(0));

    double kart[3], polar[3], phi, theta, range;
    unsigned int x, y;
    double xFactor = (double)width / 2 / M_PI;
    double yFactor = (double)height / (log(tan(MAX_ANGLE/180.0*M_PI) + (1/cos(MAX_ANGLE/180.0*M_PI)))
            - log(tan(MIN_ANGLE/180.0*M_PI) + (1/cos(MIN_ANGLE/180.0*M_PI))));
    double heightLow = log(tan(MIN_ANGLE/180.0*M_PI) + (1/cos(MIN_ANGLE/180.0*M_PI)));
    //shift all the values to positive points on image
    cv::MatIterator_<cv::Vec4f> it, end;
    DataXYZ xyz(scan->get("xyz"));
    for(unsigned int i = 0; i < xyz.size(); ++i) {
        kart[0] = xyz[i][2]/100;
        kart[1] = xyz[i][0]/-100;
        kart[2] = xyz[i][1]/100;
        toPolar(kart, polar);
        //horizontal angle of view of [0:360] and vertical of [-40:60]
        range = polar[2];
        phi = 2.0*M_PI - polar[1];
        theta = 0.5*M_PI - polar[0];
        theta = log(tan(theta) + (1/cos(theta))) - heightLow;
        x = (int) (xFactor * phi);
        if (x < 0) x = 0;
        if (x > width - 1) x = width - 1;
        y = height - 1 - (int) (yFactor * theta);
        if (y < 0) y = 0;
        if (y > height - 1) y = height - 1;
        // only map the nearest
        if (iMap.at<float>(x, y) > range || iMap.at<float>(x, y) == 0)
            iMap.at<float>(x, y) = range;
    }

    cv::Mat destMat;
    if (resize != 1.0) {
        destMat.create((int)(width*resize), (int)(height*resize), CV_32FC(1));
        cv::resize(iMap, destMat, destMat.size(), 0, 0);
    } else {
        destMat.create(width, height, CV_32FC(1));
        destMat = iMap;
    }

    for (int i = 0; i < destMat.size().width; ++i) {
        for (int j = 0; j < destMat.size().height; ++j) {
            range = destMat.at<float>(j, i);
            if (range == 0.0)
                continue;
            polar[2] = range;
            theta = (height*resize-1-i)/yFactor/resize;
            phi = j/xFactor/resize;
            polar[0] = 2.0/tan(pow(M_E, theta + heightLow)) - 0.5*M_PI;
            polar[1] = 2.0*M_PI - phi;
            toCartesian(polar, kart);
            result.push_back(cv::Vec3d(kart[1]*(-100), kart[2]*100, kart[0]*100));
        }
    }
}

/* FIXME: DOESNT WORK YET!! */
void conic_reduction(Scan *scan, double resize, unsigned int width,
       unsigned int height, std::vector<cv::Vec3d> &result)
{
    cv::Mat iMap(width, height, CV_32FC(1), cv::Scalar::all(0));

    double kart[3], polar[3], phi, theta, range, lambda, rho;
    unsigned int x, y;
    double phi_0 = 0.0;
    double lambda_0 = 0.0;
    double phi_1 = 0.0;
    double phi_2 = M_PI / 3.0; // 60 degrees
    double n = 0.5*(sin(phi_1) + sin(phi_2));
    double C = cos(phi_1)*cos(phi_1) + 2 * n * sin(phi_1);
    double rho_0 = sqrt(C - 2 * n * sin(phi_0))/n;
    double xshift = (sqrt(C+2*n)/n)*sin(n*(M_PI-lambda_0));
    double pwidth = 2.0*xshift;
    double yshift = rho_0-(sqrt(C-2*n)/n)*cos(n*(M_PI-lambda_0));
    double pheight = (sqrt(C+2*n)/n)*cos(-n*lambda_0)-(sqrt(C-2*n)/n)*cos(n*(M_PI-lambda_0));
    double xscale = width / pwidth;
    double yscale = height / pheight;
    double scale;
    if (xscale < yscale)
        scale = xscale;
    else
        scale = yscale;
    //shift all the values to positive points on image
    cv::MatIterator_<cv::Vec4f> it, end;
    DataXYZ xyz(scan->get("xyz"));
    for(unsigned int i = 0; i < xyz.size(); ++i) {
        kart[0] = xyz[i][2]/100;
        kart[1] = xyz[i][0]/-100;
        kart[2] = xyz[i][1]/100;
        toPolar(kart, polar);
        phi = polar[0] - M_PI/2.0;
        lambda = polar[1] - M_PI;
        range = polar[2];
        theta = n * (lambda - lambda_0);
        rho = sqrt(C - 2 * n * sin(phi))/n;
        x = (int) ((rho * sin(theta) + xshift)*scale);
        if (x < 0) x = 0;
        if (x > width - 1) x = width - 1;
        y = (int) ((rho * cos(theta) - rho_0 + yshift)*scale);
        if (y < 0) y = 0;
        if (y > height - 1) y = height - 1;
        // only map the nearest
        if (iMap.at<float>(x, y) > range || iMap.at<float>(x, y) == 0)
            iMap.at<float>(x, y) = range;
    }

    cv::Mat destMat;
    if (resize != 1.0) {
        destMat.create((int)(width*resize), (int)(height*resize), CV_32FC(1));
        cv::resize(iMap, destMat, destMat.size(), 0, 0);
    } else {
        destMat.create(width, height, CV_32FC(1));
        destMat = iMap;
    }

    for (int i = 0; i < destMat.size().width; ++i) {
        for (int j = 0; j < destMat.size().height; ++j) {
            range = destMat.at<float>(j, i);
            if (range == 0.0)
                continue;
            polar[2] = range;
            rho = sqrt(sqr(j/scale/resize+xshift) + sqr(i/scale/resize - yshift + rho_0));
            theta = 1/((j/scale/resize+xshift)/(i/scale/resize - yshift + rho_0));
            phi = 1/sin((C-rho*rho*n*n)/(2*n));
            lambda = lambda_0 + theta/n;
            polar[0] = phi + M_PI/2.0;
            polar[1] = lambda + M_PI;
            toCartesian(polar, kart);
            result.push_back(cv::Vec3d(kart[1]*(-100), kart[2]*100, kart[0]*100));
        }
    }
}


/**
 * Main program for reducing scans.
 * Usage: bin/scan_red -r <NR> 'dir',
 * Use -r for octree based reduction  (voxel size=<NR>)
 * and 'dir' the directory of a set of scans
 * Reduced scans will be written to 'dir/reduced'
 * 
 */
int main(int argc, char **argv)
{

  cout << "(c) Jacobs University Bremen, gGmbH, 2010" << endl << endl;
  
  if (argc <= 1) {
    usage(argv[0]);
  }

  // parsing the command line parameters
  // init, default values if not specified
  string dir;
  double red   = -1.0;
  int    start = 0,   end = -1;
  int    maxDist    = -1;
  int    minDist    = -1;
  int    octree     = 0;
  bool   scanserver = false;
  double resize     = 1.0;
  unsigned int width = 800;
  unsigned int height = 600;
  red_method rmethod = OCTREE;
  IOType type    = RIEGL_TXT;
  
  parseArgs(argc, argv, dir, red, start, end, maxDist, minDist, octree, type,
          scanserver, rmethod, resize, width, height);

  Scan::openDirectory(scanserver, dir, type, start, end);

  if(Scan::allScans.size() == 0) {
    cerr << "No scans found. Did you use the correct format?" << endl;
    exit(-1);
  }

  string reddir = dir + "reduced"; 
 
#ifdef _MSC_VER
  int success = mkdir(reddir.c_str());
#else
  int success = mkdir(reddir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
#endif
  if(success == 0) { 
    cout << "Writing scans to " << reddir << endl;
  } else if(errno == EEXIST) {
    cout << "Directory " << reddir << " exists already.  CONTINUE" << endl; 
  } else {
    cerr << "Creating directory " << reddir << " failed" << endl;
    exit(1);
  }

  for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
    Scan* scan = *it;
    scan->setRangeFilter(maxDist, minDist);
    const double* rPos = scan->get_rPos();
    const double* rPosTheta = scan->get_rPosTheta();

    const char* id = scan->getIdentifier();
    cout << "Writing Scan No. " << id;
    cout << " with " << scan->size<DataXYZ>("xyz") << " points" << endl; 
    string scanFileName;
    string poseFileName;

    scanFileName = dir + "reduced/scan" + id + ".3d";
    poseFileName = dir  + "reduced/scan" + id + ".pose";

    ofstream redptsout(scanFileName.c_str());

    std::vector<cv::Vec3d> xyz_r;

    switch (rmethod) {
        case OCTREE:
            octree_reduction(scan, red, octree, xyz_r);
            break;
        case EQUIRECTANGULAR:
            equirectangular_reduction(scan, resize, width, height, xyz_r);
            break;
        case CYLINDRICAL:
            cylindrical_reduction(scan, resize, width, height, xyz_r);
            break;
        case MERCATOR:
            /*mercator_reduction(scan, resize, width, height, xyz_r);
            break;*/
        case CONIC:
            /*conic_reduction(scan, resize, width, height, xyz_r);
            break;*/
        default:
            throw std::runtime_error(std::string("not implemented"));
    }

    for(std::vector<cv::Vec3d>::iterator it = xyz_r.begin(); it < xyz_r.end(); ++it) {
        redptsout << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
    }

    redptsout.close();
    redptsout.clear();
    
    ofstream posout(poseFileName.c_str());
    
    posout << rPos[0] << " " 
           << rPos[1] << " " 
           << rPos[2] << endl   
           << deg(rPosTheta[0]) << " " 
           << deg(rPosTheta[1]) << " " 
           << deg(rPosTheta[2]) << endl;  
    
    posout.close();
    posout.clear();
  }

  cout << endl << endl;
  cout << "Normal program end." << endl << endl;
}

