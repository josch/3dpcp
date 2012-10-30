#ifndef __NORMALS_H
#define __NORMALS_H

#include <iostream>
#include <string>
#include <fstream>
#include <errno.h>

#include <slam6d/io_types.h>
#include <slam6d/globals.icc>
#include <slam6d/scan.h>
#include "slam6d/fbr/panorama.h"
#include <scanserver/clientInterface.h>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

enum normal_method {AKNN, ADAPTIVE_AKNN, PANORAMA, PANORAMA_FAST};

void calculateNormalsAKNN(vector<Point>&,vector<Point>&, int, const double[3]);
void calculateNormalsAdaptiveAKNN(vector<Point>&,vector<Point>&,int, int, const double[3]);
void calculateNormalsPANORAMA(vector<Point>&,vector<Point>&,vector< vector< vector< cv::Vec3f > > >,const double[3]);

#endif
