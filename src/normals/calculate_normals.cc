/*
 * calculateNormals implementation
 *
 * Copyright (C) Johannes Schauer, Razvan Mihaly
 *
 * Released under the GPL version 3
 *
 */

#include "ANN/ANN.h"

#include "slam6d/point.h"
#include "slam6d/scan.h"
#include "slam6d/globals.icc"

#include <string>
using std::string;

#include <iostream>
using std::cout;
using std::endl;
using std::vector;

#include <algorithm>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

/*
 * validates input type specification
 */
void validate(boost::any& v, const std::vector<std::string>& values,
    IOType*, int) {
    if (values.size() == 0)
        throw std::runtime_error("Invalid model specification");
    string arg = values.at(0);
    try {
        v = formatname_to_io_type(arg.c_str());
    } catch (...) { // runtime_error
        throw std::runtime_error("Format " + arg + " unknown.");
    }
}

/*
 * parse commandline options, fill arguments
 */
void parse_options(int argc, char **argv, int &start, int &end,
        bool &scanserver, string &dir, IOType &iotype,
        int &maxDist, int &minDist, int &normalMethod, int &knn)
{
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "output this help message");

    po::options_description input("Input options");
    input.add_options()
        ("start,s", po::value<int>(&start)->default_value(0),
         "start at scan <arg> (i.e., neglects the first <arg> scans) "
         "[ATTENTION: counting naturally starts with 0]")
        ("end,e", po::value<int>(&end)->default_value(-1),
         "end after scan <arg>")
        ("format,f", po::value<IOType>(&iotype)->default_value(UOS),
         "using shared library <arg> for input. (chose F from {uos, uos_map, "
         "uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, "
         "riegl_txt, riegl_rgb, riegl_bin, zahn, ply})")
        ("max,M", po::value<int>(&maxDist)->default_value(-1),
         "neglegt all data points with a distance larger than <arg> 'units")
        ("min,m", po::value<int>(&minDist)->default_value(-1),
         "neglegt all data points with a distance smaller than <arg> 'units")
        ("scanserver,S", po::bool_switch(&scanserver),
         "Use the scanserver as an input method and handling of scan data")
        ("normalMethod,N", po::value<int>(&normalMethod)->default_value(0),
         "choose the method for computing normals from 0 to 4")
        ("knn,K", po::value<int>(&knn)->default_value(1),
         "select the k in kNN search")
;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(hidden);

    // options visible with --help
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input);

    // positional argument
    po::positional_options_description pd;
    pd.add("input-dir", 1);

    // process options
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(all).positional(pd).run(), vm);
    po::notify(vm);

    // display help
    if (vm.count("help")) {
        cout << cmdline_options;
        exit(0);
    }

    // add trailing slash to directory if not present yet
    if (dir[dir.length()-1] != '/') dir = dir + "/";
}

/*
void writeScan(const string& dir, 
        unsigned int scan_number,
        const vector<Point>& points) {

    stringstream ss; 
    ss << dir << "scan" << std::setw(3) << std::setfill('0') << scan_number << ".3d"; 
    ofstream scan_file;
    scan_file.open(ss.str().c_str());
    for(size_t i = 0;  i < points.size(); ++i) {
        scan_file << points[i].x << " " << points[i].y << " " << points[i].z << "\n";  
    }
    scan_file.close();

    ss.clear(); ss.str(string());
    ss << dir << "scan" << std::setw(3) << std::setfill('0') << scan_number << ".pose";
    ofstream pose_file; 
    pose_file.open(ss.str().c_str());
    pose_file << 0 << " " << 0 << " " << 0 << "\n" << 0 << " " << 0 << " " << 0 << "\n";
    pose_file.close();
}
*/

void computeNeighbors(const vector<Point>& points, int knn, double eps) 
{
	ANNpointArray point_array = annAllocPts(points.size(), 3);
	for (size_t i = 0; i < points.size(); ++i) {
		point_array[i] = new ANNcoord[3];
		point_array[i][0] = points[i].x;
		point_array[i][1] = points[i].y;
		point_array[i][2] = points[i].z;
	}

	ANNkd_tree t(point_array, points.size(), 3);
	ANNidxArray n = new ANNidx[knn];
	ANNdistArray d = new ANNdist[knn];

	for (size_t i = 0; i < points.size(); ++i) {
    vector<Point> neighbors;
		ANNpoint p = point_array[i];

		t.annkSearch(p, knn, n, d, eps);

    neighbors.push_back(points[i]);
		for (int j = 0; j < knn; ++j) {
			if ( n[j] != (int)i )
        neighbors.push_back(points[n[j]]);
		}

    Point centroid(0, 0, 0);
    for(size_t j = 0; j < neighbors.size(); ++j) {
      centroid.x += neighbors[j].x;
      centroid.y += neighbors[j].y;
      centroid.z += neighbors[j].z;
    }
    centroid.x /= (double) neighbors.size();
    centroid.y /= (double) neighbors.size();
    centroid.z /= (double) neighbors.size();

    

	}

	delete[] n;
	delete[] d;
}

int main(int argc, char **argv)
{
    // commandline arguments
    int start, end;
    bool scanserver;
    int maxDist, minDist;
    string dir;
    IOType iotype;
    int normalMethod;
    int knn;

    parse_options(argc, argv, start, end, scanserver, dir, iotype, maxDist, minDist, normalMethod, knn);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
      cerr << "No scans found. Did you use the correct format?" << endl;
      exit(-1);
    }

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
      Scan* scan = *it;

      // apply optional filtering
      scan->setRangeFilter(maxDist, minDist);
      
      
    }

    Scan::closeDirectory();

    return 0;
}
