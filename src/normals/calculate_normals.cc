/*
 * calculateNormals implementation
 *
 * Copyright (C) Johannes Schauer, Razvan Mihaly
 *
 * Released under the GPL version 3
 *
 */

#include "ANN/ANN.h"

#include "newmat/newmat.h"
#include "newmat/newmatap.h"
#include "newmat/newmatio.h"
using namespace NEWMAT;

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
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
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

void mapNormalToRGB(const Point& normal, Point& rgb) 
{
    rgb.x = 127.5 * normal.x + 127.5;
    rgb.y = 127.5 * normal.y + 127.5;
    rgb.z = 255.0 * fabs(normal.z);
} 

void writeNormals(const Scan* scan, const string& dir, 
                const vector<Point>& points, const vector<Point>& normals) 
{

    stringstream ss; 
    ss << dir << "scan" << string(scan->getIdentifier()) << ".3d"; 
    ofstream scan_file;
    scan_file.open(ss.str().c_str());
    for(size_t i = 0;  i < points.size(); ++i) {
        Point rgb;
        mapNormalToRGB(normals[i], rgb);
        scan_file << points[i].x << " " << points[i].y << " " << points[i].z << " "
                    << (unsigned int) rgb.x << " " << (unsigned int) rgb.y << " " << (unsigned int) rgb.z << "\n";  
    }
    scan_file.close();

    ss.clear(); ss.str(string());
    ss << dir << "scan" << string(scan->getIdentifier()) << ".pose";
    ofstream pose_file; 
    pose_file.open(ss.str().c_str());
    pose_file << 0 << " " << 0 << " " << 0 << "\n" << 0 << " " << 0 << " " << 0 << "\n";
    pose_file.close();
}

void computeNeighbors(const Scan* scan, const vector<Point>& points, vector<Point>& normals, int knn, double eps=1.0) 
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

    ColumnVector origin(3);
    const double *scan_pose = scan->get_rPos();
    for (int i = 0; i < 3; ++i)
        origin(i+1) = scan_pose[i];

    for (size_t i = 0; i < points.size(); ++i) {
        ColumnVector point_vector(3);
        point_vector(1) = points[i].x - origin(1);
        point_vector(2) = points[i].y - origin(2);
        point_vector(3) = points[i].z - origin(3);
        point_vector = point_vector / point_vector.NormFrobenius();

        vector<Point> neighbors;
        ANNpoint p = point_array[i];

        t.annkSearch(p, knn, n, d, eps);

        neighbors.push_back(points[i]);
        for (int j = 0; j < knn; ++j) {
            if ( n[j] != (int)i )
                neighbors.push_back(points[n[j]]);
        }

        Point centroid(0, 0, 0);
        for (size_t j = 0; j < neighbors.size(); ++j) {
            centroid.x += neighbors[j].x;
            centroid.y += neighbors[j].y;
            centroid.z += neighbors[j].z;
        }
        centroid.x /= (double) neighbors.size();
        centroid.y /= (double) neighbors.size();
        centroid.z /= (double) neighbors.size();

        Matrix S(3, 3);
        S = 0.0;
        for (size_t j = 0; j < neighbors.size(); ++j) {
            ColumnVector point_prime(3);
            point_prime(1) = neighbors[j].x - centroid.x;
            point_prime(2) = neighbors[j].y - centroid.y;
            point_prime(3) = neighbors[j].z - centroid.z;
            S = S + point_prime * point_prime.t();
        }
        // normalize S
        for (int j = 0; j < 3; ++j) 
            for (int k = 0; k < 3; ++k)
                S(j+1, k+1) /= (double) neighbors.size();

        SymmetricMatrix C;
        C << S;
        // compute eigendecomposition of C
        Matrix V(3,3); // for eigenvectors
        DiagonalMatrix D(3);   // for eigenvalues
        // the decomposition
        Jacobi(C, D, V);

#ifdef DEBUG       
        // Print the result
        cout << "The eigenvalues matrix:" << endl;
        cout << D << endl;
        cout << "The eigenvectors matrix:" << endl;
        cout << V << endl << endl;
#endif 
        ColumnVector v1(3);
        v1(1) = V(1,1);
        v1(2) = V(2,1);
        v1(3) = V(3,1);
        // consider first (smallest) eigenvector as the normal
        Real angle = (v1.t() * point_vector).AsScalar();

        // orient towards scan pose TODO: double check this
        if (angle < 0) {
            v1 *= -1.0;
        }
        normals.push_back( Point(v1(1), v1(2), v1(3)) );
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

    boost::filesystem::path boost_dir(dir + "normals/");
    boost::filesystem::create_directory(boost_dir);

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        vector<Point> points, normals;
        Scan* scan = *it;

        // apply optional filtering
        scan->setRangeFilter(maxDist, minDist);

        DataXYZ xyz = scan->get("xyz");
        DataReflectance xyz_reflectance = scan->get("reflectance");
        unsigned int nPoints = xyz.size();
        for(unsigned int i = 0; i < nPoints; ++i) {
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

            points.push_back(Point(x, y, z));
        }

        computeNeighbors(scan, points, normals, knn);
        writeNormals(scan, dir + "normals/", points, normals);
    }

    Scan::closeDirectory();

    return 0;
}
