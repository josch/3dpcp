/*
 * scan_red implementation
 *
 * Copyright (C) Dorit Borrmann, Razvan-George Mihalyi, Remus Dumitru 
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

#include "slam6d/metaScan.h"
#include "slam6d/io_utils.h"
#include "slam6d/scan.h"
#include "slam6d/fbr/fbr_global.h"
#include "slam6d/fbr/panorama.h"
#include "slam6d/fbr/scan_cv.h"

#include "scanserver/clientInterface.h"

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


#define IMAGE_HEIGHT 1000
#define IMAGE_WIDTH 3600

using namespace fbr;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

enum reduction_method {OCTREE, RANGE, INTERPOLATE};

/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const po::variables_map & vm,
                         const char *opt1, const char *opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted()
        && vm.count(opt2) && !vm[opt2].defaulted())
        throw std::logic_error(string("Conflicting options '")
                               + opt1 + "' and '" + opt2 + "'.");
}

/* Function used to check that if 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const po::variables_map & vm,
                       const char *for_what, const char *required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0
            || vm[required_option].defaulted())
            throw std::logic_error(string("Option '") + for_what +
                                   "' requires option '" +
                                   required_option + "'.");
}

/*
 * validates panorama method specification
 */
namespace fbr {
    void validate(boost::any& v, const std::vector<std::string>& values,
                  projection_method*, int) {
        if (values.size() == 0)
            throw std::runtime_error("Invalid model specification");
        string arg = values.at(0);
        if(strcasecmp(arg.c_str(), "EQUIRECTANGULAR") == 0) v = EQUIRECTANGULAR;
        else if(strcasecmp(arg.c_str(), "CYLINDRICAL") == 0) v = CYLINDRICAL;
        else if(strcasecmp(arg.c_str(), "MERCATOR") == 0) v = MERCATOR;
        else if(strcasecmp(arg.c_str(), "RECTILINEAR") == 0) v = RECTILINEAR;
        else if(strcasecmp(arg.c_str(), "PANNINI") == 0) v = PANNINI;
        else if(strcasecmp(arg.c_str(), "STEREOGRAPHIC") == 0) v = STEREOGRAPHIC;
        else if(strcasecmp(arg.c_str(), "ZAXIS") == 0) v = ZAXIS;
        else if(strcasecmp(arg.c_str(), "CONIC") == 0) v = CONIC;
        else throw std::runtime_error(std::string("projection method ") + arg + std::string(" is unknown"));
    }
}

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
 * validates reduction method specification
 */
void validate(boost::any& v, const std::vector<std::string>& values,
        reduction_method*, int) {
    if (values.size() == 0)
        throw std::runtime_error("Invalid model specification");
    string arg = values.at(0);
    if(strcasecmp(arg.c_str(), "OCTREE") == 0) v = OCTREE;
    else if(strcasecmp(arg.c_str(), "RANGE") == 0) v = RANGE;
    else if(strcasecmp(arg.c_str(), "INTERPOLATE") == 0) v = INTERPOLATE;
    else throw std::runtime_error(std::string("reduction method ") + arg + std::string(" is unknown"));
}

void parse_options(int argc, char **argv, int &start, int &end,
        bool &scanserver, int &width, int &height,
        fbr::projection_method &ptype, string &dir, IOType &iotype,
        int &maxDist, int &minDist, reduction_method &rtype, double &scale,
        double &voxel, int &octree)
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
         "Use the scanserver as an input method and handling of scan data");

    po::options_description reduction("Reduction options");
    reduction.add_options()
        ("reduction,r", po::value<reduction_method>(&rtype),
         "choose reduction method (OCTREE, RANGE, INTERPOLATE)")
        ("scale,S", po::value<double>(&scale),
         "scaling factor")
        ("voxel,v", po::value<double>(&voxel),
         "voxel size")
        ("projection,P", po::value<fbr::projection_method>(&ptype),
         "projection method or panorama image")
        ("octree,O", po::value<int>(&octree));

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(reduction).add(hidden);

    // options visible with --help
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(reduction).add(input);

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
        cout << "\nExample usage:\n" << endl;
        exit(0);
    }

#ifndef _MSC_VER
    if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
    if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif
}

void createdirectory(string dir)
{
    int success = mkdir(dir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);

    if (success == 0 || errno == EEXIST) {
        cout << "Writing to " << dir << endl;
    } else {
        cerr << "Creating directory " << dir << " failed" << endl;
        exit(1);
    }
}

cv::Mat scan2mat(Scan *source)
{
    DataXYZ xyz = source->get("xyz");
    unsigned int nPoints = xyz.size();
    cv::Mat scan(nPoints,1,CV_32FC(4));
    scan = cv::Scalar::all(0);
    cv::MatIterator_<cv::Vec4f> it;
    it = scan.begin<cv::Vec4f>();
    for(unsigned int i = 0; i < nPoints; i++){
        (*it)[0] = xyz[i][0];
        (*it)[1] = xyz[i][1];
        (*it)[2] = xyz[i][2];
        ++it;
    }
    return scan;
}

void reduce_octree(Scan *scan, vector<cv::Vec3f> &reduced_points, int octree,
        int red, double scale)
{
    scan->setReductionParameter(red, octree);
    DataXYZ xyz_r(scan->get("xyz reduced"));
    unsigned int nPoints = xyz_r.size();

    for(unsigned int j = 0; j < nPoints; j++) {
        cv::Vec3f vec(xyz_r[j][0], xyz_r[j][1], xyz_r[j][2]);
        reduced_points.push_back(vec);
    }
}

void reduce_range(Scan *scan, vector<cv::Vec3f> &reduced_points,
        int width, int height, fbr::projection_method ptype, double scale)
{
    panorama image(width, height, ptype);
    image.createPanorama(scan2mat(scan));
    image.getDescription();

    /// Resize the range image, specify desired interpolation method
    cv::Mat range_image_resized; // reflectance_image_resized;
    resize(image.getMap(), range_image_resized, cv::Size(), 
            scale, scale, cv::INTER_NEAREST);
    for(int i = 0; i < range_image_resized.rows; i++) {
        for(int j = 0; j < range_image_resized.cols; j++) {
            cv::Vec3f vec = range_image_resized.at<cv::Vec3f>(i, j);
            reduced_points.push_back(vec);
        }
    }
}

void reduce_interpolation(Scan *scan,
        vector<cv::Vec3f> &reduced_points, int width, int height,
        fbr::projection_method ptype, double scale)
{
    panorama image(width, height, ptype);
    image.createPanorama(scan2mat(scan));
    image.getDescription();

    /// Resize the range image, specify desired interpolation method
    cv::Mat range_image_resized; // reflectance_image_resized;
    string ofilename;
    stringstream ss;
    resize(image.getRangeImage(), range_image_resized, cv::Size(),
            scale, scale, cv::INTER_NEAREST);
    for(int i = 0; i < range_image_resized.rows; i++) {
        for(int j = 0; j < range_image_resized.cols; j++) {
            cv::Vec3f vec = range_image_resized.at<cv::Vec3f>(i, j);
            reduced_points.push_back(vec);
        }
    }
}

/*
 * given a vector of 3d points, write them out as uos files
 */
void write3dfile(vector<cv::Vec3f> &points, string &dir, string id)
{
    ofstream outfile(dir + "/scan" + id + ".3d");

    for (vector<cv::Vec3f>::iterator it=points.begin(); it < points.end(); it++) {
        outfile << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
    }

    outfile.close();
}

// write .pose files
// .frames files can later be generated from them using ./bin/pose2frames
void writeposefile(string &dir, const double* rPos, const double* rPosTheta, string id)
{
    ofstream posefile(dir + "/scan" + id + ".pose");
    posefile << rPos[0] << " " << rPos[1] << " " << rPos[2] << endl;
    posefile << deg(rPosTheta[0]) << " "
          << deg(rPosTheta[1]) << " "
          << deg(rPosTheta[2]) << endl;
    posefile.close();
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
    int start, end;
    bool scanserver;
    int width, height;
    int maxDist, minDist;
    fbr::projection_method ptype;
    string dir;
    IOType iotype;
    reduction_method rtype;
    double scale, voxel;
    int octree;

    parse_options(argc, argv, start, end, scanserver, width, height, ptype,
            dir, iotype, maxDist, minDist, rtype, scale, voxel, octree);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;

        scan->setRangeFilter(maxDist, minDist);

        vector<cv::Vec3f> reduced_points;

        switch (rtype) {
            case OCTREE:
                reduce_octree(scan, reduced_points, voxel, octree, scale);
                break;
            case RANGE:
                reduce_range(scan, reduced_points, width, height, ptype, scale);
                break;
            case INTERPOLATE:
                reduce_interpolation(scan, reduced_points, width, height, ptype, scale);
                break;
            default:
                break;
        }

        string reddir = dir + "reduced";
        createdirectory(reddir);

        write3dfile(reduced_points, reddir, scan->getIdentifier());
        writeposefile(reddir, scan->get_rPos(), scan->get_rPosTheta(), scan->getIdentifier());
    }
}
