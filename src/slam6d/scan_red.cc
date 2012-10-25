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

void reduction_option_dependency(const po::variables_map & vm, reduction_method stype, const char *option)
{
    if (vm.count("reduction") && vm["reduction"].as<reduction_method>() == stype) {
        if (!vm.count(option)) {
            throw std::logic_error (string("this reduction option needs ")+option+" to be set");
        }
    }
}

void reduction_option_conflict(const po::variables_map & vm, reduction_method stype, const char *option)
{
    if (vm.count("reduction") && vm["reduction"].as<reduction_method>() == stype) {
        if (vm.count(option)) {
            throw std::logic_error (string("this reduction option is incompatible with ")+option);
        }
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
        double &voxel, int &octree, bool &use_reflectance)
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
        ("reflectance,R", po::bool_switch(&use_reflectance),
         "Use reflectance when reducing points and save scan files in UOSR format");

    po::options_description reduction("Reduction options");
    reduction.add_options()
        ("reduction,r", po::value<reduction_method>(&rtype)->required(),
         "choose reduction method (OCTREE, RANGE, INTERPOLATE)")
        ("scale,S", po::value<double>(&scale),
         "scaling factor")
        ("voxel,v", po::value<double>(&voxel),
         "voxel size")
        ("projection,P", po::value<fbr::projection_method>(&ptype),
         "projection method or panorama image")
        ("octree,O", po::value<int>(&octree),
         "0 -> center\n1 -> random\nN>1 -> random N")
        ("width,w", po::value<int>(&width),
         "width of panorama")
        ("height,h", po::value<int>(&height),
         "height of panorama");

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

    // display help
    if (vm.count("help")) {
        cout << cmdline_options;
        cout << endl
            << "Example usage:" << endl
            << "\t./bin/scan_red -s 0 -e 0 -f uos --reduction OCTREE --voxel 10 --octree 0 dat" << endl
            << "\t./bin/scan_red -s 0 -e 0 -f uos --reduction RANGE --scale 0.5 --projection EQUIRECTANGULAR --width 3600 --height 1000 dat" << endl
            << "\t./bin/scan_red -s 0 -e 0 -f uos --reduction INTERPOLATE --scale 0.2 --projection EQUIRECTANGULAR --width 3600 --height 1000 dat" << endl;
        exit(0);
    }

    po::notify(vm);

    reduction_option_dependency(vm, OCTREE, "voxel");
    reduction_option_dependency(vm, OCTREE, "octree");
    reduction_option_conflict(vm, OCTREE, "scale");
    reduction_option_conflict(vm, OCTREE, "projection");
    reduction_option_conflict(vm, OCTREE, "width");
    reduction_option_conflict(vm, OCTREE, "height");

    reduction_option_conflict(vm, RANGE, "voxel");
    reduction_option_conflict(vm, RANGE, "octree");
    reduction_option_dependency(vm, RANGE, "scale");
    reduction_option_dependency(vm, RANGE, "projection");
    reduction_option_dependency(vm, RANGE, "width");
    reduction_option_dependency(vm, RANGE, "height");

    reduction_option_conflict(vm, INTERPOLATE, "voxel");
    reduction_option_conflict(vm, INTERPOLATE, "octree");
    reduction_option_dependency(vm, INTERPOLATE, "scale");
    reduction_option_dependency(vm, INTERPOLATE, "projection");
    reduction_option_dependency(vm, INTERPOLATE, "width");
    reduction_option_dependency(vm, INTERPOLATE, "height");

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

void scan2mat(Scan *source, cv::Mat &mat)
{
    DataXYZ xyz = source->get("xyz");
    unsigned int nPoints = xyz.size();
    mat.create(nPoints,1,CV_32FC(4));
    mat = cv::Scalar::all(0);
    cv::MatIterator_<cv::Vec4f> it = mat.begin<cv::Vec4f>();
    for(unsigned int i = 0; i < nPoints; i++){
        (*it)[0] = xyz[i][0];
        (*it)[1] = xyz[i][1];
        (*it)[2] = xyz[i][2];
        ++it;
    }
}

void reduce_octree(Scan *scan, vector<cv::Vec3f> &reduced_points, int octree,
        int red, bool use_reflectance)
{
    if (use_reflectance) {
      unsigned int types = PointType::USE_REFLECTANCE;
      PointType pointtype(types);
      scan->setReductionParameter(red, octree, pointtype);
    } else {
      scan->setReductionParameter(red, octree);
    }
    DataXYZ xyz_r(scan->get("xyz reduced"));

    cout << red << " " << octree  << " " << xyz_r.size() << " " << endl;
    for(unsigned int j = 0; j < xyz_r.size(); j++) {
        reduced_points.push_back(cv::Vec3d(xyz_r[j][0], xyz_r[j][1], xyz_r[j][2]));
    }
}

void reduce_range(Scan *scan, string reddir, string id, int width, int height,
        fbr::projection_method ptype, double scale, bool use_reflectance)
{
    panorama image(width, height, ptype);
    cv::Mat mat;
    scan2mat(scan, mat);
    image.createPanorama(mat);
    image.getDescription();

    /// Resize the range image, specify desired interpolation method
    cv::Mat range_image_resized; // reflectance_image_resized;
    resize(image.getRangeImage(), range_image_resized, cv::Size(),
            scale, scale, cv::INTER_NEAREST);
    // why does the panorama class have functionality to write a UOS file?
    // XXX: change recoverPointCloud to just return a point cloud
    image.recoverPointCloud(range_image_resized, reddir + "/scan" + id + ".3d");
}

void reduce_interpolation(Scan *scan,
        vector<cv::Vec3f> &reduced_points, int width, int height,
        fbr::projection_method ptype, double scale, bool use_reflectance)
{
    panorama image(width, height, ptype);
    cv::Mat mat;
    scan2mat(scan, mat);
    image.createPanorama(mat);
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
    bool use_reflectance;

    parse_options(argc, argv, start, end, scanserver, width, height, ptype,
            dir, iotype, maxDist, minDist, rtype, scale, voxel, octree, use_reflectance);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;

        scan->setRangeFilter(maxDist, minDist);

        vector<cv::Vec3f> reduced_points;

        string reddir = dir + "reduced";
        createdirectory(reddir);

        switch (rtype) {
            case OCTREE:
                reduce_octree(scan, reduced_points, octree, voxel, use_reflectance);
                write3dfile(reduced_points, reddir, scan->getIdentifier());
                break;
            case RANGE:
                reduce_range(scan, reddir, scan->getIdentifier(), width, height, ptype, scale, use_reflectance);
                break;
            case INTERPOLATE:
                reduce_interpolation(scan, reduced_points, width, height, ptype, scale, use_reflectance);
                write3dfile(reduced_points, reddir, scan->getIdentifier());
                break;
            default:
                cerr << "unknown method" << endl;
                return 1;
                break;
        }

        writeposefile(reddir, scan->get_rPos(), scan->get_rPosTheta(), scan->getIdentifier());
    }
}
