#include "slam6d/scan.h"
#include "slam6d/globals.icc"

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <string>
#include <iostream>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace cv;
using namespace std;

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

void parse_options(int argc, char **argv, int &start, int &end,
        bool &scanserver, string &dir, int &maxDist, int &minDist,
        IOType &iotype)
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
        ("scanserver,S", po::value<bool>(&scanserver)->default_value(false),
         "Use the scanserver as an input method and handling of scan data");

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

void convertToMat(Scan* scan, cv::Mat& scan_cv) 
{
  DataXYZ xyz = scan->get("xyz");
  DataReflectance xyz_reflectance = scan->get("reflectance");
  unsigned int nPoints = xyz.size();
  scan_cv.create(nPoints,1,CV_32FC(3));
  scan_cv = cv::Scalar::all(0); 
  double zMax = numeric_limits<double>::min(); 
  double zMin = numeric_limits<double>::max();
  cv::MatIterator_<cv::Vec3f> it = scan_cv.begin<cv::Vec3f>();
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
    (*it)[0] = x;
    (*it)[1] = y;
    (*it)[2] = z;
    //(*it)[3] = reflectance;
    //finding min and max of z                                      
    if (z > zMax) zMax = z;
    if (z < zMin) zMin = z;
    ++it;
  }
}

int main(int argc, char **argv)
{
    int start, end;
    bool scanserver;
    int maxDist, minDist;
    string dir;
    IOType iotype;

    parse_options(argc, argv, start, end, scanserver, dir, maxDist, minDist,
            iotype);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;
        scan->setRangeFilter(maxDist, minDist);

        Mat samples;
        convertToMat(scan, samples); //samples is 1 x number_of_points

        int cluster_count = 2;
        Mat labels;
        int attempts = 5;
        Mat centers;
        kmeans(samples, cluster_count, labels, TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001), attempts, KMEANS_PP_CENTERS, centers );

        cout << "Samples size: " << samples.size().width << " x " << samples.size().height << endl;
        cout << "Centers size: " << centers.size().width << " x " << centers.size().height << endl;
        cout << "Labels size: " << labels.size().width << " x " << labels.size().height << endl;

        ofstream cluster1, cluster2; 
        cluster1.open( (dir + "scan100.3d").c_str() );
        cluster2.open( (dir + "scan101.3d").c_str() );
        for (int i = 0; i < samples.rows; ++i) {
          int cluster_idx = labels.at<int>(i, 0); // labels is 1 x number_of_points
          switch (cluster_idx) {
            case 0: 
              {         
                cv::Vec3f entry = samples.at<cv::Vec3f>(i, 0);   
                for (int j = 0; j < 3; ++j) {          
                  cluster1 << entry(j) << " ";
                }                          
                cluster1 << endl;
              }
              break;
            case 1:
              {         
                cv::Vec3f entry = samples.at<cv::Vec3f>(i, 0);   
                for (int j = 0; j < 3; ++j) {          
                  cluster2 << entry(j) << " ";
                }                          
                cluster2 << endl;
              }
              break;
            default:
              break;
          }
        }
        cluster1.close();
        cluster2.close();
     }
}
