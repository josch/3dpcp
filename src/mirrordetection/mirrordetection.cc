#include "slam6d/scan.h"
#include "slam6d/globals.icc"

#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <string>
#include <iostream>
#include <limits>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <sys/stat.h>
#include <sys/types.h>

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
        IOType &iotype, unsigned int &cluster_count,
        unsigned int &attempts)
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

    po::options_description cluster("Cluster options");
    cluster.add_options()
        ("count,c", po::value<unsigned int>(&cluster_count)->default_value(2),
         "the number of clusters to split the scan by")
        ("attempts,a", po::value<unsigned int>(&attempts)->default_value(5),
         "how many times the algorithm is executed using different initial "
         "labelings");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(cluster).add(hidden);

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
 * convert a scan to a matrix
 */
void convertToMat(Scan* scan, cv::Mat& scan_cv)
{
  DataXYZ xyz = scan->get("xyz");
  unsigned int nPoints = xyz.size();
  scan_cv.create(nPoints,1,CV_32FC(3));
  scan_cv = cv::Scalar::all(0);
  cv::MatIterator_<cv::Vec3f> it = scan_cv.begin<cv::Vec3f>();
  for(unsigned int i = 0; i < nPoints; i++, ++it) {
    (*it)[0] = xyz[i][0];
    (*it)[1] = xyz[i][1];
    (*it)[2] = xyz[i][2];
  }
}

/*
 * create a directory
 */
void createdirectory(string segdir)
{
    int success = mkdir(segdir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);

    if (success == 0 || errno == EEXIST) {
        cout << "Writing segmentations to " << segdir << endl;
    } else {
        cerr << "Creating directory " << segdir << " failed" << endl;
        exit(1);
    }
}
/*
 * given a sample and label matrix, write them to uos files
 */
void write3dfiles(Mat &samples, Mat &labels, unsigned int cluster_count, unsigned int nPoints, string &dir)
{
    vector<ofstream*> outfiles(cluster_count);
    for (unsigned int i = 0; i < cluster_count; i++) {
        std::stringstream outfilename;
        outfilename << dir << "/scan" << std::setw(3) << std::setfill('0') << i << ".3d";
        outfiles[i] = new ofstream(outfilename.str());
    }

    for (unsigned int i = 0; i < nPoints; ++i) {
        int cluster_idx = labels.at<int>(i, 0); // labels is 1 x nPoints
        cv::Vec3f entry = samples.at<cv::Vec3f>(i, 0);
        (*(outfiles[cluster_idx])) << entry(0) << " " << entry(1) << " " << entry(2) << endl;
    }

    for (unsigned i = 0; i < cluster_count; i++) {
        outfiles[i]->close();
    }
}

// write .pose files
// .frames files can later be generated from them using ./bin/pose2frames
void writeposefiles(int num, string &dir, const double* rPos, const double* rPosTheta)
{
    for (int i = 0; i < num; i++) {
        std::stringstream posefilename;
        posefilename << dir << "/scan" << std::setw(3) << std::setfill('0') << i << ".pose";
        ofstream posefile(posefilename.str());
        posefile << rPos[0] << " " << rPos[1] << " " << rPos[2] << endl;
        posefile << deg(rPosTheta[0]) << " " << deg(rPosTheta[1]) << " " << deg(rPosTheta[2]) << endl;
        posefile.close();
    }
}

/*
 * print the found cluster centers
 */
void printcenters(Mat &centers, unsigned int cluster_count)
{
    cout << "Found centers:" << endl;
    for (unsigned int i = 0; i < cluster_count; i++) {
        cout << "\tcenter #" << i << ": ("
            << centers.at<float>(i,0) << " " << centers.at<float>(i, 1)
            << " " << centers.at<float>(i, 2) << ")" << endl;
    }
}

int main(int argc, char **argv)
{
    int start, end;
    bool scanserver;
    int maxDist, minDist;
    string dir;
    IOType iotype;
    unsigned int cluster_count, attempts;

    parse_options(argc, argv, start, end, scanserver, dir, maxDist, minDist,
            iotype, cluster_count, attempts);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    Mat samples, labels, centers;
    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;
        scan->setRangeFilter(maxDist, minDist);

        string mirrordir = dir + "mirrors" + scan->getIdentifier();
        createdirectory(mirrordir);

        convertToMat(scan, samples);

        kmeans(samples, cluster_count, labels,
                TermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 10000, 0.0001),
                attempts, KMEANS_PP_CENTERS, centers);

        printcenters(centers, cluster_count);

        write3dfiles(samples, labels, cluster_count, scan->size<DataXYZ>("xyz"), mirrordir);

        writeposefiles(cluster_count, mirrordir, scan->get_rPos(), scan->get_rPosTheta());
     }
}
