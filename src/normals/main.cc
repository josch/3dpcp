/**
 * 
 * Copyright (C) Jacobs University Bremen
 *
 * @author Vaibhav Kumar Mehta
 * @file normals.cc
 */

#include "normals/normals.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

/*
 * validates normal calculation method specification
 */
void validate(boost::any& v, const std::vector<std::string>& values,
		    normal_method*, int) {
  if (values.size() == 0)
    throw std::runtime_error("Invalid model specification");
  string arg = values.at(0);
  if(strcasecmp(arg.c_str(), "AKNN") == 0) v = AKNN;
  else if(strcasecmp(arg.c_str(), "ADAPTIVE_AKNN") == 0) v = ADAPTIVE_AKNN;
  else if(strcasecmp(arg.c_str(), "PANORAMA") == 0) v = PANORAMA;
  else if(strcasecmp(arg.c_str(), "PANORAMA_FAST") == 0) v = PANORAMA_FAST;
  else throw std::runtime_error(std::string("normal calculation method ") + arg + std::string(" is unknown"));
}

/// validate IO types
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

/// Parse commandline options
void parse_options(int argc, char **argv, int &start, int &end, bool &scanserver, int &max_dist, int &min_dist, string &dir,
                   IOType &iotype, int &k1, int &k2, normal_method &ntype,int &width,int &height)
{
  /// ----------------------------------
  /// set up program commandline options
  /// ----------------------------------
  po::options_description cmd_options("Usage: calculateNormals <options> where options are (default values in brackets)");
  cmd_options.add_options()
    ("help,?", "Display this help message")
    ("start,s", po::value<int>(&start)->default_value(0), "Start at scan number <arg>")
    ("end,e", po::value<int>(&end)->default_value(-1), "Stop at scan number <arg>")
    ("scanserver,S", po::value<bool>(&scanserver)->default_value(false), "Use the scanserver as an input method")
    ("format,f", po::value<IOType>(&iotype)->default_value(UOS),
	"using shared library <arg> for input. (chose format from [uos|uosr|uos_map|"
	"uos_rgb|uos_frames|uos_map_frames|old|rts|rts_map|ifp|"
	"riegl_txt|riegl_rgb|riegl_bin|zahn|ply])")
    ("max,M", po::value<int>(&max_dist)->default_value(-1),"neglegt all data points with a distance larger than <arg> 'units")
    ("min,m", po::value<int>(&min_dist)->default_value(-1),"neglegt all data points with a distance smaller than <arg> 'units")
    ("normal,g", po::value<normal_method>(&ntype)->default_value(AKNN), "normal calculation method "
	"(AKNN, ADAPTIVE_AKNN, PANORAMA, PANORAMA_FAST)")
    ("K1,k", po::value<int>(&k1)->default_value(20), "<arg> value of K value used in the nearest neighbor search of ANN or" 								     "kmin for k-adaptation")
    ("K2,K", po::value<int>(&k2)->default_value(20), "<arg> value of Kmax for k-adaptation")
    ("width,w", po::value<int>(&width)->default_value(1280),"width of panorama image")
    ("height,h", po::value<int>(&height)->default_value(960),"height of panorama image")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-dir", po::value<string>(&dir), "input dir");

  po::positional_options_description pd;
  pd.add("input-dir", 1);

  po::options_description all;
  all.add(cmd_options).add(hidden);

  po::variables_map vmap;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), vmap);
  po::notify(vmap);

  if (vmap.count("help")) {
    cout << cmd_options << endl << endl;
    cout << "SAMPLE COMMAND FOR CALCULATING NORMALS" << endl;
    cout << " bin/calculateNormals -s 0 -e 0 -f UOS -g AKNN -k 20 dat/" <<endl;
    cout << endl << endl;
    cout << "SAMPLE COMMAND FOR VIEWING CALCULATING NORMALS IN RGB SPACE" << endl;
    cout << " bin/show -c -f UOS_RGB dat/normals/" << endl;
    exit(-1);
  }

  // read scan path
  if (dir[dir.length()-1] != '/') dir = dir + "/";

}

/*
 * retrieve a cv::Mat with x,y,z,r from a scan object
 * functionality borrowed from scan_cv::convertScanToMat but this function
 * does not allow a scanserver to be used, prints to stdout and can only
 * handle a single scan
 */
cv::Mat scan2mat(Scan *source)
{
  DataXYZ xyz = source->get("xyz");

  DataReflectance xyz_reflectance = source->get("reflectance");

  unsigned int nPoints = xyz.size();
  cv::Mat scan(nPoints,1,CV_32FC(4));
  scan = cv::Scalar::all(0);

  cv::MatIterator_<cv::Vec4f> it;

  it = scan.begin<cv::Vec4f>();
  for(unsigned int i = 0; i < nPoints; i++){
    float x, y, z, reflectance;
    x = xyz[i][0];
    y = xyz[i][1];
    z = xyz[i][2];
    if(xyz_reflectance.size() != 0)
	 {
	   reflectance = xyz_reflectance[i];
		
	   //normalize the reflectance
	   reflectance += 32;
	   reflectance /= 64;
	   reflectance -= 0.2;
	   reflectance /= 0.3;
	   if (reflectance < 0) reflectance = 0;
	   if (reflectance > 1) reflectance = 1;
	 }	
	
    (*it)[0] = x;
    (*it)[1] = y;
    (*it)[2] = z;
    if(xyz_reflectance.size() != 0)
	 (*it)[3] = reflectance;
    else
	 (*it)[3] = 0;	
	
    ++it;
  }
  return scan;
}
/*
 * convert a matrix of float values (range image) to a matrix of unsigned
 * eight bit characters using different techniques
 */
cv::Mat float2uchar(cv::Mat &source, bool logarithm, float cutoff)
{
  cv::Mat result(source.size(), CV_8U, cv::Scalar::all(0));
  float max = 0;
  // find maximum value
  if (cutoff == 0.0) {
    // without cutoff, just iterate through all values to find the largest
    for (cv::MatIterator_<float> it = source.begin<float>();
	    it != source.end<float>(); ++it) {
	 float val = *it;
	 if (val > max) {
	   max = val;
	 }
    }
  } else {
    // when a cutoff is specified, sort all the points by value and then
    // specify the max so that <cutoff> values are larger than it
    vector<float> sorted(source.cols*source.rows);
    int i = 0;
    for (cv::MatIterator_<float> it = source.begin<float>();
	    it != source.end<float>(); ++it, ++i) {
	 sorted[i] = *it;
    }
    std::sort(sorted.begin(), sorted.end());
    max = sorted[(int)(source.cols*source.rows*(1.0-cutoff))];
    cout << "A cutoff of " << cutoff << " resulted in a max value of " << max << endl;
  }

  cv::MatIterator_<float> src = source.begin<float>();
  cv::MatIterator_<uchar> dst = result.begin<uchar>();
  cv::MatIterator_<float> end = source.end<float>();
  if (logarithm) {
    // stretch values from 0 to max logarithmically over 0 to 255
    // using the logarithm allows to represent smaller values with more
    // precision and larger values with less
    max = log(max+1);
    for (; src != end; ++src, ++dst) {
	 float val = (log(*src+1)*255.0)/max;
	 if (val > 255)
	   *dst = 255;
	 else
	   *dst = (uchar)val;
    }
  } else {
    // stretch values from 0 to max linearly over 0 to 255
    for (; src != end; ++src, ++dst) {
	 float val = (*src*255.0)/max;
	 if (val > 255)
	   *dst = 255;
	 else
	   *dst = (uchar)val;
    }
  }
  return result;
}
/// Write a pose file with the specofied name
void writePoseFiles(string dir, const double* rPos, const double* rPosTheta,int scanNumber)
{
  string poseFileName = dir + "/scan" + to_string(scanNumber, 3) + ".pose";
  ofstream posout(poseFileName.c_str());

  posout << rPos[0] << " "
	    << rPos[1] << " "
	    << rPos[2] << endl
	    << deg(rPosTheta[0]) << " "
	    << deg(rPosTheta[1]) << " "
	    << deg(rPosTheta[2]) << endl;
  posout.clear();
  posout.close();
}

/// write scan files for all segments
void writeScanFiles(string dir, vector<Point> &points, vector<Point> &normals, int scanNumber)
{
  string ofilename = dir + "/scan" + to_string(scanNumber, 3) + ".3d";
  ofstream normptsout(ofilename.c_str());
	
  for (size_t i=0; i<points.size(); ++i)
    {
	 int r,g,b;
	 r = (int)(normals[i].x * (127.5) + 127.5);
	 g = (int)(normals[i].y * (127.5) + 127.5);
	 b = (int)(fabs(normals[i].z) * (255.0)); 
	 normptsout <<points[i].x<<" "<<points[i].y<<" "<<points[i].z<<" "<<r<<" "<<g<<" "<<b<<" "<<endl;
    }
  normptsout.clear();
  normptsout.close();
}

/// =============================================
/// Main
/// =============================================
int main(int argc, char** argv)
{
  int start, end;
  bool scanserver;
  int max_dist, min_dist;
  string dir;
  IOType iotype;
  int k1, k2;
  normal_method ntype;
  int width, height;

  parse_options(argc, argv, start, end, scanserver, max_dist, min_dist,
			 dir, iotype, k1, k2, ntype, width, height);

  /// ----------------------------------
  /// Prepare and read scans
  /// ----------------------------------
  if (scanserver) {
    try {
	 ClientInterface::create();
    } catch(std::runtime_error& e) {
	 cerr << "ClientInterface could not be created: " << e.what() << endl;
	 cerr << "Start the scanserver first." << endl;
	 exit(-1);
    }
  }

  /// Make directory for saving the scan segments
  string normdir = dir + "normals";

#ifdef _MSC_VER
  int success = mkdir(normdir.c_str());
#else
  int success = mkdir(normdir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
#endif
  if(success == 0) {
    cout << "Writing segments to " << normdir << endl;
  } else if(errno == EEXIST) {
    cout << "WARN: Directory " << normdir << " exists already. Contents will be overwriten" << endl;
  } else {
    cerr << "Creating directory " << normdir << " failed" << endl;
    exit(1);
  }

  /// Read the scans
  Scan::openDirectory(scanserver, dir, iotype, start, end);
  if(Scan::allScans.size() == 0) {
    cerr << "No scans found. Did you use the correct format?" << endl;
    exit(-1);
  }
	
  cv::Mat img;

  /// --------------------------------------------
  /// Initialize and perform segmentation
  /// --------------------------------------------
  std::vector<Scan*>::iterator it = Scan::allScans.begin();
  int scanNumber = 0;

  for( ; it != Scan::allScans.end(); ++it) {
    Scan* scan = *it;

    // apply optional filtering
    scan->setRangeFilter(max_dist, min_dist);
	
    const double* rPos = scan->get_rPos();
    const double* rPosTheta = scan->get_rPosTheta();

    /// read scan into points
    DataXYZ xyz(scan->get("xyz"));
    vector<Point> points;
    points.reserve(xyz.size());
    vector<Point> normals;
    normals.reserve(xyz.size());
	
    for(unsigned int j = 0; j < xyz.size(); j++) {
	 points.push_back(Point(xyz[j][0], xyz[j][1], xyz[j][2]));
    }
	
    if(ntype == AKNN)
	 calculateNormalsAKNN(normals,points, k1, rPos);
    else if(ntype == ADAPTIVE_AKNN)	
	 calculateNormalsAdaptiveAKNN(normals,points, k1, k2, rPos);
    else 
	 {
	   // create panorama
	   fbr::panorama fPanorama(width, height, fbr::EQUIRECTANGULAR, 1, 0, fbr::EXTENDED);
	   fPanorama.createPanorama(scan2mat(scan));
	
	   // the range image has to be converted from float to uchar
	   img = fPanorama.getRangeImage();
	   img = float2uchar(img, 0, 0.0);

	   if(ntype == PANORAMA)
		calculateNormalsPANORAMA(normals,points,fPanorama.getExtendedMap(), rPos);
	   else if(ntype == PANORAMA_FAST)
		cout << "PANORAMA_FAST is not working yet" << endl;
	   //   calculateNormalsFAST(normals,points,img,fPanorama.getExtendedMap());
	 }

    // pose file (repeated for the number of segments
    writePoseFiles(normdir, rPos, rPosTheta, scanNumber);
    // scan files for all segments
    writeScanFiles(normdir, points,normals,scanNumber);

    scanNumber++;
  }

  // shutdown everything
  if (scanserver)
    ClientInterface::destroy();
  else
    Scan::closeDirectory();

  cout << "Normal program end" << endl;

  return 0;
}

