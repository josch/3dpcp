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
#include "normals/pointNeighbor.h"
#include "slam6d/scan.h"
#include "slam6d/globals.icc"
#include "slam6d/fbr/panorama.h"

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

enum normal_method {KNN_PCA, AKNN_PCA, PANO_PCA, PANO_SRI};

void normal_option_dependency(const po::variables_map & vm, normal_method ntype, const char *option)
{
    if (vm.count("normalMethod") && vm["normalMethod"].as<normal_method>() == ntype) {
        if (!vm.count(option)) {
            throw std::logic_error (string("this normal method needs ")+option+" to be set");
        }
    }
}

void normal_option_conflict(const po::variables_map & vm, normal_method ntype, const char *option)
{
    if (vm.count("normalMethod") && vm["normalMethod"].as<normal_method>() == ntype) {
        if (vm.count(option)) {
            throw std::logic_error (string("this normal method is incompatible with ")+option);
        }
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

void validate(boost::any& v, const std::vector<std::string>& values,
        normal_method*, int) {
    if (values.size() == 0)
        throw std::runtime_error("Invalid model specification");
    string arg = values.at(0);
    if(strcasecmp(arg.c_str(), "KNN_PCA") == 0) v = KNN_PCA;
    else if(strcasecmp(arg.c_str(), "AKNN_PCA") == 0) v = AKNN_PCA;
    else if(strcasecmp(arg.c_str(), "PANO_PCA") == 0) v = PANO_PCA;
    else if(strcasecmp(arg.c_str(), "PANO_SRI") == 0) v = PANO_SRI;
    else throw std::runtime_error(std::string("normal method ") + arg + std::string(" is unknown"));
}

/*
 * parse commandline options, fill arguments
 */
void parse_options(int argc, char **argv, int &start, int &end,
        bool &scanserver, string &dir, IOType &iotype,
        int &maxDist, int &minDist, normal_method &normalMethod, int &knn,
        int &kmin, int &kmax, double& alpha, int &width, int &height,
        bool &flipnormals)
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
        ;

    po::options_description normal("Normal options");
    normal.add_options()
        ("normalMethod,N", po::value<normal_method>(&normalMethod)->default_value(KNN_PCA),
         "choose the method for computing normals:\n"
         "KNN_PCA  -- use kNN and PCA\n"
         "AKNN_PCA -- use adaptive kNN and PCA\n"
         "PANO_PCA -- use panorama image neighbors and PCA\n"
         "PANO_SRI -- use panorama image neighbors and spherical range image differentiation\n")
        ("knn,K", po::value<int>(&knn),
         "select the k in kNN search")
        ("kmin,1", po::value<int>(&kmin),
         "select k_min in adaptive kNN search")
        ("kmax,2", po::value<int>(&kmax),
         "select k_max in adaptive kNN search")
        ("alpha,a", po::value<double>(&alpha),
         "select the alpha parameter for detecting an ill-conditioned neighborhood")
        ("width,w", po::value<int>(&width),
         "width of panorama")
        ("height,h", po::value<int>(&height),
         "height of panorama")
        ("flipnormals,F", po::bool_switch(&flipnormals),
         "flip orientation of normals towards scan pose")
        ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(normal).add(hidden);

    // options visible with --help
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input).add(normal);

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

    normal_option_dependency(vm, KNN_PCA, "knn");
    normal_option_conflict(vm, KNN_PCA, "kmin");
    normal_option_conflict(vm, KNN_PCA, "kmax");
    normal_option_conflict(vm, KNN_PCA, "alpha");
    normal_option_conflict(vm, KNN_PCA, "width");
    normal_option_conflict(vm, KNN_PCA, "height");

    normal_option_conflict(vm, AKNN_PCA, "knn");
    normal_option_dependency(vm, AKNN_PCA, "kmin");
    normal_option_dependency(vm, AKNN_PCA, "kmax");
    normal_option_dependency(vm, AKNN_PCA, "alpha");
    normal_option_conflict(vm, AKNN_PCA, "width");
    normal_option_conflict(vm, AKNN_PCA, "height");

    normal_option_conflict(vm, PANO_PCA, "knn");
    normal_option_conflict(vm, PANO_PCA, "kmin");
    normal_option_conflict(vm, PANO_PCA, "kmax");
    normal_option_conflict(vm, PANO_PCA, "alpha");
    normal_option_dependency(vm, PANO_PCA, "width");
    normal_option_dependency(vm, PANO_PCA, "height");

    normal_option_conflict(vm, PANO_SRI, "knn");
    normal_option_conflict(vm, PANO_SRI, "kmin");
    normal_option_conflict(vm, PANO_SRI, "kmax");
    normal_option_conflict(vm, PANO_SRI, "alpha");
    normal_option_dependency(vm, PANO_SRI, "width");
    normal_option_dependency(vm, PANO_SRI, "height");

    // add trailing slash to directory if not present yet
    if (dir[dir.length()-1] != '/') dir = dir + "/";
}

/*
 * retrieve a cv::Mat with x,y,z,r from a scan object
 * functionality borrowed from scan_cv::convertScanToMat but this function
 * does not allow a scanserver to be used, prints to stdout and can only
 * handle a single scan
 */
void scan2mat(Scan* scan, cv::Mat& scan_cv) {
    DataXYZ xyz = scan->get("xyz");
    unsigned int nPoints = xyz.size();
    scan_cv.create(nPoints,1,CV_32FC(4));
    scan_cv = cv::Scalar::all(0);
    cv::MatIterator_<cv::Vec3f> it = scan_cv.begin<cv::Vec3f>();
    for(unsigned int i = 0; i < nPoints; i++){
        (*it)[0] = xyz[i][0];
        (*it)[1] = xyz[i][1];
        (*it)[2] = xyz[i][2];
        ++it;
    }
}

/**
 * Helper function that maps x, y, z to R, G, B using a linear function
 */
void mapNormalToRGB(const Point& normal, Point& rgb)
{
    rgb.x = 127.5 * normal.x + 127.5;
    rgb.y = 127.5 * normal.y + 127.5;
    rgb.z = 255.0 * fabs(normal.z);
}

/**
 * Write normals to .3d files using the uos_rgb format
 */
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
        scan_file
            << points[i].x << " " << points[i].y << " " << points[i].z << " "
            << (unsigned int) rgb.x << " " << (unsigned int) rgb.y << " "
            << (unsigned int) rgb.z << "\n";
    }
    scan_file.close();

    ss.clear(); ss.str(string());
    ss << dir << "scan" << string(scan->getIdentifier()) << ".pose";
    ofstream pose_file;
    pose_file.open(ss.str().c_str());
    pose_file << 0 << " " << 0 << " " << 0 << "\n" << 0 << " " << 0 << " " << 0 << "\n";
    pose_file.close();
}

/**
 * Compute eigen decomposition of a point and its neighbors using the NEWMAT library
 * @param point - input points with corresponding neighbors
 * @param e_values - out parameter returns the eigenvalues
 * @param e_vectors - out parameter returns the eigenvectors
 */
void computeEigenDecomposition(const PointNeighbor& point,
        DiagonalMatrix& e_values, Matrix& e_vectors)
{
    Point centroid(0, 0, 0);
    vector<Point> neighbors = point.neighbors;

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
    // the decomposition
    Jacobi(C, e_values, e_vectors);

#ifdef DEBUG
    // Print the result
    cout << "The eigenvalues matrix:" << endl;
    cout << e_values << endl;
#endif
}

/**
 * Compute neighbors using kNN search
 * @param points - input set of points
 * @param points_neighbors - output set of points with corresponding neighbors
 * @param knn - k constant in kNN search
 * @param kmax - to be used in adaptive knn search as the upper bound on adapting the k constant, defaults to -1 for regular kNN search
 * @param alpha - to be used in adaptive knn search for detecting ill-conditioned neighborhoods
 * @param eps - parameter required by the ANN library in kNN search
 */
void computeKNearestNeighbors(const vector<Point>& points,
        vector<PointNeighbor>& points_neighbors, int knn, int kmax=-1,
        double alpha=1000.0, double eps=1.0)
{
    ANNpointArray point_array = annAllocPts(points.size(), 3);
    for (size_t i = 0; i < points.size(); ++i) {
        point_array[i] = new ANNcoord[3];
        point_array[i][0] = points[i].x;
        point_array[i][1] = points[i].y;
        point_array[i][2] = points[i].z;
    }

    ANNkd_tree t(point_array, points.size(), 3);
    ANNidxArray n;
    ANNdistArray d;

    if (kmax < 0) {
        /// regular kNN search, allocate memory for knn
        n = new ANNidx[knn];
        d = new ANNdist[knn];
    } else {
        /// adaptive kNN search, allocate memory for kmax
        n = new ANNidx[kmax];
        d = new ANNdist[kmax];
    }

    for (size_t i = 0; i < points.size(); ++i) {
        vector<Point> neighbors;
        ANNpoint p = point_array[i];

        t.annkSearch(p, knn, n, d, eps);

        neighbors.push_back(points[i]);
        for (int j = 0; j < knn; ++j) {
            if ( n[j] != (int)i )
                neighbors.push_back(points[n[j]]);
        }

        PointNeighbor current_point(points[i], neighbors);
        points_neighbors.push_back( current_point );
        Matrix e_vectors(3,3); e_vectors = 0.0;
        DiagonalMatrix e_values(3); e_values = 0.0;
        computeEigenDecomposition( current_point, e_values, e_vectors );

        if (kmax > 0) {
            /// detecting an ill-conditioned neighborhood
            if (e_values(3) / e_values(2) > alpha && e_values(2) > 0.0) {
                if (knn < kmax)
                    cout << "Increasing kmin to " << ++knn << endl;
            }
        }
    }

    delete[] n;
    delete[] d;
}

/**
 * Compute neighbors using panorama images
 * @param fPanorama - input panorama image created from the current scan
 * @param points_neighbors - output set of points with corresponding neighbors
 */
void computePanoramaNeighbors(Scan* scan,
        vector<PointNeighbor>& points_neighbors, int width, int height)
{
    cv::Mat scan_cv;
    scan2mat(scan, scan_cv);
    fbr::panorama fPanorama(width, height, fbr::EQUIRECTANGULAR, 1, 0, fbr::EXTENDED);
    fPanorama.createPanorama(scan_cv);
    cv::Mat img = fPanorama.getRangeImage();
    vector<vector<vector<cv::Vec3f> > > extended_map = fPanorama.getExtendedMap();
    for (int row = 0; row < height; ++row) {
        for (int col = 0; col < width; ++col) {
            vector<cv::Vec3f> points_panorama = extended_map[row][col];
            /// if no points found, skip pixel
            if (points_panorama.size() < 1) continue;
            /// foreach point from panorama consider all the rest as its neighbors
            for (size_t point_idx = 0; point_idx < points_panorama.size(); ++point_idx) {
                Point point;
                point.x = points_panorama[point_idx][0];
                point.y = points_panorama[point_idx][1];
                point.z = points_panorama[point_idx][2];
                vector<Point> neighbors;
                for (size_t i = 0; i < points_panorama.size(); ++i) {
                    if (i != point_idx)
                        neighbors.push_back(Point (points_panorama[i][0], points_panorama[i][1], points_panorama[i][2]) );
                }
                /// compute neighbors by examining adjacent pixels
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        if (!(i==0 && j==0) && !(row+i < 0 || col+j < 0)
                                &&  !(row+i >= height || col+j >= width) ) {
                            vector<cv::Vec3f> neighbors_panorama = extended_map[row+i][col+j];
                            for (size_t k = 0; k < neighbors_panorama.size(); ++k)
                                neighbors.push_back(Point (neighbors_panorama[k][0],
                                            neighbors_panorama[k][1],
                                            neighbors_panorama[k][2]) );
                        }
                    }
                }
                PointNeighbor current_point(point, neighbors);
                // add neighbors from range image
                for (int i = -1; i <= 1; ++i) {
                    for (int j = -1; j <= 1; ++j) {
                        if (!(i==0 && j==0) && !(row+i < 0 || col+j < 0)
                                && !(row+i >= height || col+j >= width)) {
                            current_point.range_neighbors[1+i][1+j] = img.at<float>(row+i,col+j);
                        }
                    }
                }
                current_point.range_image_row = row;
                current_point.range_image_col = col;
                points_neighbors.push_back(current_point);
            }
        }
    }
}

/**
 * Compute normals using PCA given a set of points and their neighbors
 * @param scan - pointer to current scan, used to compute the position vectors
 * @param points - input set of points with corresponding neighbors
 * @param normals - output set of normals
 */
void computePCA(const Scan* scan, const vector<PointNeighbor>& points,
        vector<Point>& normals, bool flipnormals)
{
    ColumnVector origin(3);
    const double *scan_pose = scan->get_rPos();
    for (int i = 0; i < 3; ++i)
        origin(i+1) = scan_pose[i];

    for(size_t i = 0; i < points.size(); ++i) {
        vector<Point> neighbors = points[i].neighbors;

        if (points[i].neighbors.size() < 2) {
            normals.push_back( Point(0,0,0) );
            continue;
        }

        ColumnVector point_vector(3);
        point_vector(1) = points[i].point.x - origin(1);
        point_vector(2) = points[i].point.y - origin(2);
        point_vector(3) = points[i].point.z - origin(3);
        point_vector = point_vector / point_vector.NormFrobenius();

        Matrix e_vectors(3,3); e_vectors = 0.0;
        DiagonalMatrix e_values(3); e_values = 0.0;
        computeEigenDecomposition(points[i], e_values, e_vectors);

        ColumnVector v1(3);
        v1(1) = e_vectors(1,1);
        v1(2) = e_vectors(2,1);
        v1(3) = e_vectors(3,1);
        // consider first (smallest) eigenvector as the normal
        Real angle = (v1.t() * point_vector).AsScalar();

        // orient towards scan pose
        // works better when orientation is not flipped
        if (flipnormals && angle < 0) {
            v1 *= -1.0;
        }
        normals.push_back( Point(v1(1), v1(2), v1(3)) );
    }
}

void sobel(PointNeighbor& current_point, double& dRdTheta, double& dRdPhi)
{
  dRdPhi += 10*current_point.range_neighbors[0][1];
  dRdPhi += 3 *current_point.range_neighbors[0][0];
  dRdPhi += 3 *current_point.range_neighbors[0][2];
  dRdPhi -= 10*current_point.range_neighbors[2][1];
  dRdPhi -= 3 *current_point.range_neighbors[2][0];
  dRdPhi -= 3 *current_point.range_neighbors[2][2];
    
  dRdTheta += 10*current_point.range_neighbors[1][0];
  dRdTheta += 3 *current_point.range_neighbors[0][0];
  dRdTheta += 3 *current_point.range_neighbors[2][0];
  dRdTheta -= 10*current_point.range_neighbors[1][2];
  dRdTheta -= 3 *current_point.range_neighbors[0][2];
  dRdTheta -= 3 *current_point.range_neighbors[2][2];
}

void computeSRI(const Scan* scan, const vector<PointNeighbor>& points,
        vector<Point>& normals, bool flipnormals, int width, int height)
{
    float factorx = width/360.0;
    float factory = height/100.0;
    float deg2rad = M_PI/180;
    for(size_t pidx = 0; pidx < points.size(); ++pidx) {
        PointNeighbor current_point = points[pidx];
        int i = current_point.range_image_row;
        int j = current_point.range_image_col;
        double kart[3], polar[3];
        kart[0] = current_point.point.x;
        kart[1] = current_point.point.y;
        kart[2] = current_point.point.z;
        toPolar(kart, polar);
        double rho = polar[2];
        double phi = polar[1] * 180/M_PI;
        double theta = polar[0] * 180/M_PI;
        double n[3];

        // if no point return 0 normal
        if (rho < 10) {
            for (int i = 0; i < 3; i++)
                n[i] = 0;
            normals.push_back( Point(n[0], n[1], n[2]) );
            continue;
        }

        // partial derivative values
        double dRdTheta = 0.0;
        double dRdPhi = 0.0;

        if (i != 0 && i != height-1 && j != 0 && j != width-1)
            sobel(current_point, dRdTheta, dRdPhi);

        // add more weight to derivatives
        dRdTheta *= factory*2;
        dRdPhi *= factorx*2;

        phi += (30.0*deg2rad);

        n[0] = cos(theta) * sin(phi) - sin(theta) * dRdTheta / rho / sin(phi) + 
            cos(theta) * cos(phi) * dRdPhi / rho;

        n[1] = sin(theta) * sin(phi) + cos(theta) * dRdTheta / rho / sin(phi) + 
            sin(theta) * cos(phi) * dRdPhi / rho;

        n[2] =  cos(phi) - sin(phi) * dRdPhi / rho;

        n[2] = -n[2];

        double m = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        n[0] /= m; n[1] /= m; n[2] /= m;

        normals.push_back( Point(n[0], n[1], n[2]) );
    }
}

int main(int argc, char **argv)
{
    // commandline arguments
    int start, end;
    bool scanserver;
    int maxDist, minDist;
    string dir;
    IOType iotype;
    normal_method normalMethod;
    int knn, kmin, kmax;
    double alpha;
    int width, height;
    bool flipnormals;

    parse_options(argc, argv, start, end, scanserver, dir, iotype, maxDist,
            minDist, normalMethod, knn, kmin, kmax, alpha, width, height,
            flipnormals);

    Scan::openDirectory(scanserver, dir, iotype, start, end);

    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    boost::filesystem::path boost_dir(dir + "normals/");
    boost::filesystem::create_directory(boost_dir);

    for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;
        vector<Point> points;

        // apply optional filtering
        scan->setRangeFilter(maxDist, minDist);

        DataXYZ xyz = scan->get("xyz");
        unsigned int nPoints = xyz.size();
        for(unsigned int i = 0; i < nPoints; ++i) {
            float x, y, z;
            x = xyz[i][0];
            y = xyz[i][1];
            z = xyz[i][2];
            points.push_back(Point(x, y, z));
        }

        vector<PointNeighbor> points_neighbors;
        vector<Point> normals;

        switch (normalMethod) {
            case KNN_PCA:
                computeKNearestNeighbors(points, points_neighbors, knn);
                computePCA(scan, points_neighbors, normals, flipnormals);
                break;
            case AKNN_PCA:
                computeKNearestNeighbors(points, points_neighbors, kmin, kmax, alpha);
                computePCA(scan, points_neighbors, normals, flipnormals);
                break;
            case PANO_PCA:
                computePanoramaNeighbors(scan, points_neighbors, width, height);
                computePCA(scan, points_neighbors, normals, flipnormals);
                break;
            case PANO_SRI:
                computePanoramaNeighbors(scan, points_neighbors, width, height);
                computeSRI(scan, points_neighbors, normals, flipnormals, width, height);
                break;
            default:
                break;
        }

        if (points.size() != normals.size()) {
            cerr << "got " << points.size() << " points but " << normals.size() << " normals" << endl;
            return 1;
        }

        writeNormals(scan, dir + "normals/", points, normals);
    }

    Scan::closeDirectory();

    return 0;
}
