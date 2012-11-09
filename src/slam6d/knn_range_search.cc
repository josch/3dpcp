/*
 * knn_range_search implementation
 *
 * Copyright (C) Johannes Schauer, Razvan Mihalyi
 *
 * Released under the GPL version 3
 *
 */

#include "ANN/ANN.h"

#include "slam6d/point.h"
#include "slam6d/scan.h"
#include "slam6d/globals.icc"
#include "slam6d/kd.h"

#include <string>
using std::string;

#include <iostream>
#include <fstream>
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
        int &maxDist, int &minDist, int &knn, double &range, size_t &maxpoints)
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

    po::options_description search("Search options");
    search.add_options()
        ("knn,K", po::value<int>(&knn)->default_value(0),
         "if greater than 0 do knn search, otherwise range search")
        ("range,R", po::value<double>(&range)->default_value(20.0),
         "select the max range for knn and range search")
        ("maxpoints,p", po::value<size_t>(&maxpoints),
         "maximum number of points to investigate")
        ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(search).add(hidden);

    // options visible with --help
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input).add(search);

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

    if (vm["start"].as<int>() > vm["end"].as<int>()) {
        throw std::runtime_error ("--end must be bigger or equal --start");
    }

    // add trailing slash to directory if not present yet
    if (dir[dir.length()-1] != '/') dir = dir + "/";
}

void calculateANN(double **points, size_t nPoints, int k, double range, vector<vector<double *>> &neighbors) {
    ANNpointArray pa = annAllocPts(nPoints, 3);
    for (size_t i=0; i<nPoints; ++i) {
        pa[i][0] = points[i][0];
        pa[i][1] = points[i][1];
        pa[i][2] = points[i][2];
    }

    ANNkd_tree t(pa, nPoints, 3);
    neighbors.reserve(nPoints);

    double sqradius = sqr(range);

    for (size_t i=0; i<nPoints; ++i) {
        ANNpoint p = pa[i];
        int m = t.annkFRSearch(p, sqradius, 0, NULL, NULL, 0.0);

        // make sure enough space is allocated but at least k
        m = m < k ? k : m;
        ANNidxArray nidx = new ANNidx[m];
        ANNdistArray d = new ANNdist[m];

        // if k is zero, do range search and set k to m
        if (k == 0) k = m;

        t.annkFRSearch(p, sqradius, k, nidx, d, 0.0);
        vector<double *> n;
        for (int j=0; j<k; ++j) {
            if (nidx[j] == ANN_NULL_IDX) {
                break; // ANN_NULL_IDX marks the end
            }
            if (nidx[j] >= (int) nPoints) {
                cerr << endl << "nidx[j] >= nPoints" << endl;
                continue;
            }
            n.push_back(points[nidx[j]]);
        }
        neighbors.push_back(n);
    }

    annDeallocPts(pa);
}

void calculateKdTree(double **points, size_t nPoints, int k, double range, vector<vector<double *>> &neighbors) {
    /// KDtree range search
    KDtree kd_tree(points, nPoints);

    neighbors.reserve(nPoints);

    for (size_t i=0; i<nPoints; ++i) {
        vector<double *> n;
        kd_tree.FindClosestKNNRange(points[i], sqr(range), n, k);
        // check distances of found neighbors
        for (size_t j = 0; j < n.size(); ++j) {
            if (sqrt(Dist2(points[i], n[j])) > range) {
                cerr << endl << "neighbor distance greater than radius" << endl;
            }
        }
        neighbors.push_back(n);
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
    int knn;
    double range;
    size_t maxpoints;

    parse_options(argc, argv, start, end, scanserver, dir, iotype, maxDist,
            minDist, knn, range, maxpoints);

    for (int iter = start; iter <= end; iter++) {
        Scan::openDirectory(scanserver, dir, iotype, iter, iter);
        if(Scan::allScans.size() == 0) {
            cerr << "No scans found. Did you use the correct format?" << endl;
            exit(-1);
        }
        for(ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
            Scan* scan = *it;
            scan->setRangeFilter(maxDist, minDist);

            DataXYZ xyz(scan->get("xyz"));
            vector<vector<double *>> neighborsANN;
            vector<vector<double *>> neighborsKD;

            size_t maxp = maxpoints > xyz.size() ? xyz.size() : maxpoints;

            cerr << "considering " << maxp << " datapoints" << endl;

            double **points = new double*[maxp];
            for (unsigned int i = 0; i < maxp; ++i) {
                points[i] = new double[3];
                for (unsigned int j = 0; j < 3; ++j) 
                    points[i][j] = xyz[i][j];
            }

            calculateANN(points, maxp, knn, range, neighborsANN);
            calculateKdTree(points, maxp, knn, range, neighborsKD);

            double epsilon = 1e-5;

            ofstream fout("output");

            for (size_t j = 0; j < maxp; ++j) {
                bool fail = false;
                cout << j;
                fout << "Point " << j << ": " << points[j][0] << " " << points[j][1] << " " << points[j][2] << endl;
                fout << "--------------------------------------------------" << endl;
                fout << "ANN neighbors: " << endl;
                for (size_t m = 0; m < neighborsANN[j].size(); ++m) {
                    fout << neighborsANN[j][m][0] << " " << neighborsANN[j][m][1] << " " << neighborsANN[j][m][2] << "\tDist: " << sqrt(Dist2(points[j], neighborsANN[j][m])) << endl;
                    bool found = false;
                    for (size_t n = 0; n < neighborsKD[j].size(); ++n) {
                        if (fabs(neighborsANN[j][m] - neighborsKD[j][n]) < epsilon) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        // compute distance between point and neighbor
                        double d = sqrt(Dist2(points[j], neighborsANN[j][m]));
                        cout << " (ann not kd: " << d << ")";
                        fail = true;
                    }
                }
                fout << "KD neighbors: " << endl;
                for (size_t m = 0; m < neighborsKD[j].size(); ++m) {
                    fout << neighborsKD[j][m][0] << " " << neighborsKD[j][m][1] << " " << neighborsKD[j][m][2] << "\tDist: " << sqrt(Dist2(points[j], neighborsKD[j][m])) << endl;
                    bool found = false;
                    for (size_t n = 0; n < neighborsANN[j].size(); ++n) {
                        if (fabs(neighborsANN[j][n] - neighborsKD[j][m]) < epsilon) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        // compute distance between point and neighbor
                        double d = sqrt(Dist2(points[j], neighborsKD[j][m]));
                        cout << " (kd not ann: " << d << ")";
                        fail = true;
                    }
                }
                // if there is no fail yet, compare the neighbor vector sizes
                // to make sure that the lists actually contain the same and
                // there are no duplicates
                if (!fail) {
                    if (neighborsANN[j].size() != neighborsKD[j].size()) {
                        cout << " " << neighborsANN[j].size() << " ANN neighbors ";
                        cout << "but " << neighborsKD[j].size() << " kd neighbors";
                        fail = true;
                    }
                }
                if (!fail)
                    cout << " success!";
                cout << endl;
                fout << endl << endl << endl;
            }
            fout.flush();
            fout.close();
            for (unsigned int i = 0; i < maxp; ++i) {
                delete []points[i];
            }
            delete []points;
        }
        Scan::closeDirectory();
    }

    return 0;
}
