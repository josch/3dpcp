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
using std::map;
using std::pair;

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
        ("maxpoints,p", po::value<size_t>(&maxpoints)->default_value(0),
         "maximum number of points to investigate - 0 means all points")
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

void calculateANN(double **points, size_t nPoints, int knn, double range, map<double*, vector<double *>> &neighbors) {
    int k;

    ANNpointArray pa = annAllocPts(nPoints, 3);
    for (size_t i=0; i<nPoints; ++i) {
        pa[i][0] = points[i][0];
        pa[i][1] = points[i][1];
        pa[i][2] = points[i][2];
    }

    ANNkd_tree t(pa, nPoints, 3);

    double sqradius = sqr(range);

    for (size_t i=0; i<nPoints; ++i) {
        ANNpoint p = pa[i];
        int m = t.annkFRSearch(p, sqradius, 0, NULL, NULL, 0.0);

        // make sure enough space is allocated but at least knn
        m = m < knn ? knn : m;
        ANNidxArray nidx = new ANNidx[m];
        ANNdistArray d = new ANNdist[m];

        // if k is zero, do range search and set k to m
        if (knn == 0)
            k = m;
        else
            k = knn;

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
        neighbors.insert(pair<double*, vector<double *>>(points[i], n));
    }

    annDeallocPts(pa);
}

void calculateKdTree(double **points, size_t nPoints, int k, double range, map<double*, vector<double *>> &neighbors) {
    /// KDtree search
    KDtree kd_tree(points, nPoints);

    for (size_t i=0; i<nPoints; ++i) {
        vector<double *> n;
        kd_tree.FindClosestKNNRange(points[i], sqr(range), n, k);
        // check distances of found neighbors
        for (size_t j = 0; j < n.size(); ++j) {
            if (sqrt(Dist2(points[i], n[j])) > range) {
                cerr << endl << "neighbor distance greater than radius" << endl;
            }
        }
        neighbors.insert(pair<double*, vector<double *>>(points[i], n));
    }
}

struct ndist {
    double *p;
    bool operator() (double *a, double *b) { return Dist2(a,p)<Dist2(b,p); }
} ndist;

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
    long starttime, endtime;

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

            // we need a map datatype because KDTree changes the order of its
            // input. Since we therefor cannot rely on the order we need to
            // save a direct mapping instead.
            map<double*, vector<double *>> neighborsANN;
            map<double*, vector<double *>> neighborsKD;

            size_t maxp;
            if (maxpoints == 0)
                maxp = xyz.size();
            else
                if (maxpoints > xyz.size())
                    maxp = xyz.size();
                else
                    maxp = maxpoints;

            cerr << "considering " << maxp << " datapoints" << endl;

            double **points = new double*[maxp];
            for (unsigned int i = 0; i < maxp; ++i) {
                points[i] = new double[3];
                for (unsigned int j = 0; j < 3; ++j) 
                    points[i][j] = xyz[i][j];
            }

            if (knn < 1)
                cout << "using range search (within a range of " << range << " cms)" << endl;
            else
                cout << "using knn search (within range of " << range << " cms)" << endl;

            starttime = GetCurrentTimeInMilliSec();
            calculateANN(points, maxp, knn, range, neighborsANN);
            endtime = GetCurrentTimeInMilliSec() - starttime;
            cout << "calculateANN done in " << endtime << " milliseconds!!!" << endl;

            starttime = GetCurrentTimeInMilliSec();
            calculateKdTree(points, maxp, knn, range, neighborsKD);
            endtime = GetCurrentTimeInMilliSec() - starttime;
            cout << "calculateKdTree done in " << endtime << " milliseconds!!!" << endl;

            cout << "writing debugging output to ./knn_range_search.log" << endl;
            ofstream fout("knn_range_search.log");

            double epsilon = 1.0e-15;

            cout.precision(15);

            // for each point, compare the neighbors calculated by ANN and
            // kdtree
            for (size_t j = 0; j < maxp; ++j) {
                bool fail = false;

                cout << j;
                fout << "Point " << j << ": " << points[j][0] << " " << points[j][1] << " " << points[j][2] << endl;
                fout << "--------------------------------------------------" << endl;
                fout << "ANN neighbors: " << endl;

                vector<double *> nANN = neighborsANN.find(points[j])->second;
                vector<double *> nKD = neighborsKD.find(points[j])->second;

                for (size_t m = 0; m < nANN.size(); ++m) {
                    fout << nANN[m][0] << " " << nANN[m][1] << " " << nANN[m][2] << "\tDist: " << sqrt(Dist2(points[j], nANN[m])) << endl;

                    bool found = false;
                    // try to find this ANN neighbor in the kdtree neighbors
                    for (size_t n = 0; n < nKD.size(); ++n) {
                        // just compare the pointer values
                        if (nANN[m] == nKD[n]) {
                            found = true;
                            break;
                        }
                    }

                    // if we do knn search, then it can happen that the
                    // farthest points in each nearest neighbor list are of
                    // equal distance but actually have different coordinates
                    // so if a correspondence was not found at this point and
                    // knn search is used, check for that
                    if (!found && knn) {
                        // get the greatest distances in both neighbor lists
                        // sort, as it is not guaranteed that both lists are
                        // sorted
                        ndist.p = points[j];
                        double *farthest_ann = *max_element(nANN.begin(), nANN.end(), ndist);
                        double *farthest_kd = *max_element(nKD.begin(), nKD.end(), ndist);

                        // check if the current neighbor is within an epsilon
                        // distance range from the greatest distance in the
                        // neighbor list
                        if (fabs(sqrt(Dist2(points[j], nANN[m]))-sqrt(Dist2(points[j], farthest_ann))) < epsilon) {
                            // if it is, check whether it is equally distant
                            // from the query point as the farthest point in
                            // the other neighbor list
                            if (fabs(sqrt(Dist2(points[j], nANN[m]))-sqrt(Dist2(points[j], farthest_kd))) < epsilon) {
                                // if such a point is found, than it was
                                // only not found before because it
                                // happened to be cut-off as it was of
                                // equal distance
                                found = true;
                            }
                        }
                    }

                    // still not found - error
                    if (!found) {
                        // compute and print distance between point and neighbor
                        double d = sqrt(Dist2(points[j], nANN[m]));
                        cout << " (ann not kd: " << d << ")";
                        fail = true;
                    }
                }
                fout << "KD neighbors: " << endl;
                for (size_t m = 0; m < nKD.size(); ++m) {
                    fout << nKD[m][0] << " " << nKD[m][1] << " " << nKD[m][2] << "\tDist: " << sqrt(Dist2(points[j], nKD[m])) << endl;

                    bool found = false;
                    // try to find this kdtree neighbor in the ANN neighbors
                    for (size_t n = 0; n < nANN.size(); ++n) {
                        // just compare the pointer values
                        if (nANN[n] == nKD[m]) {
                            found = true;
                            break;
                        }
                    }
                    // if we do knn search, then it can happen that the
                    // farthest points in each nearest neighbor list are of
                    // equal distance but actually have different coordinates
                    // so if a correspondence was not found at this point and
                    // knn search is used, check for that
                    if (!found && knn) {
                        // get the greatest distances in both neighbor lists
                        // sort, as it is not guaranteed that both lists are
                        // sorted
                        ndist.p = points[j];
                        double *farthest_ann = *max_element(nANN.begin(), nANN.end(), ndist);
                        double *farthest_kd = *max_element(nKD.begin(), nKD.end(), ndist);

                        // check if the current neighbor is within an epsilon
                        // distance range from the greatest distance in the
                        // neighbor list
                        if (fabs(sqrt(Dist2(points[j], nKD[m]))-sqrt(Dist2(points[j], farthest_kd))) < epsilon) {
                            // if it is, check whether it is equally distant
                            // from the query point as the farthest point in
                            // the other neighbor list
                            if (fabs(sqrt(Dist2(points[j], nKD[m]))-sqrt(Dist2(points[j], farthest_ann))) < epsilon) {
                                // if such a point is found, than it was
                                // only not found before because it
                                // happened to be cut-off as it was of
                                // equal distance
                                found = true;
                            }
                        }
                    }
                    // still not found - error
                    if (!found) {
                        // compute and print distance between point and neighbor
                        double d = sqrt(Dist2(points[j], nKD[m]));
                        cout << " (kd not ann: " << d << ")";
                        fail = true;
                    }
                }
                // if there is no fail yet, compare the neighbor vector sizes
                // to make sure that the lists actually contain the same and
                // there are no duplicates
                if (!fail) {
                    if (nANN.size() != nKD.size()) {
                        cout << " " << nANN.size() << " ANN neighbors ";
                        cout << "but " << nKD.size() << " kd neighbors";
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
