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

enum search_method {KNN, RANGE};

void search_option_dependency(const po::variables_map & vm, search_method ntype, const char *option)
{
    if (vm.count("search") && vm["search"].as<search_method>() == ntype) {
        if (!vm.count(option)) {
            throw std::logic_error (string("this search method needs ")+option+" to be set");
        }
    }
}

void search_option_conflict(const po::variables_map & vm, search_method ntype, const char *option)
{
    if (vm.count("search") && vm["search"].as<search_method>() == ntype) {
        if (vm.count(option)) {
            throw std::logic_error (string("this search method is incompatible with ")+option);
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
        search_method*, int) {
    if (values.size() == 0)
        throw std::runtime_error("Invalid model specification");
    string arg = values.at(0);
    if(strcasecmp(arg.c_str(), "KNN") == 0) v = KNN;
    else if(strcasecmp(arg.c_str(), "RANGE") == 0) v = RANGE;
    else throw std::runtime_error(std::string("search method ") + arg + std::string(" is unknown"));
}

/*
 * parse commandline options, fill arguments
 */
void parse_options(int argc, char **argv, int &start, int &end,
        bool &scanserver, string &dir, IOType &iotype,
        int &maxDist, int &minDist, search_method &searchMethod, 
        int &knn, double &range)
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
        ("search,a", po::value<search_method>(&searchMethod)->default_value(KNN),
         "choose the search method:\n"
         "KNN -- use kNN search\n"
         "RANGE -- use range search\n")
        ("knn,K", po::value<int>(&knn),
         "select the k in kNN search")
        ("range,R", po::value<double>(&range),
         "select the range in range search")
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

    search_option_dependency(vm, KNN, "knn");
    search_option_conflict(vm, KNN, "range");

    search_option_dependency(vm, RANGE, "range");
    search_option_conflict(vm, RANGE, "knn");

    if (vm["start"].as<int>() > vm["end"].as<int>()) {
        throw std::runtime_error ("--end must be bigger or equal --start");
    }

    // add trailing slash to directory if not present yet
    if (dir[dir.length()-1] != '/') dir = dir + "/";
}

void calculateNeighborsANN(vector<Point> &points, int k, vector<vector<int>> &neighbors) {
    cout<<"Total number of points: "<<points.size()<<endl;
    ANNpointArray pa = annAllocPts(points.size(), 3);
    for (size_t i=0; i<points.size(); ++i) {
        pa[i][0] = points[i].x;
        pa[i][1] = points[i].y;
        pa[i][2] = points[i].z;
    }

    ANNkd_tree t(pa, points.size(), 3);
    ANNidxArray nidx = new ANNidx[k];
    ANNdistArray d = new ANNdist[k];

    neighbors.reserve(points.size());

    for (size_t i=0; i<points.size(); ++i) {
        ANNpoint p = pa[i];
        t.annkSearch(p, k, nidx, d, 0.0);
        vector<int> n;
        for (int j=0; j<k; ++j) {
            n.push_back(nidx[j]);
        }
        neighbors.push_back(n);
    }

    annDeallocPts(pa);
}

int main(int argc, char **argv)
{
    // commandline arguments
    int start, end;
    bool scanserver;
    int maxDist, minDist;
    string dir;
    IOType iotype;
    search_method searchMethod;
    int knn;
    double range;

    parse_options(argc, argv, start, end, scanserver, dir, iotype, maxDist,
            minDist, searchMethod, knn, range);

    cout << "Using: " << searchMethod << " knn: " << knn << " range: " << range << endl;

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
            vector<Point> points;
            vector<vector<int>> neighbors;

            points.reserve(xyz.size());

            for(size_t j = 0; j < xyz.size(); j++) {
                points.push_back(Point(xyz[j][0], xyz[j][1], xyz[j][2]));
            }

            calculateNeighborsANN(points, 5, neighbors);

            if (points.size() != neighbors.size()) {
                cerr << "unequal number of points and neighbors" << endl;
                exit(-1);
            }

            for (size_t j = 0; j < points.size(); ++j) {
                cout << j << ":";
                for (size_t m = 0; m < neighbors[j].size(); ++m) {
                    cout << " " << neighbors[j][m];
                }
                cout << endl;
            }
        }
        Scan::closeDirectory();
    }

    return 0;
}
