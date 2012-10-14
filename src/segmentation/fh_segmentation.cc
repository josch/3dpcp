#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unistd.h>

using namespace std;

#include "slam6d/globals.icc"
#include "scanio/scan_io.h"

//#include "segmentation/Options.h"
#include "segmentation/segment-image.h"
#include "segmentation/image.h"
#include "segmentation/pnmfile.h"
#include "segmentation/misc.h"
#include "segmentation/FHGraph.h"

#include "scanserver/clientInterface.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct Options
{
    IOType type;
    float minDist;
    float maxDist;
    int start;
    int end;

    float sigma;
    float k;
    int min_size;

    float eps;

    std::string dir;
    std::string outdir;

    int reserve;

    int neighbors;
    float radius;

    bool scanserver;
};

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

struct Options parse_options(int argc, char* argv[])
{
    struct Options options;

    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "output this help message");

    po::options_description input("Input options");
    input.add_options()
        ("start,s", po::value<int>(&options.start)->default_value(0),
         "start at scan <arg> (i.e., neglects the first <arg> scans) "
         "[ATTENTION: counting naturally starts with 0]")
        ("end,e", po::value<int>(&options.end)->default_value(-1),
         "end after scan <arg>")
        ("format,f", po::value<IOType>(&options.type)->default_value(UOS),
         "using shared library <arg> for input. (chose F from {uos, uos_map, "
         "uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, "
         "riegl_txt, riegl_rgb, riegl_bin, zahn, ply})")
        ("max,M", po::value<float>(&options.maxDist)->default_value(-1),
         "neglegt all data points with a distance larger than <arg> 'units")
        ("min,m", po::value<float>(&options.minDist)->default_value(-1),
         "neglegt all data points with a distance smaller than <arg> 'units")
        ("scanserver,S", po::bool_switch(&options.scanserver),
         "Use the scanserver as an input method and handling of scan data");

    po::options_description output("Output options");
    output.add_options()
        ("output,o", po::value<string>(&options.outdir)->default_value("segments"),
         "output directory");

    po::options_description segment("Segmentation options");
    segment.add_options()
        ("sigma,g", po::value<float>(&options.sigma)->default_value(1),
         "the sigma value used in Gaussian smoothing")
        (",k", po::value<float>(&options.k)->default_value(1),
         "the K value used in the FH segmentation")
        ("min-size,I", po::value<int>(&options.min_size)->default_value(0),
         "the min size of a segment")
        ("radius,r", po::value<float>(&options.radius)->default_value(-1),
         "use radius search, specify the radius, if RADIUS > 0")
        ("nearest,n", po::value<int>(&options.neighbors)->default_value(-1),
         "use approximate NR-nearest neighbors search or limit the number of "
         "points returned by the radius search (it will pick NR points "
         "randomly)")
        ("size,R", po::value<int>(&options.reserve)->default_value(-1),
         "the size of pre-reserved std::vectors")
        ("epsilon,A", po::value<float>(&options.eps)->default_value(1.0),
         "the error used by the AKNN algorithm");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-dir", po::value<string>(&options.dir), "input dir");

    // all options
    po::options_description all;
    all.add(generic).add(input).add(output).add(segment).add(hidden);

    // options visible with --help
    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input).add(output).add(segment);

    // positional argument
    po::positional_options_description pd;
    pd.add("input-dir", 1);

    // process options
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), vm);
    po::notify(vm);

    // display help
    if (vm.count("help")) {
        cout << cmdline_options;
        exit(0);
    }

    if ( options.dir[options.dir.length()-1] != '/' ) options.dir += "/";
    if ( options.outdir[options.outdir.length()-1] != '/' ) options.outdir += "/";

    boost::filesystem::path boost_dir(options.outdir);
    if (!boost::filesystem::create_directory(boost_dir)) {
        cerr << "Couldn't create directory " << options.outdir << endl;
    }

    cout << endl << "Input directory is: " << options.dir << endl;
    cout << "Output directory is: " << options.outdir << endl;

    if ( options.neighbors<0 && options.radius<0 )
    {
        cerr << endl << "Underspecification! Nearest neighbors is: " << options.neighbors << " and search radius is: " << options.radius << endl;
        cout << cmdline_options;
        exit(1);
    }

    return options;
}

struct mycomp
{
    int rgb2int(rgb x) { return (x.r << 16) + (x.g << 8) + x.b; };
    bool operator()(rgb a, rgb b)
    {
        return rgb2int(a) < rgb2int(b);
    }
};

double weight1(Point a, Point b)
{
    return a.distance(b);
}

double weight2(Point a, Point b)
{
    return a.distance(b) * .5 + fabs(a.reflectance-b.reflectance) * .5;
}

bool compareSizes(vector<Point>* a, vector<Point>* b)
{
    return a->size() > b->size();
}

void writeScan(const string& dir, 
        unsigned int scan_number,
        const vector<Point>& points) {

    stringstream ss; 
    ss << dir << "scan" << std::setw(3) << std::setfill('0') << scan_number << ".3d"; 
    ofstream scan_file;
    scan_file.open(ss.str().c_str());
    for(size_t i = 0;  i < points.size(); ++i) {
        scan_file << points[i].x << " " << points[i].y << " " << points[i].z << "\n";  
    }
    scan_file.close();

    ss.clear(); ss.str(string());
    ss << dir << "scan" << std::setw(3) << std::setfill('0') << scan_number << ".pose";
    ofstream pose_file; 
    pose_file.open(ss.str().c_str());
    pose_file << 0 << " " << 0 << " " << 0 << "\n" << 0 << " " << 0 << " " << 0 << "\n";
    pose_file.close();
}

int main(int argc, char* argv[])
{
    struct Options options = parse_options(argc, argv);

    if (options.scanserver) {
        try {
            ClientInterface::create();
        } catch(std::runtime_error& e) {
            cerr << "ClientInterface could not be created: " << e.what() << endl;
            cerr << "Start the scanserver first." << endl;
            exit(-1);
        }
    }

    /* Read the points */
    vector<Point> points;
    if ( options.reserve > 0 )
        points.reserve(options.reserve);

    Scan::openDirectory(options.scanserver, options.dir, options.type, options.start, options.end);

    if (Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    for (std::vector<Scan*>::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;
        DataXYZ xyz = scan->get("xyz");
        DataReflectance xyz_reflectance = scan->get("reflectance");
        unsigned int nPoints = xyz.size();

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

            points.push_back(Point(x, y, z));
        }

        if (options.scanserver) {
            scan->clear(DATA_XYZ | DATA_REFLECTANCE);
        }
    }

    if (options.scanserver) {  
        Scan::closeDirectory();
        ClientInterface::destroy();
    }

    cerr << endl << "Loaded " << points.size() << " points" << endl;

    FHGraph g(points, weight2, options.sigma, options.eps, options.neighbors, options.radius);

    edge* e = g.getGraph();
    universe* segmented = segment_graph(g.getNumPoints(), g.getNumEdges(), e, options.k);

    for (int i=0; i<g.getNumEdges(); ++i)
    {
        int a = e[i].a;
        int b = e[i].b;

        int aa = segmented->find(a);
        int bb = segmented->find(b);

        if ( (aa!=bb) &&
                (segmented->size(aa) < options.min_size ||
                 segmented->size(bb) < options.min_size) )
            segmented->join(aa, bb);
    }
    delete[] e;

    int nr = segmented->num_sets();
    cerr << "Obtained " << nr << " segment(s)" << endl;

    vector< vector<Point>* > clouds;
    clouds.reserve(nr);
    for (int i=0; i<nr; ++i)
        clouds.push_back( new vector<Point> );
    map<int, int> c2c; // component to cloud
    int k = 0;
    for (int i=0; i<g.getNumPoints(); ++i)
    {
        int comp = segmented->find(i);
        if ( c2c.find(comp)==c2c.end() )
        {
            c2c[comp] = k++;
            clouds[c2c[comp]]->reserve(segmented->size(comp));
        }
        clouds[c2c[comp]]->push_back(g[i]);
    }
    g.dispose();

    sort(clouds.begin(), clouds.end(), compareSizes);

#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<nr; ++i)
    {
        writeScan(options.outdir, i, * (clouds[i]));
    }

    return 0;
}
