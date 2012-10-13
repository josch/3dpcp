#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <unistd.h>

using namespace std;

#include "slam6d/globals.icc"
#include "scanio/scan_io.h"

#include "segmentation/Options.h"
#include "segmentation/segment-image.h"
#include "segmentation/image.h"
#include "segmentation/pnmfile.h"
#include "segmentation/misc.h"
#include "segmentation/FHGraph.h"

#include "scanserver/clientInterface.h"

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
	Options options(argc, argv);

	if ( options.neighbors<0 && options.radius<0 )
	{
    cerr << endl << "Underspecification! Nearest neighbors is: " << options.neighbors << " and search radius is: " << options.radius << endl;
		options.usage();
		return 1;
	}

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

	if (options.do_out)
	{
	  sort(clouds.begin(), clouds.end(), compareSizes);

		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<nr; ++i)
		{
			writeScan(options.outdir, i, * (clouds[i]));
		}
	}

	return 0;
}
