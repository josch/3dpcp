#include <iostream>
#include <vector>
#include <map>
#include <unistd.h>
using namespace std;

#include "slam6d/scan_io_rxp.h"
#include "slam6d/scan_io_xyzr.h"
#include "segmentation/Options.h"
#include "segmentation/segment-image.h"
#include "segmentation/image.h"
#include "segmentation/pnmfile.h"
#include "segmentation/misc.h"
#include "segmentation/FHGraph.h"
#include <segmentation/SRI.h>
#include "slam6d/globals.icc"
#include <slam6d/scan_io_velodyne.h>
#include <segmentation/Timer.h>

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

int main(int argc, char* argv[])
{
	Options options(argc, argv);

	if ( options.neighbors<0 && options.radius<0 )
	{
		options.usage();
		return 1;
	}


	/* Read the points */
	vector<Point> points;
	if ( options.reserve > 0 )
		points.reserve(options.reserve);
	double euler[] = {0,0,0,0,0,0};
	ScanIO* scan;
	switch ( options.type )
	{
		case RXP:
			scan = new ScanIO_rxp();
			break;
		case VELODYNE:
			scan = new ScanIO_velodyne();
			break;
		default:
			cerr << "Type not supported (yet)!" << endl;
			return 1;
	}

	Timer *t = new Timer("Reading points... # ms");
	scan->readScans(options.start, options.end, options.dir,
		options.maxDist, options.minDist, euler, points
	);
	delete scan;
	delete t;
	cerr << "Loaded " << points.size() << " points" << endl;

	t = new Timer("Creating graph... # ms");
	FHGraph g(points, weight2, options.sigma, options.eps, options.neighbors, options.radius);
	delete t;

	t = new Timer("Segmenting image... # ms");
	edge* e = g.getGraph();
	universe* segmented = segment_graph(g.getNumPoints(), g.getNumEdges(), e, options.k);
	delete t;

	t = new Timer("Post processing... # ms");
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
	delete t;

	delete[] e;

	int nr = segmented->num_sets();
	cerr << "Obtained " << nr << " segment(s)" << endl;

	t = new Timer("Creating point clouds... # ms");
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
	delete t;

	g.dispose();

	t = new Timer("Saving PPM image... # ms");
	image<rgb> img(3600, 1800);
	for (int i=0; i<nr; ++i)
	{
		rgb c = random_rgb();
		vector<Point>* v = clouds[i];
		for (size_t j=0; j<v->size(); ++j)
		{
			SRI::point p((*v)[j].x, (*v)[j].y, (*v)[j].z);
			img.access[(int)p.y][(int)p.x] = c;
		}
	}
	savePPM(&img, "output.ppm");
	delete t;
	
	if (options.do_out)
	{
		t = new Timer("Sorting point clouds... # ms");
		sort(clouds.begin(), clouds.end(), compareSizes);
		delete t;
		
		t = new Timer("Dumping point clouds... # ms");
		#pragma omp parallel for schedule(dynamic)
		for (int i=0; i<nr; ++i)
		{
			ScanIO_xyzr::writeScan(options.outdir, i, * (clouds[i]));
		}
		delete t;
	}

	return 0;
}
