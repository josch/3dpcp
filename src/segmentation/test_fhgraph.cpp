#include <segment/FHGraph.h>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;

double f(Point p1, Point p2)
{
	return p1.distance(p2);
}

int main()
{
	vector<Point> p;
	p.push_back(Point(1,0,0));
	p.push_back(Point(2,0,0));
	p.push_back(Point(3,1,0));
	p.push_back(Point(1.5,2,0));
	p.push_back(Point(0,1,0));

	double w[5][5];
	w[0][0] = w[1][1] = w[2][2] = w[3][3] = w[4][4] = 0;
	w[0][1] = w[1][0] = 1.74;
	w[0][2] = w[2][0] = 1.87;
	w[0][3] = w[3][0] = 1.79;
	w[0][4] = w[4][0] = 1.92;
	w[1][2] = w[2][1] = 1.92;
	w[1][3] = w[3][1] = 1.79;
	w[1][4] = w[4][1] = 1.87;
	w[2][3] = w[3][2] = 2.03;
	w[2][4] = w[4][2] = 2.04;
	w[3][4] = w[4][3] = 2.03;
	
	FHGraph g(p, f, 0, 0, 4, -1);

	edge* e = g.getGraph();
	for (int i=0; i<g.getNumEdges(); ++i)
	{
		double err = fabs(w[e[i].a][e[i].b]-e[i].w)/w[e[i].a][e[i].b];
		cerr << e[i].a << " " << e[i].b << " " << e[i].w << " " << err << endl;
		if (  err > 0.1 )
		{
			return 1;
		}
	}
	delete[] e;
	return 0;
}
