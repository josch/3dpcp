/**
  * Point Cloud Segmentation using Felzenszwalb-Huttenlocher Algorithm
  *
  * Copyright (C) Jacobs University Bremen
  *
  * Released under the GPL version 3.
  *
  * @author Mihai-Cotizo Sima
  */

#include <segmentation/FHGraph.h>
#include <map>
#include <omp.h>
#include <algorithm>

using namespace std;



FHGraph::FHGraph(std::vector< Point >& ps, double weight(Point, Point), double sigma, double eps, int neighbors, float radius) :
    points( ps ), V( ps.size() )
{
    /*
     * 1. create adjency list using a map<int, vector<half_edge> >
     * 2. use get_neighbors(e, max_dist) to get all the edges e' that are at a distance smaller than max_dist than e
     * 3. using all these edges, compute the gaussian smoothed weight
     * 4. insert the edges in a new list
     */
    nr_neighbors = neighbors;
    this->radius = radius;

    /// convert from vector<Point> to double** (required by KDtree)
    this->points_size = points.size();
    points_ptr = new double*[points_size];
    for (size_t i = 0; i < points_size; ++i) {
        points_ptr[i] = new double[3];
        points_ptr[i][0] = points[i].x;
        points_ptr[i][1] = points[i].y;
        points_ptr[i][2] = points[i].z;
    }

    compute_neighbors_kd(weight, eps);    
    //compute_neighbors_ann(weight, eps);

    if ( sigma > 0.01 )
    {
        do_gauss(sigma);
    }
    else
    {
        without_gauss();
    }

    adjency_list.clear();
}

FHGraph::~FHGraph() 
{
    for (size_t i = 0; i < points_size; ++i)
        delete[] points_ptr[i];
    delete[] points_ptr;
}

void FHGraph::compute_neighbors_kd(double weight(Point, Point), double eps)
{
    adjency_list.reserve(points_size);
    adjency_list.resize(points_size);

    /// KDtree search
    KDtree kd_tree(points_ptr, points_size);

    for (size_t i = 0; i < points_size; ++i) {
        vector<double *> neighbors;
        /// distinction between knn and range search done within kd_tree
        kd_tree.FindClosestKNNRange(points_ptr[i], sqr(radius), neighbors, nr_neighbors);
        // check distances of found neighbors
        for (size_t j = 0; j < neighbors.size(); ++j) {
            if (Dist2(points_ptr[i], neighbors[j]) > sqr(radius))
                cerr << "neighbor distance greater than radius" << endl;

            he e;
            /// need to recover the index (in the original point array) of the current neighbor
            Point neighbor(neighbors[j][0], neighbors[j][1], neighbors[j][2]);
            e.x = std::find(points.begin(), points.end(), neighbor) - points.begin();
            if (e.x > (int) points.size()) {
                cerr << "couldn't find neighbor in initial points" << endl;
                return;
            }
            e.w = weight(points_ptr[i], neighbors[j]);
            adjency_list[i].push_back(e);
        }
    }
}

void FHGraph::compute_neighbors_ann(double weight(Point, Point), double eps)
{

    adjency_list.reserve(points.size());
    adjency_list.resize(points.size());

    ANNpointArray pa = annAllocPts(points.size(), 3);
    for (size_t i=0; i<points.size(); ++i)
    {
        pa[i] = new ANNcoord[3];
        pa[i][0] = points[i].x;
        pa[i][1] = points[i].y;
        pa[i][2] = points[i].z;
    }

    ANNkd_tree t(pa, points.size(), 3);

    if ( radius < 0 ) // Using knn search
    {
        nr_neighbors++;
        ANNidxArray n = new ANNidx[nr_neighbors];
        ANNdistArray d = new ANNdist[nr_neighbors];

        for (size_t i=0; i<points.size(); ++i)
        {
            ANNpoint p = pa[i];

            t.annkSearch(p, nr_neighbors, n, d, eps);

            for (int j=0; j<nr_neighbors; ++j)
            {
                if ( n[j] == (int)i ) continue;

                he e;
                e.x = n[j];
                e.w = weight(points[i], points[n[j]]);

                adjency_list[i].push_back(e);
            }

        }

        delete[] n;
        delete[] d;
    }
    else // Using radius search
    {
        float sqradius = radius*radius;


        ANNidxArray n;
        ANNdistArray d;
        int nret;
        int total = 0;

        const int MOD = 1000;
        int TMP = MOD;

        for (size_t i=0; i<points.size(); ++i)
        {
            ANNpoint p = pa[i];

            nret = t.annkFRSearch(p, sqradius, 0, NULL, NULL, eps);
            total += nret;

            n = new ANNidx[nret];
            d = new ANNdist[nret];
            t.annkFRSearch(p, sqradius, nret, n, d, eps);

            if ( nr_neighbors > 0 && nr_neighbors < nret )
            {
                random_shuffle(n, n+nret);
                nret = nr_neighbors;
            }

            for (int j=0; j<nret; ++j)
            {
                if ( n[j] == (int)i ) continue;

                he e;
                e.x = n[j];
                e.w = weight(points[i], points[n[j]]);

                adjency_list[i].push_back(e);
            }

            delete[] n;
            delete[] d;
            if ( TMP==0 )
            {
                TMP = MOD;
                cout << "Point " << i << "/" << V << ", or "<< (i*100.0 / V) << "%\r"; cout.flush();
            }
            TMP --;
        }
        cout << "Average nr of neighbors: " << (float) total / points.size() << endl;

    }

    annDeallocPts(pa);
}

static double gauss(double x, double miu, double sigma)
{
    double tmp = ((x-miu)/sigma);
    return exp(- .5 * tmp * tmp);
}

static void normalize(std::vector<double>& v)
{
    double s = 0;
    for (size_t i=0; i<v.size(); ++i)
        s += v[i];
    for (size_t i=0; i<v.size(); ++i)
        v[i] /= s;
}

string tostr(int x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}

void FHGraph::do_gauss(double sigma)
{
    edges.reserve( V * nr_neighbors);
    list<he>::iterator k, j;
    vector<double> gauss_weight, edge_weight;
    edge e;
#pragma omp parallel for private(j, k, e, gauss_weight, edge_weight) schedule(dynamic)
    for (int i=0; i<V; ++i)
    {
        for (j=adjency_list[i].begin();
             j!=adjency_list[i].end();
             j++)
        {

            gauss_weight.clear();
            edge_weight.clear();

            for (k=adjency_list[i].begin();
                 k!=adjency_list[i].end();
                 ++k)
            {
                gauss_weight.push_back(gauss(k->w, j->w, sigma));
                edge_weight.push_back(k->w);
            }
            for (k=adjency_list[j->x].begin();
                 k!=adjency_list[j->x].end();
                 ++k)
            {
                gauss_weight.push_back(gauss(k->w, j->w, sigma));
                edge_weight.push_back(k->w);
            }
            normalize(gauss_weight);

            e.a = i; e.b = j->x;
            e.w = 0;
            for (size_t k=0; k<edge_weight.size(); ++k)
                e.w += gauss_weight[k] * edge_weight[k];

#pragma omp critical
            {
                edges.push_back(e);
            }
        }
    }
}

void FHGraph::without_gauss()
{
    edges.reserve( V * nr_neighbors);
    list<he>::iterator j;
    edge e;

    for (int i=0; i<V; ++i)
    {
        for (j=adjency_list[i].begin();
             j!=adjency_list[i].end();
             j++)
        {
            e.a = i; e.b = j->x; e.w = j->w;
            edges.push_back(e);
        }
    }
}

edge* FHGraph::getGraph()
{
    edge* ret = new edge[edges.size()];
    for (size_t i=0; i<edges.size(); ++i)
        ret[i] = edges[i];
    return ret;
}

Point FHGraph::operator[](int index)
{
    return points[index];
}

int FHGraph::getNumPoints()
{
    return V;
}

int FHGraph::getNumEdges()
{
    return edges.size();
}

template<typename T>
void vectorFree(T& t) {
    T tmp;
    t.swap(tmp);
}

void FHGraph::dispose() {
    vectorFree(edges);
    vectorFree(points);
    vectorFree(adjency_list);
}

