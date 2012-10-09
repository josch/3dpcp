#include "segmentation/Timer.h"
#include "slam6d/globals.icc"
#include <sstream>
using namespace std;

string Timer::prefix_a = " └─";
string Timer::prefix_n = " │ ";
int Timer::depth = 0;

Timer::Timer(string format)
{
	gettimeofday(&time, 0);
	this->format = format;

	for (int i=0; i<depth; ++i)
		cout << prefix_n;
	cout << this->format.substr(0, this->format.find('#')) << endl;
	
	this->format.replace(0, format.find('#'), "");
	depth ++;
}

Timer::~Timer()
{
	struct timeval t;
	gettimeofday(&t, 0);
	float tv = (t.tv_sec-time.tv_sec) * 1000.0 + (t.tv_usec-time.tv_usec) * .001;
	
	stringstream ss;
	ss << tv;
	format.replace(format.find('#'), 1, ss.str());

	for (int i=0; i<depth-1; ++i)
		cout << prefix_n;
	cout << prefix_a;
	cout << format << endl;
	
	depth --;
}
