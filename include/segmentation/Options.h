#ifndef __OPTIONS_H_
#define __OPTIONS_H_

#include <string>

#include "slam6d/scan.h"

struct Options
{
	Options(int argc, char* argv[]);
	void usage();

	reader_type type;
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
	bool do_out;

	int reserve;

	int neighbors;
	float radius;
};

#endif // __OPTIONS_H_
