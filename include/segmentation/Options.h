#ifndef __OPTIONS_H_
#define __OPTIONS_H_

#include <string>

#include "slam6d/scan.h"
#include "slam6d/globals.icc"
#include "slam6d/io_utils.h"

struct Options
{
	Options(int argc, char* argv[]);
	void usage();

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
	bool do_out;

	int reserve;

	int neighbors;
	float radius;
};

#endif // __OPTIONS_H_
