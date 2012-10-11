#include "segmentation/Options.h"

#include <getopt.h>
#include <iostream>
#include <sstream>
#include <cstdlib>
using namespace std;

template<typename T> T parse(string X)
{
	stringstream ss(X);
	T result;
	ss >> result;
	return result;
}

Options::Options(int argc, char* argv[])
{
	type = RXP;
	minDist = -1;
	maxDist = -1;
	start = 0;
	end = -1;
	sigma = 1;
	k = 1;
	min_size = 0;
	outdir = "./segments/";
	do_out = false;
	reserve = -1;
	eps = 1.0;
  scanserver = false;
	
	neighbors = -1;
	radius = -1;

	char c;
	while ( (c=getopt(argc, argv, "hxf:m:M:s:e:S:K:I:o:r:n:R:A:")) != -1 )
		switch ( c )
		{
			case 'f':
				try {
          type = formatname_to_io_type(optarg);
        } catch (...) { // runtime_error
          cerr << "Format " << optarg << " unknown." << endl;
          abort();
        }
        break;
			case 'm':
				minDist = parse<float>(optarg);
				break;
			case 'M':
				maxDist = parse<float>(optarg);
				break;
			case 's':
				start = parse<int>(optarg);
				break;
			case 'e':
				end = parse<int>(optarg);
				break;
			case 'S':
				sigma = parse<float>(optarg);
				break;
			case 'K':
				k = parse<float>(optarg);
				break;
			case 'I':
				min_size = parse<int>(optarg);
				break;
			case 'h':
				usage();
				exit(0);
			case 'o':
				do_out = true;
				outdir = optarg;
				break;
			case 'r':
				radius = parse<float>(optarg);
				break;
			case 'n':
				neighbors = parse<int>(optarg);
				break;
			case 'R':
				reserve = parse<int>(optarg);
				cerr << "std::vector reserve hint " << reserve << endl;
				break;
			case 'A':
				eps = parse<float>(optarg);
				break;
      case 'x':
        scanserver = true;
        break;
      default:
	     abort ();
	  }

	if ( optind < argc )
		dir = argv[optind];
	else
	{
		usage();
		exit(1);
	}
	if ( dir[dir.length()-1] != '/' ) dir += "/";
	if ( outdir[outdir.length()-1] != '/' ) outdir += "/";

  boost::filesystem::path boost_dir(outdir);
  if (!boost::filesystem::create_directory(boost_dir)) {
    cerr << "Couldn't create directory " << outdir << endl;
  }

  cout << endl << "Input directory is: " << dir << endl;
  cout << "Output directory is: " << outdir << endl;
}

void Options::usage()
{
	cerr << "\n\nUsage: ./bin/fh_segmentation [OPTIONS] SCAN\n"
		 << "Available options:\n"
		 << "-f TYPE   scan type; currently supported: rxp\n"
     << "-x        enables scanserver\n"
		 << "-m DIST   the minimum distance used by the reader\n"
		 << "-M DIST   the maximum distance used by the reader\n"
		 << "-s FRAME  the starting frame\n"
		 << "-e FRAME  the ending frame\n"
		 << "-S SIGMA  the sigma value used in Gaussian smoothing\n"
		 << "-K K      the K value used in the FH segmentation\n"
		 << "-I SIZE   the min size of a segment\n"
		 << "-o DIR    the output directory (don't use for no output)\n"
		 << "-r RADIUS use radius search, specify the radius, if RADIUS > 0 [default RADIUS=-1]\n"
		 << "-n NR     use approximate NR-nearest neighbors search or limit the number of points\n"
		 << "          returned by the radius search (it will pick NR points randomly)\n"
		 << "-R SIZE   the size of pre-reserved std::vectors\n"
		 << "-A EPS    the error used by the AKNN algorithm\n"
		 << "-h        this help message\n\n";
}
