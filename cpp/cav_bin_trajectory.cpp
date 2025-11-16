#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// util functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bin_trajectory_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segments", new ParserFilename("input table with binned segments"), true);
  params.add_parser("bin_field", new ParserString("bin field"), true);
  params.add_parser("ofn", new ParserFilename("output matrix with bin reads counts per sample"), true);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int bin_trajectory_main(const char* name, int argc, char **argv)
{
  Parameters params;
  bin_trajectory_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_segments = params.get_string("ifn_segments");
  string ofn = params.get_string("ofn");
  string bin_field = params.get_string("bin_field");
  
  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn_libs, ifns, library_ids);

  int nlibs = ifns.size();

  cout << "number of libraries: " << nlibs << endl;
  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cav.load(ifns[i]);
  }

  // bin -> segments
  map< string, vector< Segment > > bins;
  read_binned_segments(ifn_segments, bin_field, bins);

  cout << "saving count matrix: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////
  // print header
  ///////////////////////////////////////////////////////////////////////////////////////

  out << bin_field;
  for (int i=0; i<nlibs; ++i)
    out << "\t" << "t" << i+1;
  out << endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // one line per bin
  ///////////////////////////////////////////////////////////////////////////////////////

  for (map< string, vector< Segment > >::iterator it=bins.begin(); it!=bins.end(); ++it) {
    string bin = (*it).first;
    vector< Segment >& segs = (*it).second;
    out << bin;
    for (int i=0; i<nlibs; ++i) {
      VariationSet& cav = cavs[i];
      int count = cav.get_segments_coverage(segs);
      out << "\t" << count;
    }
    out << endl;
  }
  
  return 0;
}
