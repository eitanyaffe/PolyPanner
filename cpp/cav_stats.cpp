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

void stats_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ofn", new ParserFilename("output read depth and max read length per lib"), true);

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

int cav_stats_main(const char* name, int argc, char **argv)
{
  Parameters params;
  stats_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string ofn = params.get_string("ofn");
  
  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn, ifns, library_ids);

  int nlibs = ifns.size();

  cout << "number of libraries: " << nlibs << endl;
  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cav.load(ifns[i]);
  }

  cout << "saving table: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  out << "lib\ttotal_reads\tmax_read_length\tmean_read_length" << endl;
  
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    out << "t" << i+1 << "\t" << cav.get_read_count() << "\t";
    out << cav.get_max_read_length() << "\t" << cav.get_mean_read_length() << endl;
  }
  
  return 0;
}
