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

#include "util.h"
#include "Params.h"
#include "VariationSet.h"
#include "RefineLocal.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specific function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dump_refine_local(ofstream& out, string contig, vector < VariationSet >& cavs,
		       int margin, int min_length, int max_length, double pseudo_count)
{
  RefineLocal refine(contig, cavs, margin, min_length, max_length, 0.1, pseudo_count, 10, true);
  refine.dump_scores(out);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dump_local_scores_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("Table with multiple COV files"), true);
  params.add_parser("contig", new ParserString("contig"), true);
  params.add_parser("ofn", new ParserFilename("output file"), true);
  params.add_parser("margin", new ParserInteger("Margin away from breakpoints", 10), false);
  params.add_parser("min_length", new ParserInteger("Min distance away from candidate coordinate", 150), false);
  params.add_parser("max_length", new ParserInteger("Max distance away from candidate coordinate", 2400), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("max_lib_count", new ParserInteger("Stop after X libs for debugging (if set to 0 then use all)", 0.1), false);

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

int dump_local_scores_main(const char* name, int argc, char **argv)
{
  Parameters params;
  dump_local_scores_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string contig = params.get_string("contig");
  double pseudo_count = params.get_double("pseudo_count");
  int margin = params.get_int("margin");
  int min_length = params.get_int("min_length");
  int max_length = params.get_int("max_length");
  int max_lib_count = params.get_int("max_lib_count");
  string ofn = params.get_string("ofn");
  
  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn, ifns, library_ids);

  int nlibs = ifns.size();
  if (max_lib_count > 0 && max_lib_count < nlibs)
    nlibs = max_lib_count;
  
  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cav.load(ifns[i]);
  }

  vector<string> contigs = cavs[0].get_contigs();
  massert(find(contigs.begin(), contigs.end(), contig) != contigs.end(), "contig not found");

  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  dump_refine_local(out, contig, cavs, margin, min_length, max_length, pseudo_count);
  
  //cout << "number of contigs: " <<  contigs.size() << endl;
  //for (unsigned int i=0; i<contigs.size(); ++i) {
  //  dump_refine_local(out, contigs[i], cavs, margin, min_length, max_length, pseudo_count);
  //}

  out.close();

  return 0;
}
