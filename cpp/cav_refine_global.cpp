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
#include "Dissolve.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specific function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void refine_contig_global(int& segment_index, ofstream& out_segments, ofstream& out_contigs,
			  string contig, vector < VariationSet >& cavs, double outlier_fraction,
			  double min_P, double pseudo_count, string weight_style, int min_center_seg_len)
{
  Dissolve dissolve(cavs.size(), pseudo_count, weight_style);
  dissolve.set_counts(contig, cavs);
  dissolve.init_once();
  dissolve.reduce_center(outlier_fraction, min_P);

  int contig_length = dissolve.get_contig_length();
  vector<DissolveSegment> segments;
  dissolve.get_segments(segments, min_center_seg_len);

  // output segments, switching to 1-based coords and []-style segments
  for (unsigned int i=0; i<segments.size(); ++i) {
    out_segments << "s" + to_string(segment_index) << "\t" << contig << "\t";
    out_segments << segments[i].start+1 << "\t" << segments[i].end << "\t" << (segments[i].is_center ? "F" : "T")
		 << endl;
    segment_index++;
  }

  // contig summary
  int nt_outlier = 0;
  int count_outlier = 0;
  for (unsigned int i=0; i<segments.size(); ++i) {
    if (!segments[i].is_center) {
      nt_outlier += segments[i].end - segments[i].start;
      count_outlier++;
    }
  }

  out_contigs << contig << "\t" << contig_length << "\t" << count_outlier << "\t" << nt_outlier << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void refine_global_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("Table with multiple COV files"), true);
  params.add_parser("ofn_segments", new ParserFilename("output segment file"), true);
  params.add_parser("ofn_contigs", new ParserFilename("output contig file"), true);
  params.add_parser("outlier_fraction", new ParserDouble("Fraction of outliers tested every round", 0.01), false);
  params.add_parser("p_value", new ParserDouble("Min Chi-Square P-value", 0.001), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("max_lib_count", new ParserInteger("Use only n first libs (0 is all)", 0), false);
  params.add_parser("weight_style", new ParserString("Weight function of center positions (uniform|marginal)", "uniform"), false);
  params.add_parser("min_center_segment_length", new ParserInteger("Min length of center segments", 100), false);

  if (argc == 1) {
    params.usage(name);
    cout << " weight_style=uniform: uniform weights" << endl;
    cout << " weight_style=marginal: weight distributes as the z-score of the marginal coverage" << endl;
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int refine_global_main(const char* name, int argc, char **argv)
{
  Parameters params;
  refine_global_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  double outlier_fraction = params.get_double("outlier_fraction");
  double min_P = params.get_double("p_value");
  double pseudo_count = params.get_double("pseudo_count");
  int max_lib_count = params.get_int("max_lib_count");
  string ofn_segments = params.get_string("ofn_segments");
  string ofn_contigs = params.get_string("ofn_contigs");
  string weight_style = params.get_string("weight_style");
  int min_center_seg_len = params.get_int("min_center_segment_length");
  
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

  ofstream out_segments(ofn_segments.c_str(), ios::out);
  massert(out_segments.is_open(), "could not open file %s", ofn_segments.c_str());
  out_segments << "segment\tcontig\tstart\tend\tis_outlier" << endl;

  ofstream out_contigs(ofn_contigs.c_str(), ios::out);
  massert(out_contigs.is_open(), "could not open file %s", ofn_contigs.c_str());
  out_contigs << "contig\tlength\toutlier_count\toutlier_nt" << endl;

  cout << "number of contigs: " <<  contigs.size() << endl;
  int segment_index = 1;
  for (unsigned int i=0; i<contigs.size(); ++i) {
    // cout << contigs[i] << endl;
    refine_contig_global(segment_index, out_segments, out_contigs, contigs[i], cavs, outlier_fraction,
			 min_P, pseudo_count, weight_style, min_center_seg_len);
    if (i % 1000 == 0)
      cout << "progress: " << i << "/" << contigs.size() << endl;
  }

  out_segments.close();
  out_contigs.close();

  return 0;
}
