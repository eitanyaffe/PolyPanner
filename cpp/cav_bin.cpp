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

#include "BinMatrix.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bin functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void init_segment(BinSegmentSide& seg, vector < VariationSet >& cavs,
		  double pseudo_count, int margin)
{
  int n_libs = cavs.size();
  string contig = seg.contig;
  int len = seg.side_length;
  int start = seg.is_left_side ? seg.start+margin : seg.end-len-margin;
  int end = start+len;

  seg.counts.resize(n_libs);

  int total_reads = 0;
  int total_samples = 0;
  for (unsigned int k=0; k<cavs.size(); ++k) {
    VariationSet& cav = cavs[k];
    int count = cav.get_segment_coverage(contig, start, end);
    seg.counts[k] = count + pseudo_count;
    total_reads += count;
    total_samples += (count > 0 ? 1 : 0);
  }
  seg.supporting_reads = total_reads;
  seg.supporting_samples = total_samples;
}

void read_segments(string fn, vector < VariationSet >& cavs,
		   vector< BinSegmentSide >& segs, vector< pair<int,int> >& seg_pairs,
		   int max_side_length, int min_seg_length, double pseudo_count, int margin,
		   int min_read_support, int min_samples,
		   string only_contig)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int id_ind = get_field_index("segment", fields);
  int contig_ind = get_field_index("contig", fields);
  int start_ind = get_field_index("start", fields);
  int end_ind = get_field_index("end", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int start = stoi(fields[start_ind]);
    int end = stoi(fields[end_ind]);

    int seg_length = end - start + 1;
    int side_length = min(seg_length, max_side_length) - 2*margin;

    if (seg_length < min_seg_length)
      continue;

    if (only_contig != "NA" && contig != only_contig)
      continue;
    
    BinSegmentSide seg_left(id, contig, start-1, end-1, true, side_length);
    BinSegmentSide seg_right(id, contig, start-1, end-1, false, side_length);

    init_segment(seg_left, cavs, pseudo_count, margin);
    init_segment(seg_right, cavs, pseudo_count, margin);

    bool select_left = seg_left.supporting_reads >= min_read_support && seg_left.supporting_samples >= min_samples;
    bool select_right = seg_right.supporting_reads >= min_read_support && seg_right.supporting_samples >= min_samples;
    
    // add pair only of both sides are selected
    if (select_left && select_right) {
      int index = segs.size();
      seg_pairs.push_back(make_pair(index, index+1));
    }
    
    if (select_left)
      segs.push_back(seg_left);
    if (select_right)
      segs.push_back(seg_right);
  }
}

void output_bins(ofstream& out, vector< BinSegmentSide >& segs, vector<int>& bins)
{
  out << "segment\tcontig\tstart\tend\tbin" << endl;

  // keep track of segment ids, since each segment has two sides
  set<string> visited_segment_ids;
  for (unsigned int i=0; i<segs.size(); ++i) {
    BinSegmentSide& seg = segs[i];
    if (visited_segment_ids.find(seg.id) == visited_segment_ids.end()) {
      out << seg.id << "\t" << seg.contig << "\t" << seg.start+1 << "\t" << seg.end+1 << "\t" << bins[i]+1 << endl;
      visited_segment_ids.insert(seg.id);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bin_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segments", new ParserFilename("input segment file"), true);
  params.add_parser("ofn", new ParserFilename("output binning table"), true);
  params.add_parser("threads", new ParserInteger("Number of threads used", 40), false);
  params.add_parser("p_value", new ParserDouble("Chi-Square P-value clustering threshold", 0.05), false);
  params.add_parser("margin", new ParserInteger("Margin away from breakpoints", 10), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("min_segment_length", new ParserDouble("Minimal segment length", 200), false);
  params.add_parser("max_lib_count", new ParserInteger("Use only n first libs (0 is all)", 0), false);
  params.add_parser("max_side_length", new ParserInteger("Maximal distance from fragment side used for coverage",
							 1000), false);
  params.add_parser("min_read_support", new ParserInteger("Minimal total number of reads supporting candidate side",
							  200), false);
  params.add_parser("min_samples", new ParserInteger("Minimal total number of samples with one or more supporting reads",
							  200), false);
  params.add_parser("max_coverage_factor", new ParserDouble("Threshold on coverage ratio between adjacent segments",
							    1.5), false);
  params.add_parser("only_contig", new ParserString("limit to contig", "NA"), false);

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

int bin_main(const char* name, int argc, char **argv)
{
  Parameters params;
  bin_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_segments = params.get_string("ifn_segments");
  string ofn = params.get_string("ofn");
  int thread_count = params.get_int("threads");
  double p_threshold = params.get_double("p_value");
  int margin = params.get_int("margin");
  double pseudo_count = params.get_double("pseudo_count");
  double min_seg_length = params.get_double("min_segment_length");
  double max_coverage_factor = params.get_double("max_coverage_factor");
  int max_side_length = params.get_int("max_side_length");
  int max_lib_count = params.get_int("max_lib_count");
  int min_read_support = params.get_int("min_read_support");
  int min_samples = params.get_int("min_samples");
  string only_contig = params.get_string("only_contig");

  bool debug_mode = only_contig != "NA";
  
  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn_libs, ifns, library_ids);

  int nlibs = ifns.size();
  if (max_lib_count > 0 && max_lib_count < nlibs)
    nlibs = max_lib_count;

  cout << "number of libraries: " << nlibs << endl;
  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cav.load(ifns[i]);
  }
  // segment sides
  vector< BinSegmentSide > segs;

  // indices of paired sides
  vector< pair<int,int> > seg_pairs;
  read_segments(ifn_segments, cavs, segs, seg_pairs,
		max_side_length, min_seg_length, pseudo_count, margin,
		min_read_support, min_samples,
		only_contig);

  massert(segs.size() > 0, "no segments match criteria");
  cout << "number of segments that match criteria: " << segs.size() << endl;

  BinMatrix binner(segs, nlibs);

  // init matrix
  cout << "initialing matrix" << endl;
  binner.init_matrix(thread_count);

  ofstream* dbg_out = NULL;
  if (debug_mode) {
    string mat_ofn = ofn + ".mat";
    dbg_out = new ofstream(mat_ofn.c_str(), ios::out);
  }

  // cluster segments
  vector<int> bins;
  cout << "binning contigs" << endl;
  int num = binner.cluster_segments(bins, seg_pairs, p_threshold, max_coverage_factor, dbg_out);
  cout << "number of bins: " << num << endl;

  // output results
  ofstream out(ofn.c_str(), ios::out);
  output_bins(out, segs, bins);
  out.close();

  if (debug_mode) {
    dbg_out->close();
    delete dbg_out;
  }
  
  return 0;
}
