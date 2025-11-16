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

void refine_bins_read_seed_bins(string fn, string omit_bin, string only_bin, map<string, string>& seg_map)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  int seg_ind = 0;
  int bin_ind = 1;
  
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string segment = fields[seg_ind];
    string bin = fields[bin_ind];
    if ((bin == omit_bin) || (only_bin != "NA" && bin != only_bin))
      continue;
    seg_map[segment] = bin;
  }
}

void init_segment_side(BinSegmentSide& seg, vector < VariationSet >& cavs,
		       double pseudo_count, int margin, string seed_bin)
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
  seg.seed_bin = seed_bin;
}

void refine_bins_read_segments(string fn, vector < VariationSet >& cavs,
			       vector< BinSegmentSide >& segs, map<string, string>& seg_map,
			       vector< pair<int,int> >& seg_pairs,
			       int max_side_length, double pseudo_count, int margin)
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

    if (seg_map.find(id) == seg_map.end())
      continue;
    
    string seed_bin = seg_map[id];
    
    int seg_length = end - start + 1;
    int side_length = min(seg_length, max_side_length) - 2*margin;

    BinSegmentSide seg_left(id, contig, start-1, end-1, true, side_length);
    BinSegmentSide seg_right(id, contig, start-1, end-1, false, side_length);

    init_segment_side(seg_left, cavs, pseudo_count, margin, seed_bin);
    init_segment_side(seg_right, cavs, pseudo_count, margin, seed_bin);

    int index = segs.size();
    
    seg_pairs.push_back(make_pair(index, index+1));
    segs.push_back(seg_left);
    segs.push_back(seg_right);
  }
}

void output_bins(ofstream& out, vector< BinSegmentSide >& segs, vector<int>& bins)
{
  // seed_bin -> count
  map <string, set<int> > ibin_map;
  for (unsigned int i=0; i<segs.size(); ++i) {
    BinSegmentSide& seg = segs[i];
    string seed_bin = seg.seed_bin;
    int bin = bins[i];
    ibin_map[seed_bin].insert(bin);
  }

  // track current number of refined bins within seed bin
  map <string, int > ibin_index;

  // map from output bin (int) to identifier (seed + index)
  map <int, string > obin_to_id;

  for (unsigned int i=0; i<segs.size(); ++i) {
    BinSegmentSide& seg = segs[i];
    string seed_bin = seg.seed_bin;
    int bin = bins[i];
    if (obin_to_id.find(bin) != obin_to_id.end())
      continue;
    massert(ibin_map.find(seed_bin) != ibin_map.end(), "seed bin %s not found", seed_bin.c_str());
    if (ibin_map[seed_bin].size() > 1) {
      int index = ibin_index[seed_bin] + 1;
      obin_to_id[bin] = "b" + seed_bin + "_" + to_string(index);
      ibin_index[seed_bin]++;
    } else {
      obin_to_id[bin] = "b" + seed_bin;
    }
  }
  
  out << "segment\tcontig\tstart\tend\tbin_org\tbin_index\tbin" << endl;

  // keep track of segment ids, since each segment has two sides
  set<string> visited_segment_ids;
  for (unsigned int i=0; i<segs.size(); ++i) {
    BinSegmentSide& seg = segs[i];
    int bin = bins[i];
    if (visited_segment_ids.find(seg.id) == visited_segment_ids.end()) {
      out << seg.id << "\t" << seg.contig << "\t" << seg.start+1 << "\t" << seg.end+1 << "\t";
      out << seg.seed_bin << "\t" << bin+1 << "\t" << obin_to_id[bin] << endl;
      visited_segment_ids.insert(seg.id);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void refine_bins_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segments", new ParserFilename("input segment file"), true);
  params.add_parser("ifn_segments_binned", new ParserFilename("input table associating segments (1st column) to bins (2nd column)"), true);
  params.add_parser("ofn", new ParserFilename("output binning table"), true);
  params.add_parser("threads", new ParserInteger("Number of threads used", 40), false);
  params.add_parser("omit_bin", new ParserString("Discard input bin with this value", ""), false);
  params.add_parser("only_bin", new ParserString("Limit to specific bin", "NA"), false);
  params.add_parser("p_value", new ParserDouble("Chi-Square P-value clustering threshold", 0.05), false);
  params.add_parser("margin", new ParserInteger("Margin away from breakpoints", 10), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("max_side_length", new ParserInteger("Maximal distance from fragment side used for coverage",
							 1000), false);

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

int refine_bins_main(const char* name, int argc, char **argv)
{
  Parameters params;
  refine_bins_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_segments = params.get_string("ifn_segments");
  string ifn_segments_binned = params.get_string("ifn_segments_binned");
  string ofn = params.get_string("ofn");
  string omit_bin = params.get_string("omit_bin");
  string only_bin = params.get_string("only_bin");
  int thread_count = params.get_int("threads");
  double p_threshold = params.get_double("p_value");
  int margin = params.get_int("margin");
  double pseudo_count = params.get_double("pseudo_count");
  int max_side_length = params.get_int("max_side_length");
  
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

  // map from segments to bins
  map<string, string> seg_map;
  
  refine_bins_read_seed_bins(ifn_segments_binned, omit_bin, only_bin, seg_map);
  cout << "number of input segments: " << seg_map.size() << endl;
  
  // segment sides
  vector< BinSegmentSide > segs;

  // indices of paired sides
  vector< pair<int,int> > seg_pairs;
  refine_bins_read_segments(ifn_segments, cavs, segs, seg_map, seg_pairs,
			    max_side_length, pseudo_count, margin);
  
  massert(segs.size() > 0, "no segments match criteria");

  bool debug = (only_bin != "NA");
  BinMatrix binner(segs, nlibs, p_threshold, debug);

  // init matrix
  cout << "initializing matrix ..." << endl;
  binner.init_matrix(thread_count);

  // cluster segments
  vector<int> bins;
  cout << "binning contigs ..." << endl;
  int num = binner.cluster_segments(bins, seg_pairs);
  cout << "number of bins: " << num << endl;

  // output results
  ofstream out(ofn.c_str(), ios::out);
  output_bins(out, segs, bins);
  out.close();

  return 0;
}
