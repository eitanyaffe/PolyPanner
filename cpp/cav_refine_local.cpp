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
#include <algorithm>

#include <thread>
#include <mutex>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"
#include "RefineLocal.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// specific function
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

mutex g_refine_local_mtx; 

void refine_contig_local(vector<string>& contigs, int from_ind, int to_ind,
			 vector < VariationSet >& cavs, 
			 int margin, int min_contig_length, int min_length, int max_length,
			 double min_P, double pseudo_count, bool only_dangles,
			 int cand_step, bool debug,
			 map < string, set<Breakpoint> >& bps, int& test_count)
{  
  for (int i=from_ind; i<to_ind; ++i) {
    string contig = contigs[i];
    set<Breakpoint> breakpoints;

    if (debug)
      cout << "refining contig " << contig << endl;
    
    RefineLocal refine(contig, cavs, margin, min_length, max_length, min_P,
		       pseudo_count, cand_step, debug);
    refine.get_breakpoints(only_dangles, breakpoints);

    g_refine_local_mtx.lock();
    bps[contig] = breakpoints;
    test_count += refine.get_test_count();
    g_refine_local_mtx.unlock();
  }
}

void output_breakpoints(map < string, set<Breakpoint> >& bps,
			map < string, int >& contig_lengths,
			string ofn_breakpoints, string ofn_segments, string ofn_contigs)
{
  ofstream out_breakpoints(ofn_breakpoints.c_str(), ios::out);
  massert(out_breakpoints.is_open(), "could not open file %s", ofn_breakpoints.c_str());
  out_breakpoints << "contig\tcoord\trequested_length\tp_value\tindex";
  out_breakpoints << "\tleft_side_start\tleft_side_end\tright_side_start\tright_side_end\tread_support" << endl;
  
  ofstream out_segments(ofn_segments.c_str(), ios::out);
  massert(out_segments.is_open(), "could not open file %s", ofn_segments.c_str());
  out_segments << "segment\tcontig\tstart\tend" << endl;

  ofstream out_contigs(ofn_contigs.c_str(), ios::out);
  massert(out_contigs.is_open(), "could not open file %s", ofn_contigs.c_str());
  out_contigs << "contig\tlength\tseg_count" << endl;

  cout << "number of contigs: " << bps.size() << endl;
  int segment_index = 1;
  for (map < string, set<Breakpoint> >::iterator it=bps.begin(); it != bps.end(); ++it) {
    string contig = (*it).first;
    set<Breakpoint>& breakpoints = (*it).second;
    int seg_start = 0;
    int seg_count = 0;
    
    // cout << "contig=" << contig << ", breakpoint_count=" << breakpoints.size() << endl;
    for (set<Breakpoint>::iterator it=breakpoints.begin(); it != breakpoints.end(); ++it) {
      Breakpoint bp = *it;
      int coord = bp.coord;
      
      // output breakpoint
      out_breakpoints << contig << "\t" << coord+1 << "\t" << bp.side_requested_length << "\t" << bp.p_value;
      out_breakpoints << "\t" << bp.index << "\t" << bp.left_side_start+1 << "\t" << bp.left_side_end+1;
      out_breakpoints << "\t" << bp.right_side_start+1 << "\t" << bp.right_side_end+1;
      out_breakpoints << "\t" << bp.read_support << endl;
      
      // output segment
      out_segments << "s" + to_string(segment_index) << "\t" << contig << "\t";
      out_segments << seg_start+1 << "\t" << coord+1 << endl;
      seg_start = coord;
      segment_index++;
      seg_count++;
    }

    int contig_length = contig_lengths[contig];
    
    // handle last segment
    if (seg_start < contig_length) {
      out_segments << "s" + to_string(segment_index) << "\t" << contig << "\t";
      out_segments << seg_start+1 << "\t" << contig_length << endl;
      segment_index++;
      seg_count++;
    }
    
    // contig summary
    out_contigs << contig << "\t" << contig_length << "\t" << seg_count << endl;
  }

  out_breakpoints.close();
  out_segments.close();
  out_contigs.close();
}

void refine_apply_BH(double alpha, int test_count, map < string, set<Breakpoint> >& bps)
{
  vector<double> pvals;
  int count = 0;
  for (map < string, set<Breakpoint> >::iterator it=bps.begin(); it != bps.end(); ++it) {
    string contig = (*it).first;
    set<Breakpoint>& breakpoints = (*it).second;
    for (set<Breakpoint>::iterator it=breakpoints.begin(); it != breakpoints.end(); ++it) {
      Breakpoint bp = *it;
      pvals.push_back(bp.p_value);
      count++;
    }
  }
  if (pvals.size() == 0) {
    cout << "No candidate sites found" << endl;
    return;
  }
  
  cout << "Applying BH, number of tests: " << test_count << ", number of candidates: " << count << endl;
  double p_min_fdr = apply_Benjamini_Hochberg(pvals, alpha, test_count);
  cout << "Min corrected p-value: " << p_min_fdr << endl;
  
  set<Breakpoint> empty_s;
  int count_accepted = 0;
  map < string, set<Breakpoint> > result;
  for (map < string, set<Breakpoint> >::iterator it=bps.begin(); it != bps.end(); ++it) {
    string contig = (*it).first;
    result[contig] = empty_s;
    set<Breakpoint>& breakpoints = (*it).second;
    for (set<Breakpoint>::iterator it=breakpoints.begin(); it != breakpoints.end(); ++it) {
      Breakpoint bp = *it;
      if (bp.p_value > p_min_fdr)
	continue;
      result[contig].insert(bp);
      count_accepted++;
    }
  }
  cout << "Number of candidate sites: " << count << endl;
  cout << "Number of sites that pass FDR filter: " << count_accepted << endl;
  bps = result;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void refine_local_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("threads", new ParserInteger("Number of threads", 4), true);
  params.add_parser("fdr", new ParserDouble("False discovery rate", 0.001), false);
  params.add_parser("ofn_segments", new ParserFilename("output segment file"), true);
  params.add_parser("ofn_contigs", new ParserFilename("output contig file"), true);
  params.add_parser("ofn_breakpoints", new ParserFilename("output breakpoint file"), true);
  params.add_parser("p_value", new ParserDouble("Min Chi-Square P-value", 0.001), false);
  params.add_parser("margin", new ParserInteger("Margin away from breakpoints", 10), false);
  params.add_parser("cand_step", new ParserInteger("window step size for non-dangle candidates", 10), false);
  params.add_parser("min_contig_length", new ParserInteger("Min contig length", 1000), false);
  params.add_parser("min_length", new ParserInteger("Min breakpoint side length", 150), false);
  params.add_parser("max_length", new ParserInteger("Max breakpoint side length", 2400), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("max_lib_count", new ParserInteger("Stop after X libs for debugging (if set to 0 then use all)", 0.1), false);
  params.add_parser("only_dangles", new ParserBoolean("Check only coords which end a trimmed read", true), false);
  params.add_parser("contig", new ParserString("Only single contig", "NA"), false);

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

int refine_local_main(const char* name, int argc, char **argv)
{
  Parameters params;
  refine_local_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  double min_P = params.get_double("p_value");
  double alpha = params.get_double("fdr");
  double pseudo_count = params.get_double("pseudo_count");
  int margin = params.get_int("margin");
  int min_contig_length = params.get_int("min_contig_length");
  int min_length = params.get_int("min_length");
  int max_length = params.get_int("max_length");
  int max_lib_count = params.get_int("max_lib_count");
  int cand_step = params.get_int("cand_step");
  string only_contig = params.get_string("contig");
  string ofn_breakpoints = params.get_string("ofn_breakpoints");
  string ofn_segments = params.get_string("ofn_segments");
  string ofn_contigs = params.get_string("ofn_contigs");
  int thread_count = params.get_int("threads");
  bool only_dangles = params.get_bool("only_dangles");

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

  vector<string> all_contigs = cavs[0].get_contigs();
  vector<string> contigs;
  map < string, int > contig_lengths;

  for (auto contig : all_contigs) {
    int clen = cavs[0].get_covs()[contig].size();
    if (clen >= min_contig_length)
      contigs.push_back(contig);
    contig_lengths[contig] = clen;
  }

  if (only_contig != "NA" && contig_lengths.find(only_contig) == contig_lengths.end()) {
    cout << "contig not found: " << only_contig << endl;
    exit(1);
  }
  
  bool debug = only_contig != "NA";
  if (debug) {
    thread_count = 1;
    contigs.resize(1);
    contigs[0] = only_contig;
  }
  
  int n_contigs = contigs.size();
  cout << "number of contigs: " <<  n_contigs << endl;

  // collect all breakpoints
  map < string, set<Breakpoint> > bps;
  int test_count = 0;
  
  thread_count = min(thread_count, n_contigs);
  int step = floor(n_contigs / thread_count);
  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : contigs.size();
    threads.push_back(thread(refine_contig_local, 
			     std::ref(contigs), from_ind, to_ind, std::ref(cavs),
			     margin, min_contig_length, min_length, max_length,
			     min_P, pseudo_count, only_dangles, cand_step, debug,
			     std::ref(bps), std::ref(test_count)));
  }
  cout << "waiting for " << thread_count << " threads to finish\n";
  for (auto& th : threads) th.join();

  refine_apply_BH(alpha, test_count, bps);
  output_breakpoints(bps, contig_lengths, ofn_breakpoints, ofn_segments, ofn_contigs);

  return 0;
}
