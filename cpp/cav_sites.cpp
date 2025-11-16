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

#include <boost/math/distributions/chi_squared.hpp>
using namespace boost::math;

#include "util.h"
#include "Params.h"
#include "VariationSet.h"
#include "Dissolve.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// sites functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void sites_read_segments(string fn, map< string, set< Segment > >& segs)
{
 cout << "reading binned segment table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int segment_ind = get_field_index("segment", fields);
  int contig_ind = get_field_index("contig", fields);
  int start_ind = get_field_index("start", fields);
  int end_ind = get_field_index("end", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string seg_id = fields[segment_ind];
    string contig = fields[contig_ind];
    int start = safe_string_to_int(fields[start_ind], "coordinate") - 1;
    int end = safe_string_to_int(fields[end_ind], "coordinate") - 1;

    Segment seg(seg_id, contig, start, end);
    segs[contig].insert(seg);
  }
}

string sites_get_seg_id(string contig, int coord, map< string, set< Segment > >& segs)
{
  if (segs.find(contig) == segs.end())
    return "none";
  set<Segment>& csegs = segs[contig];

  auto it = upper_bound(csegs.begin(), csegs.end(), coord,
  			[](int coord, const Segment& seg)
			{ return coord < seg.end; });
  if (it == csegs.end()) {
    // cout << "not found" << endl;
    return "none";
  } else {
    // cout << "segment=" << it->id << ", " << it->start << "-" << it->end << endl;
    return it->id;
  }
}

bool sites_is_internal(string contig, int coord, int margin,
		       map< string, set< Segment > >& segs)
{
  if (segs.find(contig) == segs.end())
    return false;
  set<Segment>& csegs = segs[contig];

  auto it = upper_bound(csegs.begin(), csegs.end(), coord,
  			[](int coord, const Segment& seg)
			{ return coord < seg.end; });
  if (it == csegs.end()) {
    return false;
  } else {
    return coord >= (it->start + margin) && coord <= (it->end - margin);
  }
}

void annotate_vars(string ofn, double pseudo_count,
		   vector<VariationSet>& vsets, VariationSet& vset_sum,
		   int regional_window, int margin,
		   map<string, int>& contig_map,
		   map< string, string >& bin_map,
		   map< string, set< Segment > >& segs,
		   map< string, map< int, set< Variation > > >& vars)
{
  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "vid" << "\t" << "contig" << "\t" << "coord" << "\t" << "variant" << "\t";
  out << "segment" << "\t" << "bin" << "\t";
  out << "var_count" << "\t" << "total_count" << "\t";
  out << "n_samples" << "\t" << "is_internal" << "\t";
  out << "variant_p" << "\t";
  out << "regional_p" << "\t";
  out << "comp_p" << endl;
  
  int N = vsets.size();

  int vcount = 1;
  for (map< string, map< int, set< Variation > > >::iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    map< int, set< Variation > >& vars_contig = (*it).second;
    for (map< int, set< Variation > >::iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      set< Variation > vars_coord = (*jt).second;

      string segment_id = sites_get_seg_id(contig, coord, segs);
      string bin = (bin_map.find(segment_id) != bin_map.end()) ? bin_map[segment_id] : "none";
      
      // get total coverage
      int coord_cov = vset_sum.get_coverage(contig, coord);

      // add the reference if has any coverage
      // int ref_count = vset_sum.get_ref_count(contig, coord);
      // Variation ref_var;
      // if (ref_count > 0)
      // 	vars_coord.insert(ref_var);

      // correspondance between local and regional total coverage vectors
      vector<double> loc_v(N);
      vector<double> reg_v(N);
      for (int i=0; i<N; ++i) {
	loc_v[i] = vsets[i].get_coverage(contig, coord) + pseudo_count;
	reg_v[i] = vsets[i].get_regional_coverage(contig, coord, regional_window, margin) + pseudo_count;
      }
      double p_reg = chi_square_test_p(loc_v, reg_v);
      // double rho_reg = pearson_test(loc_v, reg_v, NULL);

      // near contig side
      massert(contig_map.find(contig) != contig_map.end(), "contig not found");
      // int clen = contig_map[contig];

      // margin from segment sides
      int is_internal = sites_is_internal(contig, coord, margin, segs);
      // bool is_side = coord < contig_side_size || coord > (clen-contig_side_size);
      
      for (set< Variation >::iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt);
	int var_cov = vset_sum.get_var_count(contig, coord, var);

	// count number of unique samples
	int n_samples = 0;
	
	// variant vector
	vector<double> var_v(N);
	vector<double> comp_v(N);
	for (int i=0; i<N; ++i) {
	  int v_count = vsets[i].get_var_count(contig, coord, var);
	  if (v_count > 0)
	    n_samples++;
	  var_v[i] = v_count + pseudo_count;
	  comp_v[i] = loc_v[i] - var_v[i] + pseudo_count;
	}
	// compare variant and local total coverage 
	double p_var = chi_square_test_p(loc_v, var_v);

	// compare complement of variant with regional
	double p_comp = chi_square_test_p(comp_v, reg_v);

	string vid = "v" + to_string(vcount);
	
	// output
	out << vid << "\t" << contig << "\t" << coord+1 << "\t" << var << "\t";
	out << segment_id << "\t" << bin << "\t";
	out << var_cov << "\t" << coord_cov << "\t";
	out << n_samples << "\t" << (is_internal ? "T" : "F") << "\t";
	out << p_var << "\t" << p_reg << "\t";
	// out << rho_reg << "\t";
	out << p_comp << endl;

	vcount++;
      }
    }
  }

  out.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void sites_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("table with multiple POP files"), true);
  params.add_parser("ifn_sites", new ParserFilename("input site table"), true);
  params.add_parser("ifn_contigs", new ParserFilename("contig table"), true);
  params.add_parser("ifn_segments", new ParserFilename("segment table"), true);
  params.add_parser("ifn_segment_bins", new ParserFilename("segment-bin table"), true);
  params.add_parser("regional_window", new ParserInteger("window used for regional coverage (bp)", 1000), false);
  params.add_parser("margin", new ParserInteger("flag sites close to contig end (bp)"), 10), 
  params.add_parser("pseudo_count", new ParserDouble("pseudo-count for coverage counts", 0.1), false);
  params.add_parser("ofn", new ParserFilename("output site file"), true);

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

int sites_main(const char* name, int argc, char **argv)
{
  Parameters params;
  sites_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_sites = params.get_string("ifn_sites");
  string ifn_contigs = params.get_string("ifn_contigs");
  string ifn_segments = params.get_string("ifn_segments");
  string ifn_segment_bins = params.get_string("ifn_segment_bins");
  double pseudo_count = params.get_double("pseudo_count");
  int regional_window = params.get_int("regional_window");
  int margin = params.get_int("margin");
  string ofn = params.get_string("ofn");


  map<string, int> contig_map;
  read_contig_length(ifn_contigs, contig_map);
  cout << "number of contigs: " << contig_map.size() << endl;

  map< string, set< Segment > > segs;
  sites_read_segments(ifn_segments, segs);

  // map from segment to bin
  map< string, string > bin_map;
  read_string_map(ifn_segment_bins, "segment", "bin", bin_map);
  
  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn_libs, ifns, library_ids);
  
  map< string, map< int, set< Variation > > > vars;
  read_sites(ifn_sites, vars);
  
  vector<VariationSet> vsets(ifns.size());
  VariationSet sum_vset;
  for (unsigned int i=0; i<ifns.size(); ++i) {
    VariationSet& varset_i = vsets[i];
    varset_i.load(ifns[i]);
    if (i == 0)
      sum_vset = varset_i;
    else
      sum_vset = sum_vset + varset_i;

  }

  annotate_vars(ofn, pseudo_count, vsets, sum_vset, regional_window, margin,
		contig_map, bin_map, segs, vars);

  return 0;
}
