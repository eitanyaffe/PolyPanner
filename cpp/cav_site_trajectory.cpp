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
// util functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void site_trj_read_segments(string fn, map<string, map < int, string > >& bin_map)
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
  int bin_ind = get_field_index("bin", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int start = stoi(fields[start_ind])-1;
    int end = stoi(fields[end_ind])-1;
    string bin = fields[bin_ind];
    for (int coord=start; coord<=end; ++coord)
      bin_map[contig][coord] = bin;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void site_trajectory_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segments", new ParserFilename("input segment file, with bin field"), true);
  params.add_parser("ifn_sites", new ParserFilename("input table with sites"), true);
  params.add_parser("ofn_counts", new ParserFilename("output matrix with var counts"), true);
  params.add_parser("ofn_totals", new ParserFilename("output matrix with total counts"), true);

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

int site_trajectory_main(const char* name, int argc, char **argv)
{
  Parameters params;
  site_trajectory_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_segments = params.get_string("ifn_segments");
  string ifn_sites = params.get_string("ifn_sites");
  string ofn_counts = params.get_string("ofn_counts");
  string ofn_totals = params.get_string("ofn_totals");
  
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

  // contig -> coord -> bin
  map<string, map < int, string > > bin_map;
  site_trj_read_segments(ifn_segments, bin_map);

  map< string, map< int, set< Variation > > > keys;
  read_sites(ifn_sites, keys);

  cout << "saving var matrix: " << ofn_counts << endl;
  ofstream out_counts(ofn_counts.c_str());
  massert(out_counts.is_open(), "could not open file %s", ofn_counts.c_str());

  cout << "saving total matrix: " << ofn_totals << endl;
  ofstream out_totals(ofn_totals.c_str());
  massert(out_totals.is_open(), "could not open file %s", ofn_totals.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////
  // print header
  ///////////////////////////////////////////////////////////////////////////////////////

  out_counts << "vid\tcontig\tcoord\tvariant\tbin";
  out_totals << "vid\tcontig\tcoord\tvariant\tbin";
  for (int i=0; i<nlibs; ++i) {
    out_counts << "\t" << "t" << i+1;
    out_totals << "\t" << "t" << i+1;
  }
  out_counts << endl;
  out_totals << endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // go over all sites
  ///////////////////////////////////////////////////////////////////////////////////////

  int index = 1;
  for (map< string, map< int, set< Variation > > >::iterator it=keys.begin(); it!=keys.end(); ++it) {
    string contig = (*it).first;
    map< int, set< Variation > >& keys_contig = (*it).second;
    for (map< int, set< Variation > >::iterator jt=keys_contig.begin(); jt != keys_contig.end(); ++jt) {
      int coord = (*jt).first;
      set< Variation >& keys_coord = (*jt).second;

      // get bin, if exists
      string bin = "NA";
      if (bin_map.find(contig) != bin_map.end() && bin_map[contig].find(coord) != bin_map[contig].end())
	bin = bin_map[contig][coord];
      
      for (set< Variation >::iterator kt=keys_coord.begin(); kt != keys_coord.end(); ++kt) {
	Variation var = (*kt);

	string vid = string("v") + std::to_string(index++);
	out_counts << vid << "\t" << contig << "\t" << coord+1 << "\t" << var << "\t" << bin;
	out_totals << vid << "\t" << contig << "\t" << coord+1 << "\t" << var << "\t" << bin;
	
	for (int i=0; i<nlibs; ++i) {
	  VariationSet& cav = cavs[i];

	  int count = cav.get_var_count(contig, coord, var);
	  int total = cav.get_coverage(contig, coord);

	  out_counts << "\t" << count;
	  out_totals << "\t" << total;
	}
	out_counts << endl;
	out_totals << endl;
      }
    }
  }
  
  return 0;
}
