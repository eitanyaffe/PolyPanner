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

void vcluster_read_sites(string fn,
			 map<string, string>& vcluster2bin,
			 map< string, set< VariationKey > >& keys)
{
  cout << "reading site table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int id_ind = get_field_index("vid", fields);
  int contig_ind = get_field_index("contig", fields);
  int coord_ind = get_field_index("coord", fields);
  int var_ind = get_field_index("variant", fields);
  int bin_ind = get_field_index("bin", fields);
  int vcluster_ind = get_field_index("vcluster", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int coord = safe_string_to_int(fields[coord_ind], "coordinate") - 1;
    Variation var = Variation(fields[var_ind]);
    string vcluster = fields[vcluster_ind];
    string bin = fields[bin_ind];
    
    VariationKey vk(id, contig, coord, var);
    keys[vcluster].insert(vk);
    vcluster2bin[vcluster] = bin;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void vcluster_trajectory_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_sites", new ParserFilename("input table with sites"), true);
  params.add_parser("ofn_counts", new ParserFilename("output matrix with raw counts"), true);
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

int vcluster_trajectory_main(const char* name, int argc, char **argv)
{
  Parameters params;
  vcluster_trajectory_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
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

  map< string, set< VariationKey > > keys;
  map<string, string> vcluster2bin;
  vcluster_read_sites(ifn_sites, vcluster2bin, keys);

  cout << "saving var count matrix: " << ofn_counts << endl;
  ofstream out_counts(ofn_counts.c_str());
  massert(out_counts.is_open(), "could not open file %s", ofn_counts.c_str());

  cout << "saving total matrix: " << ofn_totals << endl;
  ofstream out_totals(ofn_totals.c_str());
  massert(out_totals.is_open(), "could not open file %s", ofn_totals.c_str());

  ///////////////////////////////////////////////////////////////////////////////////////
  // print header
  ///////////////////////////////////////////////////////////////////////////////////////

  out_counts << "vcluster\tbin";
  out_totals << "vcluster\tbin";
  for (int i=0; i<nlibs; ++i) {
    out_counts << "\t" << "t" << i+1;
    out_totals << "\t" << "t" << i+1;
  }
  out_counts << endl;
  out_totals << endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // go over all sites
  ///////////////////////////////////////////////////////////////////////////////////////

  for (map< string, set< VariationKey > >::iterator it=keys.begin(); it!=keys.end(); ++it) {
    string vcluster = (*it).first;
    set < VariationKey >& vks = (*it).second;
    string bin = vcluster2bin[vcluster];
    
    // aggregate counts across vcluster
    vector<int> counts(nlibs);
    vector<int> totals(nlibs);
    for (int i=0; i<nlibs; ++i) {
      counts[i] = 0;
      totals[i] = 0;
    }
    for (set < VariationKey >::iterator jt=vks.begin(); jt != vks.end(); ++jt) {
      VariationKey vk = (*jt);
      string contig = vk.contig;
      int coord = vk.coord;
      Variation var = vk.var;

      for (int i=0; i<nlibs; ++i) {
	VariationSet& cav = cavs[i];
	
	int count = cav.get_var_count(contig, coord, var);
	int total = cav.get_coverage(contig, coord);
	counts[i] += count;
	totals[i] += total;
      }
    }
    
    // output results
    out_counts << vcluster << "\t" << bin;
    out_totals << vcluster << "\t" << bin;
    for (int i=0; i<nlibs; ++i) {
      out_counts << "\t" << counts[i];
      out_totals << "\t" << totals[i];
    }
    out_counts << endl;
    out_totals << endl;
  }
  
  return 0;
}
