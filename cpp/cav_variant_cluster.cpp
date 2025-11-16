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

#include <unistd.h>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

#include "ClusterVariants.h"
#include "Resource.h"

typedef map< string, set< VariationKey > >::iterator VkIterator;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// util functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void variant_cluster_read_segments(string fn,
				   vector <string>& bins,
				   map<string, map < int, int > >& bin_map)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  map<string, int> bin_to_index;
  
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

    // add if needed
    if (bin_to_index.find(bin) == bin_to_index.end()) {
      bin_to_index[bin] = bin_to_index.size();
      bins.push_back(bin);
    }
    
    int bin_index = bin_to_index[bin]; 
    for (int coord=start; coord<=end; ++coord)
      bin_map[contig][coord] = bin_index;
  }
}

void variant_cluster_read_sites(string fn,
				vector <string>& bins,
				map<string, map < int, int > >& bin_map,
				map< string, set < VariationKey > >& keys)
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

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int coord = safe_string_to_int(fields[coord_ind], "coordinate") - 1;
    Variation var = Variation(fields[var_ind]);

    if (bin_map.find(contig) == bin_map.end() || bin_map[contig].find(coord) == bin_map[contig].end())
      continue;
    int bin_index = bin_map[contig][coord];
    string bin = bins[bin_index];
    
    VariationKey vk(id, contig, coord, var);
    keys[bin].insert(vk);
  }
}

void cluster_variants_single(string bin,
			     set < VariationKey >& bin_keys,
			     map< string, vector < ClusterVariantItem > >& result,
			     vector < VariationSet >& cavs,
			     int nlibs,
			     double p_threshold, double pseudo_count,
			     int thread_count)
{
  // prepare vector of ClusterVariantItem
  vector < ClusterVariantItem > cvi_vec;
  for (set < VariationKey >::iterator jt=bin_keys.begin(); jt != bin_keys.end(); ++jt) {
    VariationKey vk = (*jt);
    ClusterVariantItem cvi(vk.id, vk.contig, vk.coord, vk.var);
    cvi.counts.resize(cavs.size());
    for (unsigned int k=0; k<cavs.size(); ++k) {
      VariationSet& cav = cavs[k];
      int count = cav.get_var_count(cvi.contig, cvi.coord, cvi.var);
      cvi.counts[k] = count + pseudo_count;
    }
    cvi_vec.push_back(cvi);
  }
  
  vector<int> vclusters(cvi_vec.size());
  ClusterVariants cv(bin, p_threshold, cvi_vec, nlibs, thread_count);
  cv.cluster_variants(vclusters);
  for (unsigned int i=0; i<cvi_vec.size(); ++i) {
    cvi_vec[i].counts.resize(0);
    cvi_vec[i].cluster_index = vclusters[i]+1;
  }
  
  result[bin] = cvi_vec;
}

void cluster_variants(map< string, set < VariationKey > >& keys,
		      map< string, vector < ClusterVariantItem > >& result,
		      vector < VariationSet >& cavs,
		      int nlibs,
		      double p_threshold, double pseudo_count,
		      int thread_count)
{
  for (VkIterator it=keys.begin(); it!=keys.end(); ++it) {
    string bin = (*it).first;
    set < VariationKey >& bin_keys = (*it).second;

    cluster_variants_single(bin, bin_keys, result, cavs, nlibs,
			    p_threshold, pseudo_count, thread_count);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void variant_cluster_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segments", new ParserFilename("input segment file, with bin field"), true);
  params.add_parser("ifn_sites", new ParserFilename("input table with sites"), true);
  params.add_parser("threads", new ParserInteger("Number of threads", 10), true);
  params.add_parser("p_value", new ParserDouble("Maximal p-value for Chi-Square test", 0.01), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("ofn", new ParserFilename("output table"), true);
  
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

int variant_cluster_main(const char* name, int argc, char **argv)
{
  Parameters params;
  variant_cluster_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_segments = params.get_string("ifn_segments");
  string ifn_sites = params.get_string("ifn_sites");
  double p_threshold = params.get_double("p_value");
  double pseudo_count = params.get_double("pseudo_count");
  string ofn = params.get_string("ofn");
  int thread_count = params.get_int("threads");

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

  // vector of bins
  vector <string> bins;
  
  // contig -> coord -> bin_index
  map<string, map < int, int > > bin_map;
  
  variant_cluster_read_segments(ifn_segments, bins, bin_map);

  // bin -> set of VariantionKey
  map< string, set < VariationKey > > keys;
  variant_cluster_read_sites(ifn_sites, bins, bin_map, keys);

  // result
  map< string, vector < ClusterVariantItem > > vars;

  cluster_variants(keys, vars, cavs, nlibs, p_threshold, pseudo_count, thread_count);

  //////////////////////////////////////////////////////////////////////////////////////////
  // save results
  //////////////////////////////////////////////////////////////////////////////////////////

  cout << "saving table: " << ofn << endl;
  ofstream out(ofn.c_str());
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  out << "vid\tcontig\tcoord\tvariant\tbin\tvcluster" << endl;

  for (map< string, vector < ClusterVariantItem > >::iterator it=vars.begin(); it != vars.end(); ++it) {
    string bin = (*it).first;
    vector < ClusterVariantItem >& bin_vars = (*it).second;
    
    // cout << "bin: " << bin << ", n=" << nclusters << endl;
    for (unsigned int i=0; i<bin_vars.size(); ++i) {
      ClusterVariantItem& cvi = bin_vars[i];
      out << cvi.id << "\t" << cvi.contig << "\t" << cvi.coord+1 << "\t" << cvi.var << "\t";
      out << bin << "\t" << bin << "_v" << cvi.cluster_index << endl;
    }
  }
  
  return 0;
}
