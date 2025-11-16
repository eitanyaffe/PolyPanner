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
#include "Filter.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void cov_mat_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_cavs", new ParserFilename("Table with multiple POP files"), true);
  params.add_parser("ifn_segs", new ParserFilename("Segment table"), true);
  params.add_parser("ifn_fasta", new ParserFilename("Contig fasta"), true);
  params.add_parser("actual_nts", new ParserBoolean("Actual nts or random nts used to avoid TNF usage)", true), true);
  params.add_parser("ofn_fasta", new ParserFilename("Output fasta file"), true);
  params.add_parser("ofn_mat", new ParserFilename("Output matrix with segments, bp-coverage and bp-variance"), true);
  params.add_parser("ofn_lib_map", new ParserFilename("Output library index to library ID mapping file"), false);
  params.add_parser("min_segment_length", new ParserInteger("Minimum segment length", 1000), false);

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

int cov_mat_main(const char* name, int argc, char **argv)
{
  Parameters params;
  cov_mat_init_params(name, argc, argv, params);

  string ifn_cavs = params.get_string("ifn_cavs");
  string ifn_segs = params.get_string("ifn_segs");
  string ifn_fasta = params.get_string("ifn_fasta");
  bool use_random_nts = !params.get_bool("actual_nts");
  string ofn_fa = params.get_string("ofn_fasta");
  string ofn_mat = params.get_string("ofn_mat");
  string ofn_lib_map = params.get_string("ofn_lib_map");
  int min_segment_length = params.get_int("min_segment_length");

  map<string, string> fasta;
  load_fasta(ifn_fasta, fasta);

  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn_cavs, ifns, library_ids);

  int nlibs = ifns.size();
  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cav.load(ifns[i]);
  }

  cout << "reading table: " << ifn_segs << endl;
  ifstream in(ifn_segs.c_str());
  massert(in.is_open(), "could not open file %s", ifn_segs.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int seg_ind = get_field_index("segment", fields);
  int contig_ind = get_field_index("contig", fields);
  int start_ind = get_field_index("start", fields);
  int end_ind = get_field_index("end", fields);

  // coverage matrix
  ofstream mat_out(ofn_mat.c_str(), ios::out);
  massert(mat_out.is_open(), "could not open file %s", ofn_mat.c_str());
  mat_out << "segment";
  for (int i=0; i<nlibs; ++i)
    mat_out << "\tcov_" << i+1 << "\t" << "var_" << i+1;
  mat_out << endl;

  // segment fasta
  ofstream fa_out(ofn_fa.c_str(), ios::out);
  massert(fa_out.is_open(), "could not open file %s", ofn_fa.c_str());
  
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string segment = fields[seg_ind];
    string contig = fields[contig_ind];
    int start = atoi(fields[start_ind].c_str())-1;
    int end = atoi(fields[end_ind].c_str())-1;
    int len = end-start+1;

    if (len < min_segment_length)
      continue;

    // get segment nts
    string nts;
    if (use_random_nts) {
      nts = rand_nts(len);
    } else {
      massert(fasta.find(contig) != fasta.end(), "contig %s not found in fasta file", contig.c_str());
      string cfasta = fasta[contig];
      massert(start >= 0 && end < (int)cfasta.length(), "segment coords out of range");
      nts = cfasta.substr(start, len);
    }
      
    fa_out << ">" << segment << endl;
    fa_out << nts << endl;
    
    // get coverage
    mat_out << segment;
    for (int i=0; i<nlibs; ++i) {
      double mean, var;
      cavs[i].get_segment_stats(contig, start, end, mean, var);
      mat_out << "\t" << mean  << "\t" << var;
    }
    mat_out << endl;
  }
  
  // write library mapping file if requested
  if (!ofn_lib_map.empty()) {
    cout << "writing library mapping to " << ofn_lib_map << endl;
    ofstream lib_map_out(ofn_lib_map.c_str());
    massert(lib_map_out.is_open(), "could not open file %s", ofn_lib_map.c_str());
    
    lib_map_out << "lib_index\tlib_id" << endl;
    for (int i = 0; i < nlibs; ++i) {
      lib_map_out << "lib_" << (i+1) << "\t" << library_ids[i] << endl;
    }
    
    lib_map_out.close();
  }
  
  return 0;
}
