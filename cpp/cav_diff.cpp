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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void diff_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn1", new ParserFilename("first POP filename"), true);
  params.add_parser("ifn2", new ParserFilename("second POP filename"), true);

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

void diff_analyze(VariationSet& varset1, VariationSet& varset2)
{
  // check basic stats
  vector<string> contigs1 = varset1.get_contigs();
  vector<string> contigs2 = varset2.get_contigs();
  
  set<string> contigs1_set(contigs1.begin(), contigs1.end());
  set<string> contigs2_set(contigs2.begin(), contigs2.end());
  
  bool contigs_same = (contigs1_set == contigs2_set);
  
  // collect variant keys
  map< string, map< int, set< Variation > > > keys1, keys2;
  varset1.collect_var_keys(keys1);
  varset2.collect_var_keys(keys2);
  
  // count differences
  int variants_only_in_1 = 0;
  int variants_only_in_2 = 0;
  int variants_count_diff = 0;
  
  // check all variants from both files
  set<string> all_contigs;
  for (auto& p : keys1) all_contigs.insert(p.first);
  for (auto& p : keys2) all_contigs.insert(p.first);
  
  for (const string& contig : all_contigs) {
    bool has_contig1 = keys1.find(contig) != keys1.end();
    bool has_contig2 = keys2.find(contig) != keys2.end();
    
    set<int> all_coords;
    if (has_contig1) {
      for (auto& p : keys1[contig]) all_coords.insert(p.first);
    }
    if (has_contig2) {
      for (auto& p : keys2[contig]) all_coords.insert(p.first);  
    }

    for (int coord : all_coords) {
      bool has_coord1 = has_contig1 && keys1[contig].find(coord) != keys1[contig].end();
      bool has_coord2 = has_contig2 && keys2[contig].find(coord) != keys2[contig].end();
      
      set<Variation> all_vars;
      if (has_coord1) {
        for (auto& var : keys1[contig][coord]) all_vars.insert(var);
      }
      if (has_coord2) {
        for (auto& var : keys2[contig][coord]) all_vars.insert(var);
      }

      for (const Variation& var : all_vars) {
        int count1 = 0, count2 = 0;
        
        if (has_coord1 && keys1[contig][coord].find(var) != keys1[contig][coord].end()) {
          count1 = varset1.get_vars()[contig][coord][var];
        }
        if (has_coord2 && keys2[contig][coord].find(var) != keys2[contig][coord].end()) {
          count2 = varset2.get_vars()[contig][coord][var];
        }

        if (count1 > 0 && count2 == 0) variants_only_in_1++;
        else if (count1 == 0 && count2 > 0) variants_only_in_2++;
        else if (count1 != count2) variants_count_diff++;
      }
    }
  }
  
  // print brief summary
  if (contigs_same) {
    cout << "contigs: same" << endl;
  } else {
    cout << "contigs: different (" << contigs1.size() << " vs " << contigs2.size() << ")" << endl;
  }
  
  if (variants_only_in_1 == 0 && variants_only_in_2 == 0 && variants_count_diff == 0) {
    cout << "variants: identical" << endl;
  } else {
    cout << "variants: different";
    if (variants_only_in_1 > 0) cout << " (" << variants_only_in_1 << " only in file1)";
    if (variants_only_in_2 > 0) cout << " (" << variants_only_in_2 << " only in file2)";
    if (variants_count_diff > 0) cout << " (" << variants_count_diff << " count differences)";
    cout << endl;
  }
}

int diff_main(const char* name, int argc, char **argv)
{
  Parameters params;
  diff_init_params(name, argc, argv, params);

  string ifn1 = params.get_string("ifn1");
  string ifn2 = params.get_string("ifn2");

  VariationSet varset1;
  varset1.load(ifn1);

  VariationSet varset2;
  varset2.load(ifn2);

  diff_analyze(varset1, varset2);

  return 0;
} 