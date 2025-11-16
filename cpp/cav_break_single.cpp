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
// dissolve functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void save_round(string ofn, Dissolve& dissolve)
{
  vector < double > chivalues = dissolve.get_chivalues();
  vector < double > pvalues = dissolve.get_pvalues();
  vector < bool > is_center = dissolve.get_is_center();
  vector < double > weights = dissolve.get_weights();
  vector < double > margs = dissolve.get_marginals();

  massert(pvalues.size() == chivalues.size(), "size of vectors must match");

  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  out << "coord\tchivalue\tPvalue\tweight\tmarginal\tis_center" << endl;
  for (unsigned int i=0; i<pvalues.size(); ++i)
    out << i+1 << "\t" << chivalues[i] << "\t" << pvalues[i] << "\t" << weights[i] << "\t" << margs[i] << "\t" << (is_center[i] ? "T" : "F") << endl;

  out.close();
}

void save_counts(string ofn_prefix, string contig, vector < VariationSet >& covs)
{
  string ofn = ofn_prefix + ".mat";
  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  out << "coord";
  for (unsigned int i=0; i<covs.size(); ++i)
    out << "\t" << "l"+to_string(i+1);
  out << endl;

  int contig_length = covs[0].get_covs()[contig].size();
  for (int coord=0; coord<contig_length; ++coord) {
    out << coord;
    for (unsigned int i=0; i<covs.size(); ++i) {
      map<string, vector<int> >& set_cov = covs[i].get_covs();
      massert(set_cov.find(contig) != set_cov.end(), "contig not found");
      out << "\t" << set_cov[contig][coord];
    }
  out << endl;
  }

  out.close();
}

void dissolve_contig(string ofn_prefix, string contig, vector < VariationSet >& covs,
		     double outlier_fraction, double min_P, double pseudo_count, string weight_style)
{
  Dissolve dissolve(covs.size(), pseudo_count, weight_style);
  dissolve.set_counts(contig, covs);
  dissolve.init_once();

  for (unsigned i=0; i<1000; ++i) {
    string ofn = ofn_prefix + "." + to_string(i+1);
    save_round(ofn, dissolve);
    if (dissolve.reduce_center_round(outlier_fraction, min_P) == 0)
      break;
  }
  string ofn = ofn_prefix + ".final";
  save_round(ofn, dissolve);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void break_single_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("Table with multiple COV files"), true);
  params.add_parser("ofn_prefix", new ParserFilename("output file"), true);
  params.add_parser("contig", new ParserFilename("limit to a single contig", ""), false);
  params.add_parser("outlier_fraction", new ParserDouble("Fraction of outliers tested every round", 0.01), false);
  params.add_parser("p_value", new ParserDouble("Min Chi-Square P-value", 0.001), false);
  params.add_parser("pseudo_count", new ParserDouble("Add pseudo-count", 0.1), false);
  params.add_parser("max_lib_count", new ParserInteger("Use only n first libs (0 is all)", 0), false);
  params.add_parser("weight_style", new ParserString("Weight function of center positions (uniform|marginal)", "uniform"), false);

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

int break_single_main(const char* name, int argc, char **argv)
{
  Parameters params;
  break_single_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string contig = params.get_string("contig");
  double outlier_fraction = params.get_double("outlier_fraction");
  double min_P = params.get_double("p_value");
  double pseudo_count = params.get_double("pseudo_count");
  int max_lib_count = params.get_int("max_lib_count");
  string ofn_prefix = params.get_string("ofn_prefix");
  string weight_style = params.get_string("weight_style");

  vector< string > ifns;
  vector< string > library_ids;
  read_library_table(ifn, ifns, library_ids);

  int nlibs = ifns.size();
  if (max_lib_count > 0 && max_lib_count < nlibs)
    nlibs = max_lib_count;

  vector < VariationSet > cavs(nlibs);
  for (int i=0; i<nlibs; ++i) {
    VariationSet& cav = cavs[i];
    cov.load(ifns[i]);
  }

  save_counts(ofn_prefix, contig, covs);
  dissolve_contig(ofn_prefix, contig, covs, outlier_fraction, min_P, pseudo_count, weight_style);

  return 0;
}
