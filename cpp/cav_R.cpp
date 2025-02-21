#include <Rcpp.h>
#include "VariationSet.h"
#include "util.h"

using namespace Rcpp;
using namespace std;

// echo PKG_CPPFLAGS=-Wall -Wno-write-strings -std=c++0x > ~/.R/Makevars

/*
library("Rcpp")
Sys.setenv(PKG_CPPFLAGS="-Wall -Wno-write-strings -std=c++0x")
sourceCpp("cav_R.cpp", verbose=T)
cav = cav_load("/relman03/work/users/eitany/tempo/subjects/AAB/assembly/megahit/cav_v4/sets_N1/sets/s1/set.cav")
df = data.frame(contig="s1", from=1, to=10)
cav_query(cav, df)
cc = cav_contigs(cav)
*/

// [[Rcpp::export]]
Rcpp::XPtr<VariationSet> cav_load(string fn)
{
  VariationSet* varset = new VariationSet;
  varset->load(fn);
  Rcpp::XPtr<VariationSet> result(varset, true);
  return result;
}

// [[Rcpp::export]]
vector<string> cav_contigs(Rcpp::XPtr<VariationSet> xptr)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  return varset->get_contigs();
}

// [[Rcpp::export]]
DataFrame cav_query(Rcpp::XPtr<VariationSet> xptr, DataFrame df)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> coord_r;
  vector<string> var_r;
  vector<int> count_r;
  vector<int> total_r;
  vector<int> cc_start_r;
  vector<int> cc_end_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    varset->get_contig_data(contig, from, to, contig_r, coord_r, var_r, count_r, total_r, cc_start_r, cc_end_r);
  }

  // switch back to 1-based
  for (unsigned int i=0; i<coord_r.size(); ++i)
    coord_r[i]++;

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("coord") = coord_r,
				       Named("var") = var_r,
				       Named("count") = count_r,
				       Named("total") = total_r,
				       Named("cc_start") = cc_start_r,
				       Named("cc_end") = cc_end_r);

  return result;
}

/*
library("Rcpp")
Sys.setenv(PKG_CPPFLAGS="-Wall -Wno-write-strings -std=c++0x")
sourceCpp("cav_R.cpp", verbose=T)
cav = cav_load("/relman03/work/users/eitany/tempo/subjects/AAB/assembly/megahit/cav_v4/sets_N1/sets/s1/hosts.cav")
df = data.frame(contig="s1", from=1, to=1000)
cav_summary(cav, df, 2, 0.8)
 */

// [[Rcpp::export]]
DataFrame cav_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df, int min_cov, double max_seg_freq)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];
  StringVector vtype_v = df["vtype"];
  
  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> seg_r;
  vector<int> cov_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    VariType vtype = Variation::string_to_type(as<string>(vtype_v[i]));
  
    int seg, cov;
    varset->get_summary(contig, vtype, from, to, min_cov, max_seg_freq, seg, cov);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    seg_r.push_back(seg);
    cov_r.push_back(cov);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("seg.count") = seg_r,
				       Named("xcov") = cov_r);
  return result;
}

// [[Rcpp::export]]
DataFrame cav_cov_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df, int min_cov)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> cov_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    int cov;
    varset->get_cov_summary(contig, from, to, min_cov, cov);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    cov_r.push_back(cov);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("xcov") = cov_r);
  return result;
}

// [[Rcpp::export]]
DataFrame cav_freq_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> count_r;
  vector<double> freq_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    int count;
    double freq;
    varset->get_freq_summary(contig, from, to, count, freq);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    count_r.push_back(count);
    freq_r.push_back(freq);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("count") = count_r,
				       Named("freq") = freq_r);
  return result;
}

// [[Rcpp::export]]
int cav_regional_cov(Rcpp::XPtr<VariationSet> xptr,
			   string contig, int coord, int regional_window, int contig_side_size)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  return varset->get_regional_coverage(contig, coord, regional_window, contig_side_size);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// annotate sites
////////////////////////////////////////////////////////////////////////////////////////////////////

bool sites_is_internal(string contig, int coord, int margin,
		       Rcpp::XPtr<map< string, set< Segment > > > segs_p)
{
  map< string, set< Segment > >& segs = *segs_p;
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

// [[Rcpp::export]]
Rcpp::XPtr<map< string, set< Segment > > > cav_build_segments(DataFrame df)
{
  StringVector segment_v = df["segment"];
  StringVector contig_v = df["contig"];
  IntegerVector start_v = df["start"];
  IntegerVector end_v = df["end"];

  map< string, set< Segment > >* segs = new map< string, set< Segment > >;
  for (unsigned int i=0; i<segment_v.size(); ++i) {
    string seg_id = as<string>(segment_v[i]);
    string contig = as<string>(contig_v[i]);
    int start = start_v[i] - 1;
    int end = end_v[i] - 1;
    
    Segment seg(seg_id, contig, start, end);
    (*segs)[contig].insert(seg);
  }

  Rcpp::XPtr<map< string, set< Segment > > > result(segs, true);
  return result;
}

// [[Rcpp::export]]
DataFrame cav_annotate_sites(Rcpp::XPtr<VariationSet> varset,
			     List vsets,
			     double pseudo_count,
			     int regional_window, int margin,
			     Rcpp::XPtr<map< string, set< Segment > > > segs_p,
			     DataFrame df)
{
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> coord_r;
  vector<string> var_r;
  vector<int> count_r;
  vector<int> total_r;
  vector<int> n_samples_r;
  vector<bool> is_internal_r;
  vector<double> variant_p_r;
  vector<double> regional_p_r;
  vector<double> comp_p_r;

  map< string, map< int, map <Variation, int> > >& vars = varset->get_vars();
  
  int N = vsets.length();
  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    // int from = from_v[i] - 1;
    // int to = to_v[i] - 1;

    if (vars.find(contig) == vars.end())
      continue;
    
    map< int, map <Variation, int> >& vars_contig = vars[contig];
    for (map< int, map <Variation, int> >::iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      map <Variation, int> vars_coord = (*jt).second;

      // get total coverage
      int coord_cov = varset->get_coverage(contig, coord);

      // correspondance between local and regional total coverage vectors
      vector<double> loc_v(N);
      vector<double> reg_v(N);
      for (int i=0; i<N; ++i) {
	Rcpp::XPtr<VariationSet> v = vsets[i];	
	loc_v[i] = v->get_coverage(contig, coord) + pseudo_count;
	reg_v[i] = v->get_regional_coverage(contig, coord, regional_window, margin) + pseudo_count;
      }
      double p_reg = chi_square_test_p(loc_v, reg_v);

      // margin from segment sides, as defined in the input
      // int is_internal = coord >= (from + margin) && coord <= (to - margin);
      int is_internal = sites_is_internal(contig, coord, margin, segs_p);
      
      for (map <Variation, int>::iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int var_cov = varset->get_var_count(contig, coord, var);

	// count number of unique samples
	int n_samples = 0;
	
	// variant vector
	vector<double> var_v(N);
	vector<double> comp_v(N);
	for (int i=0; i<N; ++i) {
	  Rcpp::XPtr<VariationSet> v = vsets[i];	
	  int v_count = v->get_var_count(contig, coord, var);
	  if (v_count > 0)
	    n_samples++;
	  var_v[i] = v_count + pseudo_count;
	  comp_v[i] = loc_v[i] - var_v[i] + pseudo_count;
	}
	// compare variant and local total coverage 
	double p_var = chi_square_test_p(loc_v, var_v);

	// compare complement of variant with regional
	double p_comp = chi_square_test_p(comp_v, reg_v);

	// output
	contig_r.push_back(contig);
	coord_r.push_back(coord+1);
	var_r.push_back(var.to_string());
	count_r.push_back(var_cov);
	total_r.push_back(coord_cov);
	n_samples_r.push_back(n_samples);
	is_internal_r.push_back(is_internal);
	variant_p_r.push_back(p_var);
	regional_p_r.push_back(p_reg);
	comp_p_r.push_back(p_comp);
      }
    }
  }

  DataFrame result = DataFrame::create(Named("contig") = contig_r,
				       Named("coord") = coord_r,
				       Named("var") = var_r,
				       Named("count") = count_r,
				       Named("total") = total_r,
				       Named("n_samples") = n_samples_r,
				       Named("is_internal") = is_internal_r,
				       Named("variant_p") = variant_p_r,
				       Named("regional_p") = regional_p_r,
				       Named("comp_p") = comp_p_r);
  
  return result;
}

// [[Rcpp::export]]
Rcpp::XPtr<VariationSet> cav_add(List vsets)
{
  VariationSet* vset = new VariationSet;
  Rcpp::XPtr<VariationSet> v1 = vsets[0];

  *vset = *v1;
  int N = vsets.length();
  for (int i=1; i<N; ++i) {
    Rcpp::XPtr<VariationSet> v = vsets[i];	
    *vset = *vset + *v;
  }
  Rcpp::XPtr<VariationSet> result(vset, true);
  return result;
}

