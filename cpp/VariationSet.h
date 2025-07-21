#ifndef __VARIATIONSET__
#define __VARIATIONSET__

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include "Variation.h"

using namespace std;

typedef unsigned long ulong;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variation Set
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VariationSet {
 private:
  static const char* m_magicref;

  // contig -> coordinate -> variation -> count
  map< string, map< int, map <Variation, int> > > m_vars;

  // total coverage
  // contig -> coordinate -> count
  map<string, vector<int> > m_covs;

  // variants classified as seqeuncing noise
  // contig -> coordinate -> count
  map<string, vector<int> > m_covs_noise;
  
  // count reads that precisely start at coord
  map<string, vector<int> > m_covs_start;
  
  // count cumsum of reads that start at coord or left of it
  map<string, vector<int> > m_covs_cs;

  // keep max read, to warn if query segment is too short
  int m_read_count;
  int m_max_read_length;

  // track total length for efficient mean calculation
  long long m_total_read_length;
  double m_mean_read_length;

  void save_map_count_vector(boost::iostreams::filtering_ostream& out, map<string, vector<int> >& map_count);
  void load_map_count_vector(boost::iostreams::filtering_istream& in, map<string, vector<int> >& map_count);
  
  void save_cav(boost::iostreams::filtering_ostream& out);
  void load_cav(boost::iostreams::filtering_istream& in);

  void add_var_map(const map< string, map< int, map <Variation, int> > >& vars);

  // cumsum vector initialized manaully after adding all reads
  void init_cov_cs();
 public:
  VariationSet();
  VariationSet(const map<string, int>& contig_map);
  VariationSet(string fn);

  // direct functions for lazy coders
  map< string, map< int, map <Variation, int> > >& get_vars();
  map<string, vector<int> >& get_covs();
  map<string, vector<int> >& get_covs_cs();
  map<string, vector<int> >& get_covs_noise();
  
  // constant access
  const map< string, map< int, map <Variation, int> > >& get_vars() const;
  const map<string, vector<int> >& get_covs() const;
  
  // read stats
  int get_read_count();
  int get_max_read_length();
  double get_mean_read_length();
  
  // get contigs
  vector<string> get_contigs();

  // get number of variants
  void get_total_var_counts(int& count, int& unique_count);

  // get support for variant
  int get_contig_length(string contig) const;
  
  // get support for variant
  int get_var_count(string& contig, int& coord, Variation& var);

  // get support for ref sequence
  int get_ref_count(const string contig, const int coord) const;

  // get variations at coord
  void get_vars(const string contig, const int coord, vector<Variation>& vars);

  // get total coverage for coord
  int get_coverage(const string contig, const int coord);

  // get total coverage around coordinate
  int get_regional_coverage(const string contig, const int coord,
			    int regional_window, int contig_side_size);
  
  // get number of reads classified as having seq noise
  int get_coverage_noise(const string contig, const int coord);
  
  // get total coverage for segment
  int get_segment_coverage(const string contig, const int start, int end);

  // get total coverage of segment collection
  int get_segments_coverage(vector< Segment >& segs);
  
  // get mean/sd coverage of segment collection
  void get_mean_sd_segments_coverage(vector< Segment >& segs, double& cov_mean, double& cov_sd);
  
  // get major allele
  Variation get_major(const string contig, const int coord);
  
  // get allele counts
  void get_counts(const string contig, const int coord, vector<int>& counts);

  // get allele counts with variant identity
  void get_var_counts(const string contig, const int coord, bool ascending,
		      vector < pair < Variation, int > >& var_counts) const;
  
  // get all contig/coord/var keys
  void collect_var_keys(map< string, map< int, set< Variation > > >& keys);

  // get all contig/coord keys
  void collect_coord_keys(map< string, set< int > >& keys);

  // track read total count and max length
  void add_read(string contig, int coord, int length);
  
  // call after all read loaded, to compute mean read length
  void reads_done();
  
  // I/O
  void save(string fn);
  void load(string fn);

  friend VariationSet operator+(VariationSet const &, VariationSet const &);

  // for plotting, detailed
  void get_contig_data(string contig, int from, int to,
		       vector<string>& contig_r,
		       vector<int>& coord_r,
		       vector<string>& var_r,
		       vector<int>& count_r,
		       vector<int>& total_r,
		       vector<int>& cc_start_r,
		       vector<int>& cc_end_r);

  // for plotting, binning
  void get_summary(string contig, VariType vtype,
		   int& from, int& to, int min_cov, double max_seg_freq,
		   int& seg_count_r, int& median_cov_r);

  // for plotting, binning, only coverage
  void get_cov_summary(string contig, int& from, int& to, int min_cov, int& median_cov_r);

  // for plotting, binning, mean frequency
  void get_freq_summary(string contig, int& from, int& to,
			int& count_r, double& median_freq_r);
  
  // to reduce space, resrict to specified contigs
  void restrict(vector<string>& contigs);

  // restrict to sites with minimal contig length
  void restrict(map< string, map< int, set <Variation> > >& sites);

  // get mean and sd of coverage of segment
  void get_segment_stats(string contig, int start, int end,
			 double& mean_r, double& var_r);
};

#endif
