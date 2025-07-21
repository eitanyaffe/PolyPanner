#include <cstring>
#include <numeric>

#include "VariationSet.h"
#include "util.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variation Set
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

VariationSet::VariationSet() : m_read_count(0), m_max_read_length(0), m_total_read_length(0) {}

VariationSet::VariationSet(const map<string, int>& contig_map) : VariationSet()
{
  for (map<string, int>::const_iterator it=contig_map.begin(); it!=contig_map.end(); ++it) {
    string contig = (*it).first;
    int length = (*it).second;

    vector<int>& coverage_contig = m_covs[contig];
    coverage_contig.resize(length);

    vector<int>& coverage_contig_start = m_covs_start[contig];
    coverage_contig_start.resize(length);

    vector<int>& coverage_contig_cs = m_covs_cs[contig];
    coverage_contig_cs.resize(length);
    
    vector<int>& noise_contig = m_covs_noise[contig];
    noise_contig.resize(length);
  }
}

VariationSet::VariationSet(string fn)  : VariationSet()
{
  load(fn);
}

map< string, map< int, map <Variation, int> > >& VariationSet::get_vars()
{
  return m_vars;
}

const map< string, map< int, map <Variation, int> > >& VariationSet::get_vars() const
{
  return m_vars;
}

map<string, vector<int> >& VariationSet::get_covs()
{
  return m_covs;
}

const map<string, vector<int> >& VariationSet::get_covs() const
{
  return m_covs;
}

map<string, vector<int> >& VariationSet::get_covs_cs()
{
  return m_covs_cs;
}

map<string, vector<int> >& VariationSet::get_covs_noise()
{
  return m_covs_noise;
}

int VariationSet::get_read_count()
{
  return m_read_count;
}

int VariationSet::get_max_read_length()
{
  return m_max_read_length;
}

double VariationSet::get_mean_read_length()
{
  return m_mean_read_length;
}

vector<string> VariationSet::get_contigs()
{
  vector<string> result;
  for (auto const& element : m_covs)
    result.push_back(element.first);
  return result;
}

void VariationSet::collect_var_keys(map< string, map< int, set< Variation > > >& keys)
{
  map< string, map< int, map <Variation, int> > >& vars = get_vars();
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    const map< int, map <Variation, int> >& vars_contig = (*it).second;
    map< int, set <Variation> >& result_vars_contig = keys[contig];

    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      const map <Variation, int>& vars_coord = (*jt).second;
      set <Variation >& result_vars_coord = result_vars_contig[coord];
      for (map <Variation, int>::const_iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	result_vars_coord.insert(var);
      }
    }
  }
}

void VariationSet::collect_coord_keys(map< string, set< int > >& keys)
{
  map< string, map< int, map <Variation, int> > >& vars = get_vars();
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    const map< int, map <Variation, int> >& vars_contig = (*it).second;
    set <int>& keys_contig = keys[contig];

    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      keys_contig.insert(coord);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Save Variation Set
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VariationSet::save_cav(boost::iostreams::filtering_ostream& out)
{
  // write number of contigs
  size_t size = m_vars.size();
  out.write(reinterpret_cast<const char*>(&size), sizeof(size));

  for (map< string, map< int, map <Variation, int> > >::iterator it=m_vars.begin(); it!=m_vars.end(); ++it) {
    string contig = (*it).first;
    map< int, map <Variation, int> >& vars_contig = (*it).second;

    // write contig id length
    size_t size_contig = contig.size();
    out.write(reinterpret_cast<const char*>(&size_contig), sizeof(size_contig));

    // write contig id
    out.write(contig.c_str(), sizeof(char)*contig.size());

    // write number of coords per contig
    size_t size_vars_contig = vars_contig.size();
    out.write(reinterpret_cast<const char*>(&size_vars_contig), sizeof(size_vars_contig));

    for (map< int, map <Variation, int> >::iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      map <Variation, int>& vars_coord = (*jt).second;

      // write coord
      out.write(reinterpret_cast<const char*>(&coord), sizeof(int));

      // write number of vars per coord
      size_t size_vars_coord = vars_coord.size();
      out.write(reinterpret_cast<const char*>(&size_vars_coord), sizeof(size_vars_coord));

      for (map <Variation, int>::iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;

	var.save(out);
	out.write(reinterpret_cast<const char*>(&count), sizeof(int));
      }
    }
  }
}

void VariationSet::load_cav(boost::iostreams::filtering_istream& in)
{
  // load number of contigs
  size_t n_contigs;
  in.read(reinterpret_cast<char*>(&n_contigs), sizeof(n_contigs));

  for (unsigned i=0; i<n_contigs; ++i) {
    string contig;

    // read contig id size
    size_t size_contig;
    in.read(reinterpret_cast<char*>(&size_contig), sizeof(size_contig));

    // read contig id
    contig.resize(size_contig);
    in.read(&contig[0], sizeof(char)*contig.size());

    map< int, map <Variation, int> >& vars_contig = m_vars[contig];

    // read number of coords per contig
    size_t n_coords;
    in.read(reinterpret_cast<char*>(&n_coords), sizeof(n_coords));

    for (unsigned j=0; j<n_coords; ++j) {

      // read coord
      int coord;
      in.read(reinterpret_cast<char*>(&coord), sizeof(int));

      map <Variation, int>& vars_coord = vars_contig[coord];

      // read number of vars per coord
      size_t n_vars;
      in.read(reinterpret_cast<char*>(&n_vars), sizeof(n_vars));

      for (unsigned k=0; k<n_vars; ++k) {
	Variation var;
	int count;

	var.load(in);
	in.read(reinterpret_cast<char*>(&count), sizeof(int));

	vars_coord[var] = count;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// coverage load/save (VariationSet)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void VariationSet::save_map_count_vector(boost::iostreams::filtering_ostream& out,
					 map<string, vector<int> >& map_count)
{  
  size_t size = map_count.size();
  out.write(reinterpret_cast<const char*>(&size), sizeof(size));
  
  for (map<string, vector<int> >::iterator it=map_count.begin(); it!=map_count.end(); ++it) {
    string contig = (*it).first;
    vector<int>& coverage_contig = (*it).second;

    // write contig id length
    size_t size_contig = contig.size();
    out.write(reinterpret_cast<const char*>(&size_contig), sizeof(size_contig));

    // write contig id
    out.write(contig.c_str(), sizeof(char)*contig.size());

    // write coverage vector length
    size_t length = coverage_contig.size();
    out.write(reinterpret_cast<const char*>(&length), sizeof(length));

    // write coverage vector values
    out.write(reinterpret_cast<const char*>(&coverage_contig[0]), coverage_contig.size()*sizeof(int));
  }
}

void VariationSet::load_map_count_vector(boost::iostreams::filtering_istream& in,
					 map<string, vector<int> >& map_count)
{
  // cout << "loading coverage vectors" << endl;
  // read number of contigs
  size_t n_contigs;
  in.read(reinterpret_cast<char*>(&n_contigs), sizeof(n_contigs));

  for (unsigned i=0; i<n_contigs; ++i) {
    // read contig
    string contig;

    size_t size_contig;
    in.read(reinterpret_cast<char*>(&size_contig), sizeof(size_contig));

    contig.resize(size_contig);
    in.read(&contig[0], sizeof(char)*contig.size());

    vector<int>& coverage_contig = map_count[contig];

    // read coverage vector length
    size_t length;
    in.read(reinterpret_cast<char*>(&length), sizeof(length));

    // read coverage vector values
    coverage_contig.resize(length);
    in.read(reinterpret_cast<char*>(&coverage_contig[0]), coverage_contig.size()*sizeof(int));
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Load/Save Variation Set
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static const char VariationSet_magicref[] = {0x01,0x03,0x02,0x04};
const char* VariationSet::m_magicref = VariationSet_magicref;

void VariationSet::save(string fn)
{
  ofstream file(fn.c_str(), ios::out | ios::binary);
  massert(file.is_open(), "could not open file %s", fn.c_str());

  // adding compression filter if needed
  boost::iostreams::filtering_ostream out;
  cout << "saving POP file: " << fn << endl;
  out.push(boost::iostreams::gzip_compressor());
  out.push(file);
  
  // save magic number
  out.write(m_magicref, 4);

  out.write(reinterpret_cast<const char*>(&m_max_read_length), sizeof(m_max_read_length));
  out.write(reinterpret_cast<const char*>(&m_read_count), sizeof(m_read_count));
  out.write(reinterpret_cast<const char*>(&m_mean_read_length), sizeof(m_mean_read_length));
  
  save_map_count_vector(out, m_covs);
  save_map_count_vector(out, m_covs_cs);
  save_map_count_vector(out, m_covs_noise);
  save_cav(out);
  
  boost::iostreams::close(out);
  file.close();
}

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

void VariationSet::load(string fn)
{
  ifstream file(fn.c_str());
  massert(file.is_open(), "could not open file %s", fn.c_str());

  // adding decompression filter if needed
  boost::iostreams::filtering_istream in;
  cout << "reading POP file: " << fn << endl;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(file);
  
  char magic[4];
  in.read(magic, 4);
  massert((memcmp(magic, m_magicref, sizeof(magic)) == 0), "magic number not found, check file format");

  in.read(reinterpret_cast<char*>(&m_max_read_length), sizeof(m_max_read_length));
  in.read(reinterpret_cast<char*>(&m_read_count), sizeof(m_read_count));
  in.read(reinterpret_cast<char*>(&m_mean_read_length), sizeof(m_mean_read_length));
  
  load_map_count_vector(in, m_covs);
  load_map_count_vector(in, m_covs_cs);
  load_map_count_vector(in, m_covs_noise);
  load_cav(in);

  boost::iostreams::close(in);
  file.close();
}

void VariationSet::add_var_map(const map< string, map< int, map <Variation, int> > >& vars)
{
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    const map< int, map <Variation, int> >& vars_contig = (*it).second;
    map< int, map <Variation, int> >& this_vars_contig = m_vars[contig];

    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      const map <Variation, int>& vars_coord = (*jt).second;
      map <Variation, int>& this_vars_coord = this_vars_contig[coord];
      for (map <Variation, int>::const_iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;
	this_vars_coord[var] += count;
      }
    }
  }
}

void add_contig_vectors(set<string>& contigs,
			const map< string, vector<int> >& cmap1,
			const map< string, vector<int> >& cmap2,
			map< string, vector<int> >& cmap_r)
{
  set<string> contigs1, contigs2;
  for(map<string, vector<int> >::const_iterator it=cmap1.begin(); it != cmap1.end(); ++it)
    contigs1.insert((*it).first);
  for(map<string, vector<int> >::const_iterator it=cmap2.begin(); it != cmap2.end(); ++it)
    contigs2.insert((*it).first);
  massert(contigs == contigs1 && contigs == contigs2, "contig sets must be idnetical");

  for (set<string>::iterator it = contigs.begin(); it != contigs.end(); ++it) {
    string contig = (*it);

    massert(cmap1.find(contig) != cmap1.end(), "contig %s not found in contig map", contig.c_str());
    massert(cmap2.find(contig) != cmap2.end(), "contig %s not found in contig map", contig.c_str());
  
    const vector<int>& vec1 = cmap1.at(contig);
    const vector<int>& vec2 = cmap2.at(contig);
    massert(vec1.size() == vec2.size(),
	    "contig vectors must be equal size for contig %s", contig.c_str());
    
    vector<int>& vec_r = cmap_r[contig];
    vec_r.resize(vec1.size());
    for (unsigned i=0; i<vec1.size(); ++i)
      vec_r[i] = vec1[i] + vec2[i];
  }
}

VariationSet operator+(VariationSet const &v1, VariationSet const &v2)
{
  VariationSet result;

  // add cav maps
  result.add_var_map(v1.m_vars);
  result.add_var_map(v2.m_vars);

  result.m_read_count = v1.m_read_count + v2.m_read_count;
  result.m_max_read_length = max(v1.m_max_read_length, v2.m_max_read_length);

  set<string> contigs;
  for(map<string, vector<int> >::const_iterator it=v1.m_covs.begin(); it != v1.m_covs.end(); ++it)
    contigs.insert(it->first);

  add_contig_vectors(contigs, v1.m_covs, v2.m_covs, result.m_covs);
  add_contig_vectors(contigs, v1.m_covs_cs, v2.m_covs_cs, result.m_covs_cs);
  add_contig_vectors(contigs, v1.m_covs_noise, v2.m_covs_noise, result.m_covs_noise);

  return result;
}

int VariationSet::get_contig_length(string contig) const
{
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  return m_covs.at(contig).size();
}

int VariationSet::get_var_count(string& contig, int& coord, Variation& var)
{
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  if (var.is_ref())
    return get_ref_count(contig, coord);

  // contig not found
  if (m_vars.find(contig) == m_vars.end())
    return 0;
  map< int, map <Variation, int> >& contig_vars = m_vars[contig];

  // coord not found
  if (contig_vars.find(coord) == contig_vars.end())
    return 0;
  map <Variation, int>& coord_vars = contig_vars[coord];

  // var not found
  if (coord_vars.find(var) == coord_vars.end())
    return 0;

  return coord_vars[var];
}

inline void assert_vcoord(const vector<int>& v, int coord)
{
  massert(coord >= 0 && coord < (int)v.size(),
	  "coordinate %d out of range, length of vector: %d", coord, v.size());
}

int VariationSet::get_coverage(const string contig, const int coord)
{
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  vector<int>& cov = m_covs[contig];
  assert_vcoord(cov, coord);
  return cov[coord];
}

int VariationSet::get_regional_coverage(const string contig, const int coord,
					int regional_window, int contig_side_size)
{
  if (m_covs.find(contig) == m_covs.end()) {
    cout << "warning: contig not found: " << contig << endl;
    return -1;
  }
  int clen = m_covs.at(contig).size();
  int start = max(contig_side_size, coord-regional_window/2);
  int end = min(clen-contig_side_size, coord+regional_window/2);
  if (end>start)
    return get_segment_coverage(contig, start, end);
  else
    return -1;
} 

int VariationSet::get_coverage_noise(const string contig, const int coord)
{
  massert(m_covs_noise.find(contig) != m_covs_noise.end(), "contig not found");
  vector<int>& cov = m_covs_noise[contig];
  assert_vcoord(cov, coord);
  return cov[coord];
}

// call after all read loaded, to compute mean read length
void VariationSet::reads_done()
{
  if (m_read_count > 0) {
    m_mean_read_length = (double)m_total_read_length / m_read_count;
  } else {
    m_mean_read_length = 0.0;
  }
  init_cov_cs();
}

void VariationSet::add_read(string contig, int coord, int length)
{
  m_total_read_length += length;
  m_read_count++;
  if (length > m_max_read_length)
    m_max_read_length = length;

  massert(m_covs_start.find(contig) != m_covs_start.end(), "contig not found");
  vector<int>& cov_start = m_covs_start[contig];
  assert_vcoord(cov_start, coord);
  cov_start[coord]++;
}

void VariationSet::init_cov_cs()
{
  for (map<string, vector<int> >::iterator it=m_covs_cs.begin(); it != m_covs_cs.end(); ++it) {
    string contig = (*it).first;
    vector<int>& cov_cs = (*it).second;
    
    massert(m_covs_start.find(contig) != m_covs_start.end(), "contig not found");
    vector<int>& cov_start = m_covs_start[contig];
    int cs = 0;
    for (unsigned int coord=0; coord<cov_cs.size(); ++coord) {
      int count = cov_start[coord];
      cs += count;
      cov_cs[coord] = cs;
    }
  }
}

int VariationSet::get_segment_coverage(const string contig, const int start, int end)
{
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  vector<int>& cov = m_covs[contig];

  massert(m_covs_cs.find(contig) != m_covs_cs.end(), "contig not found");
  vector<int>& cov_cs = m_covs_cs[contig];

  int len = end-start+1;
  massert(end > start && len >= m_max_read_length,
	  "query segment length %d is larger than max read length %d, start=%d, end=%d",
	  len, m_max_read_length, start, end);
  
  assert_vcoord(cov, end);
  assert_vcoord(cov_cs, start);
  assert_vcoord(cov_cs, end);
  
  return cov_cs[end] - cov_cs[start] - cov[end];
}

int VariationSet::get_segments_coverage(vector< Segment >& segs)
{
  int result = 0;
  for (unsigned int i=0; i<segs.size(); ++i)
    result += get_segment_coverage(segs[i].contig, segs[i].start, segs[i].end);
  return result;
}

void VariationSet::get_mean_sd_segments_coverage(vector< Segment >& segs, double& cov_mean, double& cov_sd)
{
  vector<double> covs;
  for (unsigned int i=0; i<segs.size(); ++i) {
    string contig = segs[i].contig;
    int start = segs[i].start;
    int end = segs[i].end;
    vector<int>& cc = m_covs[contig];
    for (int coord=start; coord<=end; ++coord)
      covs.push_back(cc[coord]);
  }
  double cov_var;
  mean_and_var(covs, cov_mean, cov_var);
  cov_sd = sqrt(cov_var);
}

void VariationSet::get_vars(const string contig, const int coord, vector<Variation>& vars)
{
  if (m_vars.find(contig) == m_vars.end())
    return;
  map< int, map <Variation, int> >& contig_vars = m_vars[contig];
  if (contig_vars.find(coord) == contig_vars.end())
    return;
  map <Variation, int>& coord_vars = contig_vars[coord];
  for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt)
    vars.push_back((*xt).first);
}

int VariationSet::get_ref_count(const string contig, const int coord) const
{
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  const vector<int>& cov = m_covs.at(contig);
  assert_vcoord(cov, coord);
  int total_count = cov[coord];

  // contig not found
  if (m_vars.find(contig) == m_vars.end())
    return total_count;
  const map< int, map <Variation, int> >& contig_vars = m_vars.at(contig);

  // coord not found
  if (contig_vars.find(coord) == contig_vars.end())
    return total_count;
  const map <Variation, int>& coord_vars = contig_vars.at(coord);

  int total_var = 0;
  for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt) {
    Variation var = (*xt).first;
    int count = (*xt).second;
    total_var += count;
  }

  return (total_count - total_var);
}

void VariationSet::get_counts(const string contig, const int coord, vector<int>& counts)
{
  counts.clear();
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  int ref_count = get_ref_count(contig, coord);
  counts.push_back(ref_count);
  if (m_vars.find(contig) == m_vars.end())
    return;
  map< int, map <Variation, int> >& contig_vars = m_vars[contig];
  map <Variation, int>& coord_vars = contig_vars[coord];
  for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt)
    counts.push_back((*xt).second);
  sort(counts.begin(), counts.end(), greater<int>());
}

bool comp_var_count_ascending (pair < Variation, int > i, pair < Variation, int > j)
{
  return (i.second < j.second);
}

bool comp_var_count_descending (pair < Variation, int > i, pair < Variation, int > j)
{
  return (i.second > j.second);
}

void VariationSet::get_var_counts(const string contig, const int coord, bool ascending,
				  vector < pair < Variation, int > >& counts) const
{
  counts.clear();
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  int ref_count = get_ref_count(contig, coord);
  if (ref_count > 0)
    counts.push_back(make_pair(Variation(), ref_count));
  if (m_vars.find(contig) == m_vars.end())
    return;
  const map< int, map <Variation, int> >& contig_vars = m_vars.at(contig);
  if (contig_vars.find(coord) == contig_vars.end())
    return;
  const map <Variation, int>& coord_vars = contig_vars.at(coord);
  for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt)
    counts.push_back(*xt);
  if (ascending)
    sort(counts.begin(), counts.end(), comp_var_count_ascending);
  else
    sort(counts.begin(), counts.end(), comp_var_count_descending);
}

Variation VariationSet::get_major(const string contig, const int coord)
{
  Variation result;
  massert(m_covs.find(contig) != m_covs.end(), "contig not found");
  vector<int>& cov = m_covs[contig];
  assert_vcoord(cov, coord);

  // contig not found
  if (m_vars.find(contig) == m_vars.end())
    return result;
  map< int, map <Variation, int> >& contig_vars = m_vars[contig];

  // coord not found
  if (contig_vars.find(coord) == contig_vars.end())
    return result;

  // init for reference count
  int ref_count = get_ref_count(contig, coord);
  int max_count = ref_count;

  // find maximal variant
  map <Variation, int>& coord_vars = contig_vars[coord];
  for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt) {
    Variation var = (*xt).first;
    int count = (*xt).second;
    if (count > max_count) {
      max_count = count;
      result = var;
    }
  }

  return result;
}

void VariationSet::restrict(vector<string>& contigs)
{
  map< string, map< int, map <Variation, int> > > vars;
  map<string, vector<int> > covs, covs_cs, covs_noise;
  for (unsigned int i=0; i<contigs.size(); ++i) {
    string contig = contigs[i];
    vars[contig] = m_vars[contig];
    covs[contig] = m_covs[contig];
    covs_cs[contig] = m_covs_cs[contig];
    covs_noise[contig] = m_covs_noise[contig];
  }
  m_vars = vars;
  m_covs = covs;
  m_covs_cs = covs_cs;
  m_covs_noise = covs_noise;
}

void VariationSet::restrict(map< string, map< int, set <Variation> > >& sites)
{
  map< string, map< int, map <Variation, int> > > vars;
  map<string, vector<int> > covs_noise(m_covs_noise);

  for (map< string, map< int, map <Variation, int> > >::const_iterator it=m_vars.begin(); it!=m_vars.end(); ++it) {
    string contig = (*it).first;
    const map< int, map <Variation, int> >& vars_contig = (*it).second;

    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      const map <Variation, int>& vars_coord = (*jt).second;
      for (map <Variation, int>::const_iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;
	if (sites.find(contig) != sites.end() &&
	    sites[contig].find(coord) != sites[contig].end() &&
	    sites[contig][coord].find(var) != sites[contig][coord].end()) {
	  vars[contig][coord][var] = count;
	} else {
	  covs_noise[contig][coord] += count;
	}
      }
    }
  }

  m_vars = vars;
  m_covs_noise = covs_noise;
}

void VariationSet::get_total_var_counts(int& count, int& unique_count)
{
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=m_vars.begin(); it!=m_vars.end(); ++it) {
    string contig = (*it).first;
    const map< int, map <Variation, int> >& vars_contig = (*it).second;
    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      const map <Variation, int>& vars_coord = (*jt).second;
      for (map <Variation, int>::const_iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int vcount = (*xt).second;
	count += vcount;
	unique_count += 1;
      }
    }
  }
}

void VariationSet::get_contig_data(string contig, int from, int to,
				   vector<string>& contig_r,
				   vector<int>& coord_r,
				   vector<string>& var_r,
				   vector<int>& count_r,
				   vector<int>& total_r,
				   vector<int>& cc_start_r,
				   vector<int>& cc_end_r)
{
  if (m_covs.find(contig) == m_covs.end())
    return;
  
  vector<int>& contig_cov = m_covs[contig];
  map< int, map <Variation, int> >& contig_vars = m_vars[contig];

  for (int coord=from; coord<=to; ++coord) {
    if (coord < 0 || coord >= (int)contig_cov.size()) continue;

    int total_cov = contig_cov[coord];
    int cumsum = 0;

    // start with variations
    if (contig_vars.find(coord) != contig_vars.end()) {
      map <Variation, int>& coord_vars = contig_vars[coord];
      for (map <Variation, int>::const_iterator xt=coord_vars.begin(); xt!=coord_vars.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;
	contig_r.push_back(contig);
	coord_r.push_back(coord);
	var_r.push_back(var.to_string());
	count_r.push_back(count);
	total_r.push_back(total_cov);
	cc_start_r.push_back(cumsum);
	cc_end_r.push_back(cumsum+count);
	cumsum += count;
      }
    }

    // add ref
    int ref_cov = total_cov - cumsum;
    if (ref_cov != 0 || total_cov == 0) {
      contig_r.push_back(contig);
      coord_r.push_back(coord);
      var_r.push_back("none");
      count_r.push_back(ref_cov);
      total_r.push_back(total_cov);
      cc_start_r.push_back(cumsum);
      cc_end_r.push_back(cumsum+ref_cov);
    }
  }
}

void VariationSet::get_summary(string contig, VariType vtype,
			       int& from, int& to, int min_cov, double max_seg_freq,
			       int& seg_count_r, int& median_cov_r)
{
  seg_count_r = 0;
  median_cov_r = 0;

  // skip if no contig found
  if (m_covs.find(contig) == m_covs.end()) {
    cout << "warning: contig not found: " << contig << endl;
    cout << "number of contigs: " << m_covs.size() << endl;
    cout << "first contig: " << (*m_covs.begin()).first << endl;
    from = -1;
    to = -1;
    return;
  }

  vector<int>& contig_cov = m_covs[contig];
  map< int, map <Variation, int> >& vars_contig = m_vars[contig];

  // fit from left side
  if (from < 0)
    from = 0;
  if (to < 0)
    to = 0;

  // fit from right side
  if (to >= (int)contig_cov.size())
    to = contig_cov.size()-1;
  if (from >= (int)contig_cov.size())
    from = contig_cov.size()-1;

  // skip if segment is invalid 
  if (from > to) {
    from = -2;
    to = -2;
    return;
  }
  
  vector<int> acounts;

  // count sites
  for (int coord=from; coord<=to; ++coord) {
    map< int, map <Variation, int> >::const_iterator it = vars_contig.find(coord);
    if (it == vars_contig.end())
      continue;

    Variation var = get_major(contig, coord);
    get_counts(contig, coord, acounts);
    if (acounts.size() < 2 || (vtype != vtNone && !var.has_type(vtype)))
      continue;
    int total = get_coverage(contig, coord);
    double freq2 = (double)acounts[1] / (double)total;
    if (total >= min_cov && (freq2 >= (1-max_seg_freq) || !var.is_ref()))
      seg_count_r++;
  }

  // median coverage
  vector<int> counts(contig_cov.begin() + from, contig_cov.begin() + to + 1);
  median_cov_r = median<int>(counts);
}

void VariationSet::get_freq_summary(string contig, int& from, int& to,
				    int& count_r, double& median_freq_r)
{
  count_r = 0;
  median_freq_r = 0;
  
  // skip if no contig found
  if (m_covs.find(contig) == m_covs.end()) {
    from = -1;
    to = -1;
    return;
  }

  vector<int>& contig_cov = m_covs[contig];
  map< int, map <Variation, int> >& vars_contig = m_vars[contig];

  // fit from left side
  if (from < 0)
    from = 0;
  if (to < 0)
    to = 0;

  // fit from right side
  if (to >= (int)contig_cov.size())
    to = contig_cov.size()-1;
  if (from >= (int)contig_cov.size())
    from = contig_cov.size()-1;

  // skip if segment is invalid 
  if (from > to) {
    from = -2;
    to = -2;
    return;
  }
  
  vector<int> acounts;
  vector<double> freqs;

  // count sites
  for (int coord=from; coord<=to; ++coord) {
    map< int, map <Variation, int> >::const_iterator it = vars_contig.find(coord);
    if (it == vars_contig.end())
      continue;

    Variation var = get_major(contig, coord);
    get_counts(contig, coord, acounts);
    if (acounts.size() < 2)
      continue;
    int total = get_coverage(contig, coord);
    double freq = 1.0 - (double)acounts[0]/(double)total;
    count_r++;
    freqs.push_back(freq);
  }

  if (freqs.size() > 0)
    median_freq_r = median<double>(freqs);
}

void VariationSet::get_cov_summary(string contig, int& from, int& to,
				   int min_cov, int& median_cov_r)
{
  median_cov_r = 0;

  // skip if no contig found
  if (m_covs.find(contig) == m_covs.end()) {
    from = -1;
    to = -1;
    return;
  }

  vector<int>& contig_cov = m_covs[contig];

  // fit from left side
  if (from < 0)
    from = 0;
  if (to < 0)
    to = 0;

  // fit from right side
  if (to >= (int)contig_cov.size())
    to = contig_cov.size()-1;
  if (from >= (int)contig_cov.size())
    from = contig_cov.size()-1;

  // skip if segment is invalid 
  if (from > to) {
    from = -2;
    to = -2;
    return;
  }

  // median coverage
  vector<int> counts(contig_cov.begin() + from, contig_cov.begin() + to + 1);
  median_cov_r = median<int>(counts);
}

// get mean and sd of coverage for segment
void VariationSet::get_segment_stats(string contig, int start, int end,
				     double& mean_r, double& var_r)
{
  massert(m_covs.find(contig) != m_covs.end(), "contig %s not found", contig.c_str());
  vector<int>& contig_cov = m_covs[contig];
  int max_coord = (int)contig_cov.size()-1;
  massert(start >= 0 && end >= 0 && start <= end && start <= max_coord && end <= max_coord, "coords out of range");

  vector<int> v(contig_cov.begin() + start, contig_cov.begin() + end + 1);

  double sum = accumulate(v.begin(), v.end(), 0.0);
  mean_r = sum / v.size();
  
  vector<double> diff(v.size());
  transform(v.begin(), v.end(), diff.begin(),
	    bind2nd(minus<double>(), mean_r));
  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  var_r = sq_sum / v.size();
}
