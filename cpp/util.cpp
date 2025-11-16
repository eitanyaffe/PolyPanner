#include "util.h"
#include <numeric>

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
using namespace boost::math;

template <typename T> T median(vector<T>& v)
{
  size_t n = v.size() * 0.5;
  if (n >= v.size()) n = v.size()-1;
  nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}
template double median<double>(vector<double> &v);
template int median<int>(vector<int> &v);

string reverse_complement(string seq) {
  string result = seq;
  int N = seq.length();
  for (int i=0; i<N; i++) {
    char c = seq[N-i-1], r;
    switch( c )
      {
      case 'A': r = 'T'; break;
      case 'G': r = 'C'; break;
      case 'C': r = 'G'; break;
      case 'T': r = 'A'; break;
      default: r = c;
      }
    result[i] = r;
  }
  return(result);
}

void split_string(string in, vector<string> &fields, char delim)
{
  fields.resize(0);
  string field;
  for (unsigned int i=0; i<in.size(); ++i) {
    char c = in[i];
    if (c == -1)
      break;
    if(c == '\r') {
      continue;
    }
    if(c == '\n') {
      fields.push_back(field);
      field.resize(0);
      break;
    }
    if(c == delim) {
      fields.push_back(field);
      field.resize(0);
    } else {
      field.push_back(c);
    }
  }

  if (field.length() > 0)
    fields.push_back(field);
}

void split_line(istream &in, vector<string> &fields, char delim)
{
  fields.resize(0);
  string field;
  while(in) {
    char c = in.get();
    if (c == -1)
      break;
    if(c == '\r') {
      continue;
    }
    if(c == '\n') {
      fields.push_back(field);
      field.resize(0);
      break;
    }
    if(c == delim) {
      fields.push_back(field);
      field.resize(0);
    } else {
      field.push_back(c);
    }
  }

  if (field.length() > 0)
    fields.push_back(field);
}

void mexit(char *fmt, ...)
{
  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

void massert(bool cond, char *fmt, ...)
{
  if (cond)
    return;

  fprintf(stderr, "Error: ");

  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  va_end(argp);

  fprintf(stderr, "\n");
  exit(-1);
}

int safe_string_to_int(const string str, const string title)
{
    try {
      return stoi(str);
    }
    catch (std::invalid_argument const &e) {
      mexit("%s %s could not be converted to integer", title.c_str(), str.c_str());
    }
    return 0;
}

int get_field_index(string field, const vector<string>& titles)
{
  int result = -1;
  for (unsigned int i = 0; i<titles.size(); i++)
    if (titles[i] == field)
      result = i;
  if (result == -1) {
    cout << "unknown field " << field << endl;
    exit(1);
  }
  return result;
}

char index2char(int i) {
  switch (i)
    {
    case 0: return('A');
    case 1: return('C');
    case 2: return('G');
    case 3: return('T');
    }
  massert(0, "unknown character index");
  return 0;
}

int char2index(char c) {
  switch (c)
    {
    case 'A': return(0);
    case 'C': return(1);
    case 'G': return(2);
    case 'T': return(3);
    }
  massert(0, "unknown character: %c", c);
  return 0;
}

void read_sites(string fn, map< string, map< int, set< Variation > > >& keys)
{
  cout << "reading site table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_ind = get_field_index("contig", fields);
  int coord_ind = get_field_index("coord", fields);
  int var_ind = get_field_index("variant", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[contig_ind];
    int coord = safe_string_to_int(fields[coord_ind], "coordinate") - 1;
    Variation var = Variation(fields[var_ind]);
    keys[contig][coord].insert(var);
  }
}

void read_library_table(string fn, vector< string >& ifns, vector< string >& library_ids)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  ifns.clear();
  library_ids.clear();

  // parse header
  split_line(in, fields, delim);
  int fn_ind = get_field_index("fn", fields);
  int lib_id_ind = get_field_index("lib_id", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string ifn = fields[fn_ind];
    ifns.push_back(ifn);
    
    // extract library ID if available, otherwise generate from filename
    if (lib_id_ind >= 0 && lib_id_ind < static_cast<int>(fields.size())) {
      library_ids.push_back(fields[lib_id_ind]);
    } else {
      // fallback: use filename basename without extension
      string lib_id = ifn;
      size_t last_slash = lib_id.find_last_of("/\\");
      if (last_slash != string::npos) {
        lib_id = lib_id.substr(last_slash + 1);
      }
      size_t last_dot = lib_id.find_last_of(".");
      if (last_dot != string::npos) {
        lib_id = lib_id.substr(0, last_dot);
      }
      library_ids.push_back(lib_id);
    }
  }
}

double slope_origin(const vector<double>& x, const vector<double>& y)
{
  const auto s_xx = inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const auto s_xy = inner_product(x.begin(), x.end(), y.begin(), 0.0);
  const auto a    = s_xy / s_xx;
  return a;
}

void parse_cigar(string in, vector<pair < char, int> > & result)
{
  result.resize(0);
  string field;
  int length = 0;
  for (unsigned int i=0; i<in.size(); ++i) {
    char c = in[i];
    if (isdigit(c))
      length = length*10 + (c - '0');
    else {
      result.push_back(make_pair(c, length));
      length = 0;
    }
  }
}

int get_of_files(string idir)
{
  DIR *dir;
  struct dirent *ent;
  massert((dir = opendir (idir.c_str())) != NULL, "could not open input directory");
  int file_count = 0;
  while ((ent = readdir (dir)) != NULL) {
    string ifn = ent->d_name;
    if (ifn.find("fastq") == string::npos || ifn.find("~") != string::npos)
      continue;
    file_count++;
  }
  closedir (dir);
  return file_count;
}

double apply_Benjamini_Hochberg(vector<double> pvals, double alpha, int test_count)
{
  massert(pvals.size() > 0, "no pvalues supplied");
  sort(pvals.begin(), pvals.end());
  double K = pvals.size()-1;
  while(K > 0 && pvals[K] > ((alpha*K) / test_count))
    K--;
  return pvals[K];
}

void load_fasta(string fn, map<string, string>& fasta)
{
  massert(fn != "", "assembly file not defined");
  cout << "loading assembly file (fasta): " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());

  string contig = "";
  string seq = "";
  while(1) {
    string line;
    getline(in, line);
    bool eof = in.eof();
    if (eof)
      break;
    if (!eof && !in.good()) {
      cerr << "Error reading line, rdstate=" << in.rdstate() << endl;
      exit(-1);
    }
    massert(line.length() > 0, "empty line in fasta file");
    if (line[0] == '>') {
      if (contig != "")
	fasta[contig] = seq;
      contig = line.substr(1);
      if (contig.find(' ') != string::npos)
	contig = contig.substr(0, contig.find(' '));
      seq = "";
    } else {
      seq += line;
    }
  }
  if (contig != "")
    fasta[contig] = seq;
  in.close();
}

void save_fasta(string fn, map<string, string>& fasta)
{
  cout << "saving fasta file: " << fn << endl;
  ofstream out(fn.c_str());
  massert(out.is_open(), "could not open file %s", fn.c_str());
  for (map<string, string>::iterator it=fasta.begin(); it!=fasta.end(); ++it) {
    string name = (*it).first;
    string seq = (*it).second;
    out << ">" << name << endl;
    out << seq << endl;
  }
}

string rand_nts(int len)
{
  string result;
  char nts[4] = {'A', 'C', 'G', 'T'};
  for (int i=0; i<len; ++i)
    result += nts[(int)rand()%4];
  return result;
}

template <typename T> double pearson_test(vector<T>& x, vector<T>& y, double* p_value)
{
  const auto N = x.size();
  const auto s_x = accumulate(x.begin(), x.end(), 0.0) / N;
  const auto s_y = accumulate(y.begin(), y.end(), 0.0) / N;
			      
  const auto s_xx = inner_product(x.begin(), x.end(), x.begin(), 0.0) / N;
  const auto s_yy = inner_product(y.begin(), y.end(), y.begin(), 0.0) / N;

  const auto s_xy = inner_product(x.begin(), x.end(), y.begin(), 0.0) / N;

  auto rho = (s_xy - s_x*s_y) / (sqrt(s_xx-s_x*s_x) * sqrt(s_yy-s_y*s_y));

  if (p_value) {
    if (abs(rho) < 1) {
      double stat = rho * sqrt((N-2)/(1-rho*rho));
      students_t t_distrib(N-2);
      *p_value = 2 * cdf(t_distrib, -abs(stat));
    } else {
      *p_value = 0;
    }
  }
  
  return rho;
}
template double pearson_test(vector<int>& x, vector<int>& y, double* p_value);
template double pearson_test(vector<double>& x, vector<double>& y, double* p_value);

unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

template <typename T> double chi_square_test(vector < T >& c1, vector < T >& c2, double* p_value)
{
  massert(c1.size() == c2.size(), "vector length not equal");
  int n = c1.size();

  // get totals
  double tot1 = 0;
  double tot2 = 0;
  for (int i=0; i<n; ++i) {
    tot1 += c1[i];
    tot2 += c2[i];
  }
  double tot = tot1 + tot2;
  double f1 = tot1/tot;
  double f2 = tot2/tot;

  // get lib freq
  vector <double> freq_lib(n);
  for (int i=0; i<n; ++i)
    freq_lib[i] = (c1[i] + c2[i]) / tot;

  double result = 0;
  for (int i=0; i<n; ++i) {
    double exp1 = tot * f1 * freq_lib[i];
    double exp2 = tot * f2 * freq_lib[i];
    result += ((c1[i] - exp1) * (c1[i] - exp1)) / exp1;
    result += ((c2[i] - exp2) * (c2[i] - exp2)) / exp2;
  }

  if (p_value) {
    chi_squared chi_squared(n-1);
    *p_value = (1 - cdf(chi_squared, result));
  }

  return result;
}
template double chi_square_test(vector <double>& c1, vector <double>& c2, double* p_value);

template <typename T> double chi_square_test_p(vector < T >& c1, vector < T >& c2)
{
  double result = 0;
  chi_square_test(c1, c2, &result);
  return result;
}
template double chi_square_test_p(vector <double>& c1, vector <double>& c2);


void read_binned_segments(string fn, string bin_field, map< string, vector< Segment > >& bins)
{
 cout << "reading binned segment table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int seg_ind = get_field_index("segment", fields);
  int contig_ind = get_field_index("contig", fields);
  int start_ind = get_field_index("start", fields);
  int end_ind = get_field_index("end", fields);
  int bin_ind = get_field_index(bin_field, fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string seg_id = fields[seg_ind];
    string contig = fields[contig_ind];
    int start = safe_string_to_int(fields[start_ind], "coordinate") - 1;
    int end = safe_string_to_int(fields[end_ind], "coordinate") - 1;
    string bin = fields[bin_ind];

    Segment seg(seg_id, contig, start, end);
    bins[bin].push_back(seg);
  }
}

void read_string_map(string fn, string key_field, string value_field, map< string, string >& result)
{
 cout << "reading binned segment table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int key_ind = get_field_index(key_field, fields);
  int value_ind = get_field_index(value_field, fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string key = fields[key_ind];
    string value = fields[value_ind];

    result[key] = value;
  }
}

void mean_and_var(const vector<double>& vec, double& mean, double& var) {
    const size_t N = vec.size();
    if (N == 1) {
      mean = vec[0];
      var = 0;
      return;
    }

    // Calculate the mean
    mean = std::accumulate(vec.begin(), vec.end(), 0.0) / N;

    // Now calculate the variance
    auto variance_func = [&mean, &N](double accumulator, const double& val) {
        return accumulator + ((val - mean)*(val - mean) / (N - 1));
    };
    var = std::accumulate(vec.begin(), vec.end(), 0.0, variance_func);
}

void read_contig_length(string fn, map<string, int>& contig_map)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_ind = get_field_index("contig", fields);
  int length_ind = get_field_index("length", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[contig_ind];
    int length = (int)atof(fields[length_ind].c_str());
    contig_map[contig] = length;
  }
}
