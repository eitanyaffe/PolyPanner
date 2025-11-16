#ifndef __UTIL__
#define __UTIL__

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <numeric>

#include "Variation.h"

using namespace std;

template <typename T> T median(vector<T> &v);

void split_string(string in, vector<string> &fields, char delim);
void split_line(istream &in, vector<string> &fields, char delim);
void massert(bool cond, char *fmt, ...);
void mexit(char *fmt, ...);
int get_field_index(string field, const vector<string>& titles);

// error message: "<title> <str> could not be converted an integer"
int safe_string_to_int(const string str, const string title);

char index2char(int i);
int char2index(char c);

string reverse_complement(string seq);

void read_sites(string fn, map< string, map< int, set< Variation > > >& keys);

void read_library_table(string fn, vector< string >& ifns, vector< string >& library_ids);

// read contig lengths
void read_contig_length(string fn, map<string, int>& contig_map);

// passes through origin
double slope_origin(const vector<double>& x, const vector<double>& y);

void parse_cigar(string in, vector<pair < char, int> > & result);

int get_of_files(string idir);

// sorts vector and returns threshold p-value given FDR 
double apply_Benjamini_Hochberg(vector<double> pvals, double fdr, int test_count);

// load/save fasta to a map contig->nts
void load_fasta(string fn, map<string, string>& fasta);
void save_fasta(string fn, map<string, string>& fasta);

// generate random nts
string rand_nts(int len);

// pearson test
template <typename T> double pearson_test(vector<T>& v1, vector<T>& v2, double* p_value);

// chi-square test
template <typename T> double chi_square_test(vector<T>& v1, vector<T>& v2, double* p_value);
template <typename T> double chi_square_test_p(vector<T>& v1, vector<T>& v2);

unsigned long long getTotalSystemMemory();

void read_binned_segments(string fn, string bin_field, map< string, vector< Segment > >& bins);

// map string->string
void read_string_map(string fn, string key_field, string value_field, map< string, string >& result);

void read_segment_table(string fn, vector< Segment > & segs);

void mean_and_var(const vector<double>& vec, double& mean, double& var);

#endif
