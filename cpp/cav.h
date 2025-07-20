int construct_main(const char* name, int argc, char **argv);
int merge_main(const char* name, int argc, char **argv);
int combine_main(const char* name, int argc, char **argv);
int filter_main(const char* name, int argc, char **argv);
int restrict_main(const char* name, int argc, char **argv);
int diff_main(const char* name, int argc, char **argv);

// refine output of binner such as metaBAT2
int refine_bins_main(const char* name, int argc, char **argv);

// extract site trajectory
int site_trajectory_main(const char* name, int argc, char **argv);

// extract mean vcluster trajectory
int vcluster_trajectory_main(const char* name, int argc, char **argv);

// extract mean bin trajectory
int bin_trajectory_main(const char* name, int argc, char **argv);

// cluster variants within bins
int variant_cluster_main(const char* name, int argc, char **argv);

// lib stats
int cav_stats_main(const char* name, int argc, char **argv);

// create bin fasta files
int fasta_main(const char* name, int argc, char **argv);

int cov_mat_main(const char* name, int argc, char **argv);

int info_main(const char* name, int argc, char **argv);

int refine_global_main(const char* name, int argc, char **argv);
int refine_local_main(const char* name, int argc, char **argv);
int dump_local_scores_main(const char* name, int argc, char **argv);
int sites_main(const char* name, int argc, char **argv);
int bin_main(const char* name, int argc, char **argv);

int read_query_main(const char* name, int argc, char **argv);

int view_main(const char* name, int argc, char **argv);
int dump_main(const char* name, int argc, char **argv);


// deprecated functions
// see metagenomics/dep/nlv

/* int compare_main(const char* name, int argc, char **argv); */
/* int divergence_main(const char* name, int argc, char **argv); */
/* int query_main(const char* name, int argc, char **argv); */
/* int query_nts_main(const char* name, int argc, char **argv); */
/* int coverage_main(const char* name, int argc, char **argv); */
/* int segregation_main(const char* name, int argc, char **argv); */
/* int restrict_main(const char* name, int argc, char **argv); */
/* int sites_main(const char* name, int argc, char **argv); */


