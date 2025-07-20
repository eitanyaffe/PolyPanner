#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cav.h"

using namespace std;

void usage(const char* name)
{
  fprintf(stderr, "polypanner: Co-abundance variant inference algorithm\n");
  fprintf(stderr, "usage: %s <command> [options]\n", name);
  fprintf(stderr, "commands:\n");
  fprintf(stderr, "  construct: Construct POP from mapped reads\n");
  fprintf(stderr, "  info: Show basic info and stats for POP file\n");
  fprintf(stderr, "  filter: Infer sequencing error rates and true segregating sites\n");
  fprintf(stderr, "  restrict: Restrict POP to selected variants\n");
  fprintf(stderr, "  refine_global: Breakdown contigs into segments, using a global coverage distribution per contig\n");
  fprintf(stderr, "  refine: Breakdown contigs into segments, on dangle coords with a transition in coverage distribution\n");
  fprintf(stderr, "  fasta: Create fasta files for bins\n");
  fprintf(stderr, "  site_trajectory: Extract site trajectories\n");
  fprintf(stderr, "  bin_trajectory: Mean bin trajectories\n");
  fprintf(stderr, "  vcluster_trajectory: Mean vcluster trajectories\n");
  fprintf(stderr, "  stats: Save table with read depth and max read length\n");
  fprintf(stderr, "  variant_cluster: Cluster variants within genome bin\n");
  fprintf(stderr, "  cov_matrix: Generate coverage matrix for segments\n");
  fprintf(stderr, "  refine_bins: Refine bins\n");
  fprintf(stderr, "  sites: Annotate variants\n");
  fprintf(stderr, "  merge: Merge multiple POPs (adding)\n");
  fprintf(stderr, "  combine: Combine POPs from different subjects\n");
  fprintf(stderr, "  restrict: Restrict POP to a set of contigs\n");
  fprintf(stderr, "  diff: Compare two POP files and report differences\n");
  fprintf(stderr, "  dump: Dump single POP into tab-delimited tables (for debugging)\n");
  fprintf(stderr, "  dump_scores: Dump local refinement scores along contigs (for debugging)\n");
  fprintf(stderr, "  read_query: Extract reads covering specific segment (for debugging)\n");
  //  fprintf(stderr, "  compare: Out table comparing multiple POPs\n");
  //  fprintf(stderr, "  diverge: Identify diverging major alleles between two POPs\n");
  //  fprintf(stderr, "  segregation: Identify segtrgating positions in a single POP\n");
  //  fprintf(stderr, "  query: Extract counts for a set of contig/coord/variation\n");
  //  fprintf(stderr, "  query_nts: Extract sub_A/sub_C/sub_G/sub_T/ref counts for a set of contig/coord\n");
  //  fprintf(stderr, "  coverage: Extract total median coverage for segment sets\n");
  //  fprintf(stderr, "  view: Print to screen data for single contig\n");
}

int main(int argc, char **argv)
{
  if (argc == 1) {
    usage(argv[0]);
    exit(1);
  }
  string command(argv[1]);
  string name =  string(argv[0]) + " " + command;

  int rc = 0;
  if (command == "construct") {
    rc = construct_main(name.c_str(), argc-1, argv+1);
  } else if (command == "info") {
    rc = info_main(name.c_str(), argc-1, argv+1);
  } else if (command == "filter") {
    rc = filter_main(name.c_str(), argc-1, argv+1);
  } else if (command == "restrict") {
    rc = restrict_main(name.c_str(), argc-1, argv+1);
  } else if (command == "dump") {
    rc = dump_main(name.c_str(), argc-1, argv+1);
  } else if (command == "merge") {
    rc = merge_main(name.c_str(), argc-1, argv+1);
  } else if (command == "combine") {
    rc = combine_main(name.c_str(), argc-1, argv+1);
  } else if (command == "cov_matrix") {
    rc = cov_mat_main(name.c_str(), argc-1, argv+1);
  } else if (command == "refine") {
    rc = refine_local_main(name.c_str(), argc-1, argv+1);
  } else if (command == "fasta") {
    rc = fasta_main(name.c_str(), argc-1, argv+1);
  } else if (command == "stats") {
    rc = cav_stats_main(name.c_str(), argc-1, argv+1);
  } else if (command == "vcluster_trajectory") {
    rc = vcluster_trajectory_main(name.c_str(), argc-1, argv+1);
  } else if (command == "bin_trajectory") {
    rc = bin_trajectory_main(name.c_str(), argc-1, argv+1);
  } else if (command == "site_trajectory") {
    rc = site_trajectory_main(name.c_str(), argc-1, argv+1);
  } else if (command == "variant_cluster") {
    rc = variant_cluster_main(name.c_str(), argc-1, argv+1);
  } else if (command == "refine_bins") {
    rc = refine_bins_main(name.c_str(), argc-1, argv+1);
  } else if (command == "refine_global") {
    rc = refine_global_main(name.c_str(), argc-1, argv+1);
  } else if (command == "sites") {
    rc = sites_main(name.c_str(), argc-1, argv+1);
  } else if (command == "dump_scores") {
    rc = dump_local_scores_main(name.c_str(), argc-1, argv+1);
  } else if (command == "read_query") {
    rc = read_query_main(name.c_str(), argc-1, argv+1);
  } else if (command == "diff") {
    rc = diff_main(name.c_str(), argc-1, argv+1);
  } else {
    printf("unknown command: %s\n", command.c_str());
    usage(argv[0]);
    exit(1);
  }

  return rc;
}
