# PolyPanner: Detecting polymorphic sites in complex microbial communities

PolyPanner is a tool designed for the detection of polymorphic variants in metagenomic assembled genomes of complex bacterial populations. It leverages repeat sampling (either spatial or temporal) to improve metagenomic assembly quality and identify high confidence polymorphic variants, called dynamic variants. Dynamic variants are polymorphic variants with a non-constant frequency over samples, assessed while taking into account stochasticity in read depth. The remaining variants, called spurious variants, reflect homology between distinct populations (ortholog spurious variants) or homology within a population due to duplicated genes (paralog spurious variants) and are considered noise for our purpose of tracking competing alleles within bacterial populations. All formal definitions and a complete description of the implementation are detailed [here]( https://www.biorxiv.org/content/10.1101/2023.09.04.556257v3).

As input, PolyPanner receives a metagenomic co-assembly (FASTA format) and a set of aligned shotgun DNA libraries (mapped with bwa to the assembly). PolyPanner transforms the read alignments to single-nucleotide coverage vectors that represent sample-specific read counts of perfect and mismatch alignments at each base pair in the assembly. Tasks performed by PolyPanner are:
1. The removal of sequencing errors.
2. The refinement of assembly contigs into segments, through the introduction of breakpoints where there are abrupt transitions in coverage. These transitions likely reflect chimeric assembly breakpoints where the two sides of the breakpoint represent different populations.
3. The refinement of genome bins, through the removal of segments that differ in their coverage profiles.
4. The identification of dynamic variants.

The first two tasks are performed prior to genome binning, which is done by an external tool.  After binning PolyPanner proceeds to refine genome bins and to identify dynamic variants.

PolyPanner was developed by Eitan Yaffe (eitan.yaffe@gmail.com).

## Installation

Clone the PolyPanner repository to a local directory of your choice:
```
git clone https://github.com/eitanyaffe/PolyPanner.git
```

Use pre-compiled PolyPanner binaries or compile PolyPanner from code, as follows.

### Use pre-compiled binaries

Binaries for macOS (Ventura 13.3.1) or Ubuntu (20.04.1) are part of the v1.0 release. You may need to run ```chmod +x $BINARY``` if they are not executable after downloading them.

MacOS: https://github.com/eitanyaffe/PolyPanner/releases/download/v1.0.0/polypanner-v1.0.0-macOS-i386.macos

Linux: https://github.com/eitanyaffe/PolyPanner/releases/download/v1.0.0/polypanner-v1.0.0-linux-x86-64.linux

PolyPanner binaries depend on gsl.

On MacOS run:
```
brew install gsl
```

On Ubuntu run:
```
apt-get update
apt-get install libgsl0-dev
```

If you are getting an error message such as ```error while loading shared libraries: libgsl.so.23``` you are likely using a new OS version (e.g. Ubuntu 22.04 instead of 20.04), and you need to compile PolyPanner from code.

MacOS may not let you run PolyPanner because it can't verify the developer. If that is the case, go to System Settings > Privacy & Security > Security and confirm a security bypass.

### Compile from code

Compilation was tested on macOS (Ventura 13.3.1) or Ubuntu (20.04.1). You will need to install gcc, gsl and boost. Compilation should take under an hour.

On MacOS run:
```
brew install g++ gsl boost
```

On Ubuntu run:
```
apt-get update
apt-get install git build-essential libgsl0-dev libboost-all-dev
```

Clone the repository from github, and compile PolyPanner:
```
cd PolyPanner
make
```
### Dependencies

There are several helper scripts written in Perl and R. These scripts were tested on Perl versions 5.30.3 and 5.30.0, and on R versions 4.3.1 and 4.4.1, but they are likely compatible with other comparable versions.

### Validate installation

To validate the installation you can run the following unit test which should take under 20 seconds to complete. Docker is required to run metaBAT2 during the unit test. The unit test runs on a small mock dataset composed of 4 DNA libraries, each with ~20k paired reads (2x150nt). 

```
make test
```

A recommended way to start using PolyPanner is by customizing the commands of the unit test. The bash commands can be printed by peforming a dry-run:

```
make test -n
```

## Workflow overview


The standard PolyPanner workflow is detailed below. The syntax of PolyPanner commands and utility scripts, including user-defined parameters, is documented [here](docs/syntax.md).

### 1. Input

The workflow input is composed of:
- A co-assembly with a contig fasta file and a contig table (fields: contig and length).
- A set of mapped paired DNA libraries in the SAM format (mapped separately for R1 and R2).

See files under the directory ```input/test``` for an example of input files.

### 2. POP construction

The following three steps are done for each library: (1) SAM files are first converted to a tabular format using the ```utils/parse_bwa_sam.pl``` script, (2) read sides are paired using the ```utils/pair_reads.pl``` script, and (3) paired reads are transformed into POP files, which represent mapped reads in an efficient manner that allows quick queries on each library. POP files are internally compressed using gzip.

### 3. Removal of sequencing errors

Libraries are merged into a single library using the ``polypanner merge`` command. Sequencing errors are identified and removed using the ``polypanner filter`` command. Libraries are restricted to true segregating sites using the ``polypanner restrict`` command. A tab-delimited table with paths to the restricted libraries is generated using the ``utils/make_pop_table.r`` script.

### 4. Assembly refinement

Contigs are refined into segments using the ``polypanner refine`` command.

### 5. Genome binning 

Genome binning is performed by an external binning tool. To use metaBAT2, a segment coverage matrix (input for metaBAT2) is generated using the ``polypanner cov_matrix``. The output of the binning step is required to be a tab-delimited segment-bin table  with 2 columns: segment and bin. The workflow can support other binning methods, just make sure to generate the segment-bin table as specified. 

### 6. Genome refinement 

The segment-bin table is processed using the ``utils/bin_summary.r`` script. Bins are refined using the ``polypanner refine_bins`` command.

### 7. Site annotation

Sites are annotated using the ``polypanner sites`` command.

### 8. Output

Genome bin trajectories are generated using the ``polypanner bin_trajectory`` command. Site trajectories are generated using the ``polypanner site_trajectory``.

Other output files include the refined segment-bin table and the annotated site table. 

## Output formats

### Site table

Fields:
- vid: Variant identifer.
- contig/coord: The contig and coordinate within contig.
- variant: String identifier of the variant.
- segment: The segment.
- bin: The genome bin.
- var_count: Number of reads supporting the variant.
- total_count: Number of reads supporting the position.
- n_samples: Number of samples that have reads supporting the variant.
- is_internal: Is variant close to segment edge.
- variant_p: p-value of the test comparing the variant and the local coverages.
- regional_p: p-value of the test comparing the local and regional coverages.
- comp_p: p-value of the test comparing the complementary and regional coverages.

Dynamic sites are defined as sites for which is_internal equals T, variant_p<X, regional_p>X and comp_p<X, where X=0.05.

See ```output/test/sites_annotated``` in the test output for an example.

### Segment table

Fields:
- segment: segment identifier.
- contig/start/end: The contig and coordinates of segment within contig.
- bin_org: Original bin assigned by binner.
- bin: Final bin identifier, after bin trimming.

See ```output/test/seg_bin_final``` in the test output for an example.
