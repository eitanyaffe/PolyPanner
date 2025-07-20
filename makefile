all: bin/polypanner

###############################################################################################
# flags
###############################################################################################

# Ubuntu/Debian users: install required packages with:
# sudo apt-get update  
# sudo apt-get install build-essential libboost-all-dev libgsl-dev

UNAME_S:=$(shell uname -s)

ifeq ($(UNAME_S),Linux)
CFLAGS=-Wall -Wno-write-strings -std=c++14 -fext-numeric-literals
LDFLAGS=-pthread -lgsl -lgslcblas -lboost_iostreams
endif

ifeq ($(UNAME_S),Darwin)
UNAME_M:=$(shell uname -m)

# Apple Silicon (ARM64) - uses Homebrew in /opt/homebrew
ifeq ($(UNAME_M),arm64)
CFLAGS=-Wall -Wno-write-strings -std=c++14 -Wno-pragmas \
-I/opt/homebrew/opt/boost/include -I/opt/homebrew/opt/gsl/include
LDFLAGS=-pthread -lgsl -lgslcblas -lboost_iostreams \
-L/opt/homebrew/opt/boost/lib -L/opt/homebrew/opt/gsl/lib
endif

# Intel x86_64 - uses Homebrew in /usr/local  
ifeq ($(UNAME_M),x86_64)
CFLAGS=-Wall -Wno-write-strings -std=c++14 -Wno-pragmas \
-I/usr/local/opt/boost/include -I/usr/local/opt/gsl/include
LDFLAGS=-pthread -lgsl -lgslcblas -lboost_iostreams \
-L/usr/local/opt/boost/lib -L/usr/local/opt/gsl/lib
endif

endif

CC=g++

###############################################################################################
# rules
###############################################################################################

INSTALL_DIR=/makeshift-mnt/bin

HFILES=$(addprefix cpp/,\
Variation.h VariationSet.h cav.h Params.h util.h BinMatrix.h Filter.h Resource.h \
RefineLocal.h)

OBJ=$(addprefix obj/,\
cav_diff.o \
cav_sites.o \
cav_fasta.o \
cav_stats.o \
cav_bin_trajectory.o \
cav_vcluster_trajectory.o \
cav_site_trajectory.o \
Resource.o \
ClusterVariants.o \
cav_variant_cluster.o \
cav_refine_bins.o \
BinMatrix.o \
Filter.o \
Variation.o VariationSet.o Params.o util.o \
Dissolve.o RefineLocal.o \
cav.o cav_construct.o cav_filter.o \
cav_merge.o \
cav_combine.o \
cav_refine_local.o cav_refine_global.o \
cav_dump_local_scores.o \
cav_read_query.o \
cav_dump.o \
cav_restrict.o \
cav_info.o \
cav_cov_mat.o)

obj:
	@mkdir -p $@
bin:
	@mkdir -p $@

obj/%.o: cpp/%.cpp $(HFILES)
	$(CC) $(CFLAGS) -c -o $@ $<

bin/polypanner: obj bin $(OBJ)
	$(CC) -o $@ $(OBJ) $(LDFLAGS)

install: bin/polypanner
	cp bin/polypanner $(INSTALL_DIR)

clobber:
	rm -rf bin/polypanner obj

###############################################################################################
# unit test
###############################################################################################

# shared test variables
include mk/base.mak

# get input files
include mk/input.mak

# construct PP files
include mk/construct.mak

# remove sequencing errors
include mk/seq_errors.mak

# refine assembly
include mk/refine.mak

# infer genomes
include mk/genomes.mak

# infer sites
include mk/sites.mak

# ouput trajectotries
include mk/trajectories.mak

test:
	@echo "###################################"
	@echo "######### UNIT TEST START #########"
	@echo "###################################\n"
	@echo "\n######### 1. constructing PP files"
	$(MAKE) construct_all
	@echo "\n######### 2. removing seqeuncing errors"
	$(MAKE) merge filter restrict_all make_lib_table
	@echo "\n######### 3. refining assembly"
	$(MAKE) refine
	@echo "\n######### 4.inferring genomes"
	$(MAKE) cov_matrix metaBAT post_metaBAT refine_bins
	@echo "\n######### 5. inferring sites"
	$(MAKE) sites
	@echo "\n######### 6. output coverage trajectories"
	$(MAKE) bin_trajectory site_trajectory
	@echo "\n##########################################"
	@echo "######### UNIT TEST END: SUCCESS #########"
	@echo "##########################################"
