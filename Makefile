# header


# makefile configuration settings

# make flags, with optional default values
# COMPILERTYPE = (GCC|CLANG|AOMP|CRAY|INTEL)
# GPUTYPE = (NVIDIA|AMD)
# OPTLEVEL = (0|1|2|3)
# PROFILER = (OFF|ON)
COMPILERTYPE ?= GCC
GPUTYPE ?= NVIDIA
OPTLEVEL ?= 2
PROFILER ?= OFF

# set optmisation flags
OPTFLAGS = -O$(OPTLEVEL)

# set profiling flags
ifeq ($(PROFILER), ON)
	OPTFLAGS += -pg -g
endif


# compiler selection

# define c compilers
GCC = gcc
CLANG = clang
AOMP = aompcc
CRAYCC = cc
INTELCC = icc

# define c++ compilers
GCCCXX = g++
CLANGCXX = clang++
AOMPCXX = aompcc
CRAYCXX = CC
INTELCXX = icpc

# define fortran compilers
GCCFORT = gfortran
CLANGFORT = flang
CRAYFORT = ftn
INTELFORT = ifort

# c compilation flags
GCCCFLAGS = -std=c11

# fortran compilation flags
GCCFFLAGS = -cpp  -dM -ffixed-line-length-none \
	-Wall -Wextra -Wconversion -pedantic -fimplicit-none -fcheck=all
INTELFFLAGS = -cpp -extend-source -D_INTELFTN

# define cuda compilers
NCC = nvcc
NCXX = nvcc
CUDA_FLAGS = -DUSECUDA

# define hip compilers
HCC = hipcc
HCXX = hipcc
HIP_FLAGS = -DUSEHIP
OCL_FLAGS = -lOpenCL -DUSEOPENCL

# define mpi compilers
MPICC = mpicc
MPICXX = mpicxx
MPIFORT = mpif90

# openmp compilation flags
GCCOMP_FLAGS = -fopenmp -DUSEOPENMP
GCCOMPTARGET_FLAGS = -fopenmp -DUSEOPENMPTARGET
INTELOMP_FLAGS = -qopenmp -DUSEOPENMP
INTELMPTARGET_FLAGS = -qopenmp -DUSEOPENMPTARGET
GCCOACC_FLAGS = -fopenacc -fopt-info-optimized-omp -DUSEOPENACC
INTELOACC_FLAGS = -qopenacc -fopt-info-optimized-omp -DUSEOPENACC

# select default compilers, but change compilers if specified
CC = $(GCC)
CXX = $(GCCCXX)
FORT = $(GCCFORT)
GPUCC = $(NCC)
GPUCXX = $(NCXX)
FFLAGS = $(GCCFFLAGS)
OMP_FLAGS = $(GCCOMP_FLAGS)
OMPTARGET_FLAGS = $(GCCOMPTARGET_FLAGS)
OACC_FLAGS = $(GCCOACC_FLAGS)
CFLAGS = $(GCCCFLAGS)

ifeq ($(COMPILERTYPE), CLANG)
	CC = $(CLANG)
	CXX = $(CLANGCXX)
	FORT = $(CLANGFORT)
endif

ifeq ($(COMPILERTYPE), AOMP)
	CC = $(AOMP)
	CXX = $(AOMPCXX)
endif

ifeq ($(COMPILERTYPE), CRAY)
	CC = $(CRAYCC)
	CXX = $(CRAYCXX)
	FORT = $(CRAYFORT)
	FFLAGS = -eZ -ffree
endif

ifeq ($(COMPILERTYPE), INTEL)
	CC = $(INTELCC)
	CXX = $(INTELCXX)
	FORT = $(INTELFORT)
	FFLAGS = $(INTELFFLAGS)
	OMP_FLAGS = $(INTELOMP_FLAGS)
	OMPTARGET_FLAGS = $(INTELOMPTARGET_FLAGS)
endif

ifeq ($(GPUTYPE), AMD)
	GPUCC = $(HCC)
	GPUCXX = $(HCXX)
endif

# select compiler used for openmp
OMPCC = $(CC)
OMPCXX = $(CXX)
OMPFORT = $(FORT)

# select compilers used for openacc
OACCCC = $(CC)
OACCCXX = $(CXX)
OACCFORT = $(FORT)

# select compilers used for opencl
OCLC = $(CCGPU)
OCLCXX = $(CXXGPU)

# flags common to all c, c++, fortran compilers
COMMONFLAGS = $(OPTFLAGS)


# file expansions and related variables

# define source files (from contents of src/)
SRCS := $(wildcard src/*.f*) $(wildcard src/*.F*)

# define the base-names and file-suffixes of the source files
BSNS := $(notdir $(basename $(SRCS)))
SUFS := $(suffix $(SRCS))

# define the object and module files for each source file
OBJS := $(addprefix obj/,$(addsuffix .o, $(BSNS)))
MODS := $(addprefix mod/,$(addsuffix .mod, $(BSNS)))

# define binary targets to be made
BINS :=


# make commands

# make obj/, mod/, bin/ directories if they dont already exist
.PHONY : dirs
dirs :
	[ -d obj ] || mkdir obj
	[ -d mod ] || mkdir mod
	[ -d bin ] || mkdir bin

# remove contents of obj/, mod/, bin/ directories
.PHONY : clean
clean :
	rm -f obj/*
	rm -f mod/*
	rm -f bin/*

# implicit rule for arbitrary fortan targets
# 	obj/%.o : src/%.f
# 	obj/%.o : src/%.f90
# 	obj/%.o : src/%.f95
# 	...
obj/%.o : $(firstword $(addprefix src/%,$(SUFS)))
	@echo "> $@ : $^"
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $< -o $@ -J mod/

# explicit target dependencies
obj/parameters.o : src/debug.h
obj/basis.o : src/debug.h obj/parameters.o

# implicit rule for binary targets
bin/% : $(OBJS)
	@echo "> $@ : $^"
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -o $@ $(OBJS)


# info commands

# define formatting characters for info output
NULL :=
TAB := $(NULL)  $(NULL)

# all information targets
.PHONY : info_all
info_all : info_config info_build info_make info_dirs info_files

# information about makefile configuration
.PHONY : info_config
info_config :
	@echo "Makefile configuration options (with defaults):"
	@echo "> COMPILERTYPE ?= [GCC] | CLANG | AOMP | CRAY"
	@echo "> GPUTYPE      ?= [NVIDIA] | AMD"
	@echo "> OPTLEVEL     ?= 0 | 1 | [2] | 3"
	@echo "> PROFILER     ?= [OFF] | ON"

# information about current build given the commands make was passed
.PHONY : info_build
info_build :
	@echo "Compilers currently selected for (C, C++, FORTRAN) codes:"
	@echo "> CPU         := ($(CC), $(CXX), $(FORT))"
	@echo "> MPI         := ($(MPICC), $(MPICXX), $(MPIFORT))"
	@echo "> GPU         := ($(GPUCC), $(GPUCXX), )"
	@echo "> GPU-OpenMP  := ($(OMPCC), $(OMPCXX), $(OMPFORT))"
	@echo "> GPU-OpenACC := ($(OACCCC), $(OACCCXX), $(OACCFORT))"

# information about contents of src/, obj/, mod/, bin/ directories
.PHONY : info_dirs
info_dirs :
	@echo "Directory contents:"
	@echo "> src/* = $(wildcard src/*)"
	@echo "> obj/* = $(wildcard obj/*)"
	@echo "> mod/* = $(wildcard mod/*)"
	@echo "> bin/* = $(wildcard bin/*)"

# information about source, base name, and object files
.PHONY : info_files
info_files :
	@echo "Source files (and derived files):"
	@echo "> SRCS := $(SRCS)"
	@echo "> BSNS := $(BSNS)"
	@echo "> SUFS := $(SUFS)"
	@echo "> OBJS := $(OBJS)"
	@echo "> MODS := $(MODS)"
	@echo "> BINS := $(BINS)"

# information about current make commands available
.PHONY : info_make
info_make :
	@echo "Make commands:"
	@echo "> clean       : Removes contents of obj/, mod/, bin/."
	@echo "> dirs        : Creates obj/, mod/, bin/ if they dont exist."
	@echo "> info_all    : All info_* commands."
	@echo "> info_config : Configuration options for makefile."
	@echo "> info_build  : The current compilers selected in makefile."
	@echo "> info_dirs   : The contents of relevant directories."
	@echo "> info_files  : The file expansions used in makefile."
	@echo "> info_make   : The make commands available."
