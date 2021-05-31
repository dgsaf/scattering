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
	GPUCCXX = $(HCXX)
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

# define formatting characters for info output
NULL :=
TAB := $(NULL)  $(NULL)

#
.PHONY : allinfo
allinfo : configinfo buildinfo makecommands

# information about current configuration
.PHONY : configinfo
configinfo :
	$(info Compiler options:)
	$(info > Compiler can be selected with COMPILERTYPE=(GCC|CLANG|AOMP|CRAY))
	$(info > GPU compiler can be selected with GPUTYPE=(NVIDIA|AMD))
	$(info > Optimisation level can be selected with OPTLEVEL=(0|1|2|3))
	$(info > Profiling can be turned on/off with PROFILER=(OFF|ON))
	$(info )

# information about current build given the commands make was passed
.PHONY : buildinfo
buildinfo :
	$(info Current compilers selected:)
	$(info > Compiling with ($(CC)|$(CXX)|$(FORT)) for CPU focused codes)
	$(info > Compiling with ($(MPICC)|$(MPICXX)|$(MPIFORT)) \
	for MPI-CPU focused codes)
	$(info > Compiling with ($(GPUCC)|$(GPUCXX)) for GPU focused codes)
	$(info > Compiling with ($(OMPCC)|$(OMPCXX)|$(OMPFORT)) \
	for OpenMP directive GPU focused codes)
	$(info > Compiling with ($(OACCCC)|$(OACCCXX)|$(OACCFORT)) \
	for OpenACC directive GPU focused codes)
	$(info )

# information about current make commands available
.PHONY : makecommands
makecommands :
	$(info Make commands:)
	$(info > Make is configured so that the following can be compiled \
	if provided this argument:)
	$(info )

#
.PHONY : dirs
dirs :
	[ -d obj ] || mkdir obj
	[ -d mod ] || mkdir mod
	[ -d bin ] || mkdir bin

#
.PHONY : clean
clean :
	rm -f obj/*
	rm -f mod/*
	rm -f bin/*

# source files
SRCS := $(wildcard src/*.f*)

# general make commands for (f|f90|f95|f03|f08) files
# # src/*.f
# obj/%.o : src/%.f
# 	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $< -o $@ -J mod/

# # src/*.f90
# obj/%.o : src/%.f90
# 	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $< -o $@ -J mod/

# src/*.(f*)
obj/%.o : src/%.$(wildcard f*)
	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $< -o $@ -J mod/

# specific make target for file with dependency
# obj/example1.o : src/example1.f90 obj/example2.o
# 	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c $^ -o $@ -J mod/

# # bin/*
# bin/main : src/main.f90 \
# 	obj/rsg.o obj/intp.o obj/wigner.o obj/io.o obj/integrate.o obj/laguerre.o

# 	$(FORT) $(COMMONFLAGS) $(FFLAGS) -c src/potential_curves.f90 \
# 	-o obj/potential_curves.o -I mod/

# 	$(FORT) $(COMMONFLAGS) $(FFLAGS) -o bin/potential_curves \
# 	obj/potential_curves.o obj/rsg.o obj/wigner.o obj/intp.o \
# 	obj/io.o obj/integrate.o obj/laguerre.o
