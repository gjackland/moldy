######################################################
#
# MARTENSITES: Makefile for the MOLDY code
# 
# Mixed MPI/OpenMP free- and fixed-format Fortran code
#
# Instructions:
#	- Choose a material module $(MATERIAL)
#	- Choose a compiler and compile options
#	- "make clean"
#	- "make"
#
# If adding source files, make sure they are included
# in the "F77SRC" or "F90SRC" variables below. Add
# the appropriate dependencies in the bottom section
# (with care). There is no need to edit the targets
# section.
#
# Author: K.D'Mellow, Dec 2008.
#
######################################################

######################################################
#
# Build directory where compilation output files will
# be stored
#
######################################################
BUILD_DIR=.


######################################################
#
# Material Choice Section (with CPP directives)
# - Choose only *one* material.
# - Choose pre-processor prefix if needed (eg xlf).
#
######################################################
# Material Choice and Location of Executable
#MATERIAL= -DIRON_CARBON
MATERIAL= -DZIRCONIUM
#MATERIAL= -DGENERIC_ATVF
#MATERIAL= -DLENNARDJONES
#MATERIAL= -DMORSE 


 EXE=    $(BUILD_DIR)/moldy.zr

#Preprocessor prefix: uncomment if using xlf
#MATERIAL= -WF,$(MATERIAL)



######################################################
#
# Compiler and Flags Section
# (use only one compiler)
#
######################################################
#FC=     pgf90
#FC=     pgf90 -fast -O4
#FC=     pgf90 -fast -O4 -Minline=levels:1 -Minfo=inline
#FC=     pgf90 -fast -Mipa=fast,inline -Minfo=inline -tp core2-64
#FC=     pgf90 -fast -Mipa=fast,inline -Minfo=inline -tp x64
#FC=     pgf90 -fast -Mipa=fast,inline -tp x64
#FC=     pgf90 -O2 -Munroll=c:1 -Mnoframe -Mlre  -Mvect=sse -Mscalarsse -Mcache_align -Mflushz
#FC=     pgf90 -O2 -Munroll=c:1 -Mnoframe -Mlre -Mscalarsse -Mcache_align -Mflushz -Mbounds
#FC=     pgf90
#FC=     mpif90 -O4
#FC=     xlf90 -O4 -qfree
#FC=     ifort -O2
#FC=     f95 -O2
#FC=     gfortran -fbounds-check
FC=     gfortran -O4  #-g -ffpe-trap=invalid,zero,overflow

# C Preprocessing Directives 
CPPFLAGS= $(MATERIAL) $(DEBUG)

#OpenMP flags
#OMPFLAGS= -fopenmpcd 
#OMPFLAGS= -mp -Minfo=mp
#OMPFLAGS= -mp
OMPFLAGS = -fopenmp # For gfortran
#OMPFLAGS -openmp # For ifort

#Debugging
#DEBUGGING= -g

#	Any kind of profiling
#PROFILING= -pg

#Any libraries to be linked
#LIB= 

FFLAGS=     $(DEBUGGING) $(PROFILING) $(OMPFLAGS) $(CPPFLAGS)
# Or if just want to preprocess
#FFLAGS = -F -g $(PROFILING) $(OMPFLAGS) $(CPPFLAGS)

LFLAGS= $(FFLAGS) $(LIB)



######################################################
#
# Source Files Section - add files as appropriate
#
######################################################
#Makefile is integrated into make dependencies, so specify its name
MF=Makefile

#Fortran 77 source files
F77SRC=

#Fortran 90 source files
F90SRC= \
	constants_m.F90 \
	particles_m.F90 \
	debug_m.F90 \
	system_m.F90 \
	neighbourlist_m.F90 \
	dynamics_m.F90 \
	io_m.F90 \
	linkcell_m.F90 \
	lookup_m.F90 \
	matrix_m.F90 \
	moldin.F90\
	params_m.F90 \
	quench_m.F90 \
	random_m.F90 \
	analysis_m.F90 \
	parrinellorahman_m.F90 \
	thermostat_m.F90 \
	timers_m.F90 \
	utilityfns_m.F90 \
	iron_carbon.F90 \
	zirconium.F90 \
	generic_atvf.F90 \
        morse.F90 \
        CdAu.F90 \
        pairpot.F90 \
	potential_m.F90


#######################################################
#
# Makefile Targets: No need to edit this section
#
#######################################################
.SUFFIXES:
.SUFFIXES: .F90 .f .mod .o

OBJ=    $(F77SRC:.f=.o) $(F90SRC:.F90=.o)
SRC=    $(F77SRC) $(F90SRC) $(INC)

.F90.o:
	$(FC) $(FFLAGS) -c $<

.F90.mod:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

all:    $(EXE)

$(EXE): $(OBJ)
	$(FC) $(LFLAGS) -o $@ $(OBJ)

$(OBJ): $(MF) $(INC)

tar:
	tar cvf $(EXE).tar $(MF) $(SRC)

clean:
	rm -f $(OBJ) $(EXE) *.oo *.ipa *.ipo *.mod core



######################################################
#
# Dependencies section (manual at present).
# Modify this with care as appropriate to satisfy
# use-based dependencies within modules.
#
######################################################
moldin.o: timers_m.o thermostat_m.o random_m.o potential_m.o params_m.o  system_m.o  \
	 constants_m.o analysis_m.o parrinellorahman_m.o  dynamics_m.o \
	neighbourlist_m.o  particles_m.o linkcell_m.o \
	quench_m.o io_m.o

thermostat_m.o: random_m.o potential_m.o params_m.o  system_m.o  constants_m.o \
	 parrinellorahman_m.o  particles_m.o

analysis_m.o: thermostat_m.o


random_m.o: debug_m.o particles_m.o

dynamics_m.o: debug_m.o particles_m.o  system_m.o matrix_m.o parrinellorahman_m.o thermostat_m.o neighbourlist_m.o timers_m.o

io_m.o: params_m.o random_m.o particles_m.o parrinellorahman_m.o utilityfns_m.o analysis_m.o

params_m.o: constants_m.o utilityfns_m.o

quench_m.o: params_m.o  system_m.o  constants_m.o \
	analysis_m.o parrinellorahman_m.o dynamics_m.o neighbourlist_m.o  \
	linkcell_m.o particles_m.o potential_m.o

neighbourlist_m.o:  linkcell_m.o parrinellorahman_m.o

potential_m.o: params_m.o  constants_m.o particles_m.o  lookup_m.o parrinellorahman_m.o morse.o iron_carbon.o  zirconium.o  CdAu.o generic_atvf.o  pairpot.o

timers_m.o: debug_m.o constants_m.o

lookup_m.o: constants_m.o 

parrinellorahman_m.o analysis_m.o particles_m.o linkcell_m.o: constants_m.o params_m.o

analysis_m.o: parrinellorahman_m.o system_m.o dynamics_m.o

parrinellorahman_m.o: matrix_m.o lookup_m.o

linkcell_m.o: parrinellorahman_m.o

matrix_m.o: utilityfns_m.o
