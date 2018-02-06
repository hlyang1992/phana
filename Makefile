.SUFFIXES : .o .cpp
# compiler and flags
CC     = icpc
LINK   = $(CC) -static
# LINK   = $(CC)
CFLAGS = -O3 $(DEBUG) $(UFLAG)
#
OFLAGS = -O3 $(DEBUG)
INC    = $(MKL_INCS) $(TCINC) $(SPGINC)
LIB    = $(MKL_LIBS) $(TCLIB) $(SPGLIB)

# share 
# MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
# static 
MKL_LIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

MKL_INCS = -I$(MKLROOT)/include
#
# cLapack library needed
# LPKINC = -I/opt/libs/clapack/3.2.1/include
# LPKLIB = -L/opt/libs/clapack/3.2.1/lib -lclapack -lblas -lf2c #-lm
#
# Tricubic library needed
TCINC = -I/home/hlyang/root/opt/tricubic/include
# TCLIB = -L/home/hlyang/root/opt/tricubic/lib -ltricubic
TCLIB = /home/hlyang/root/opt/tricubic/lib/libtricubic.a
#
# spglib 1.8.2, used to get the irreducible q-points
# if UFLAG is not set, spglib won't be used.

# UFLAG  = -DUseSPG
# SPGINC = -I/opt/libs/spglib/1.8.2/include
# SPGLIB = -L/opt/libs/spglib/1.8.2/lib -lsymspg

# if spglib other than version 1.8.2 is used, please 
# modify file phonon.cpp, instruction can be found by searching 1.8.2

# Debug flags
# DEBUG = -g -DDEBUG
#====================================================================
ROOT   = phana
# executable name
EXE    = $(ROOT)
#====================================================================
# source and rules
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

#====================================================================
all:  ver ${EXE}

${EXE}: $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${ROOT}.tar; tar -czvf ${ROOT}.tar.gz *.cpp  *.h Makefile README

ver:
	@echo "#define VERSION `git log|grep '^commit'|wc -l`" > version.h

#====================================================================
.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
