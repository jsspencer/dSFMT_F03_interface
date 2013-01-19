SHELL=/bin/bash

VPATH = ../dSFMT-src-2.2:../dSFMT_F03_interface
CPATH = ../dSFMT-src-2.2

CC = gcc
CFLAGS = -O3 -I $(CPATH)

FC = gfortran
FFLAGS = -O3

CPPFLAGS = -DDSFMT_MEXP=19937

.PHONY:
.PHONY: tests clean

.SUFFIXES:
.SUFFIXES: .c .f90 .F90

# Run tests
tests: testc.x testf.x
	./testc.x > testc.out
	./testf.x > testf.out
	diff -w testc.out testf.out

# Programs
testc.x: dSFMT.o testc.o
	$(CC) $^ -o $@

testf.x: dSFMT.o dSFMT_utils.o dSFMT_interface.o testf.o
	$(FC) $^ -o $@

# Tidy
clean:
	rm -f *.o *.mod *.x *.out

# Objects
%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Dependencies (not complete!)
dSFMT.o: dSFMT.h dSFMT-params11213.h dSFMT-params1279.h dSFMT-params132049.h dSFMT-params19937.h dSFMT-params216091.h dSFMT-params2203.h dSFMT-params4253.h dSFMT-params44497.h dSFMT-params521.h dSFMT-params86243.h dSFMT-params.h
dSFMT_utils.o testc.o: dSFMT.h
testf.o: dSFMT_interface.o
