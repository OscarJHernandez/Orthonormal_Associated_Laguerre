# This is the makefile program for compiling the Eigen Value solver

# The compiler
FC = gfortran -g -fbounds-check -O2 
#FC = ifort -g -parallel -traceback -O1 -heap-arrays
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
FCFLAGS += -Wall -fbounds-check

# libraries needed for linking
LDFLAGS = -llapack

# List of executables to be built within the package
PROGRAMS = quadrature.f mod_integrate.f90 mod_Laguerre.f90
# The List of objects that will be built
_OBJS=  quadrature.o mod_integrate.o mod_Laguerre.o
OBJS = $(patsubst %,Objects/%,$(_OBJS))
#PROG= gen_H_FB
PROG= main

# Instructions for building all Objects needed:
Objects/%.o: Programs/%.f90
	$(FC)  -o  $@  -c  $<

Objects/%.o: Programs/%.f
	$(FC)  -o  $@  -c  $<

# Instructions for building the final executable program
$(PROG): $(OBJS)
	$(FC) $(FCFLAGS) $(OBJS) -o $(PROG) $(PROG).f90 $(LDFLAGS)

clean:
	rm -f Objects/*.o *.mod *.MOD

super_clean:
	rm -f Objects/*.o *.mod *.MOD
	rm -f Programs/*.o *.mod *.MOD
	rm -f QN_Excited*
	rm -f QN_Ground*
	rm -f Eigen*
	rm -f H_out*
	rm -f Matrix*
	
