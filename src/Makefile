FC=gfortran
#FCFLAGS=-c -O2 -static -Wall -fcheck=all -g -fbacktrace -fbounds-check -mtune=native -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace#-fcheck=bounds
FCFLAGS=-O2 -static -Wall -fcheck=all -g -fbacktrace -fbounds-check -mtune=native -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace#-fcheck=bounds
#FCFLAGS=-fbounds-check -Wuninitialized -O2 -static
#LDFLAGS= -static -O2
#OBJECTS=mod_struc.o mod_inp.o mod_out.o bpip_gauss.o

EXE=../exe/bpip_gauss.exe
MAIN_SCRIPT=bpip_gauss.f90

all:
	$(FC) $(MAIN_SCRIPT) $(FCFLAGS) -o ${EXE}

#all: $(OBJECTS)
#	$(FC) $(LDFLAGS) $(OBJECTS) -o ${EXE}
#
#%.o: %.f90
#	$(FC) $(FCFLAGS) $< -o $@
#
#clean:
#	rm -rf *.o *.mod *.exe
