FC=gfortran
FCFLAGS=-c -O2 -static -Wall -fcheck=all -g -fbacktrace -fbounds-check -mtune=native -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace#-fcheck=bounds
#FCFLAGS=-fbounds-check -Wuninitialized -O2 -static
LDFLAGS= -static -O2
EXE=../exe/bpip_gauss.exe

OBJECTS=mod_struc.o mod_inp.o mod_out.o bpip_gauss.o

all: $(OBJECTS)
	$(FC) $(LDFLAGS) $(OBJECTS) -o ${EXE}

%.o: %.f90
	$(FC) $(FCFLAGS) $< -o $@

clean:
	rm -rf *.o *.mod *.exe

