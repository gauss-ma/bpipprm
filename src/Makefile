FC=gfortran

DEBUG=#-g -pg -fbacktrace
FCFLAGS=-O2 -ffree-line-length-none $(DEBUG) -Wall #-fcheck=bounds #-mtune=native #-fcheck=all-fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal

EXE=../exe
new_code=bpipv25252.f90
v04274_code=Bpipprm_v04274/Bpipprm.for
v19191_code=Bpipprm_v19191/bpipprm_19191_drft.for

all:
	$(FC) $(new_code) $(FCFLAGS) -o ${EXE}/bpipprm.exe
v04274:
	$(FC) $(v04274_code) $(FCFLAGS) -o ${EXE}/Bpipprm.exe
v19191:
	$(FC) $(v19191_code) $(FCFLAGS) -o ${EXE}/bpipprm_v19191_drft.exe

