FC=gfortran

FCFLAGS=-O2 -static -Wall -fcheck=all -g -fbacktrace -fbounds-check -mtune=native -ffpe-trap=invalid,zero,overflow,underflow,denormal -fbacktrace#-fcheck=bounds

EXE=../exe
new_code=bpip_gauss.f90
v04274_code=Bpipprm_v04274/Bpipprm.for
v19191_code=Bpipprm_v19191/bpipprm_19191_drft.for


all:
	$(FC) $(new_code) $(FCFLAGS) -o ${EXE}/bpipprm_gauss.exe
v04274:
	$(FC) $(v04274_code) $(FCFLAGS) -o ${EXE}/Bpipprm.exe
v19191:
	$(FC) $(v19191_code) $(FCFLAGS) -o ${EXE}/bpipprm_v19191_drft.exe

