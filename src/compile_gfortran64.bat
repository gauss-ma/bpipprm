rem @echo off
setlocal
set COMPILE_FLAGS=-fbounds-check -Wuninitialized -O2 -static
set LINK_FLAGS= -static -O2

gfortran -m64 -c %COMPILE_FLAGS% mod_struc.f90 
gfortran -m64 -c %COMPILE_FLAGS% mod_inp.f90
gfortran -m64 -c %COMPILE_FLAGS% mod_out.o 

gfortran -m64 -o bpip_gauss.exe %LINK_FLAGS% mod_struc.o mod_inp.o mod_out.o bpip_gauss.o 

del *.o
del *.mod

