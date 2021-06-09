@REM                                                                    + + +
@echo off

setlocal

set COMPILE_FLAGS=/O2 /check:format /Qipo /Qprec-div- /QaxSSE2 /trace
set LINK_FLAGS=/O2 /Qipo /check:format /Qprec-div- /QaxSSE2


ifort /compile_only %COMPILE_FLAGS% mod_struc.f90
ifort /compile_only %COMPILE_FLAGS% mod_inp.f90
ifort /compile_only %COMPILE_FLAGS% mod_out.f90
ifort /compile_only %COMPILE_FLAGS% bpip_gauss.f90

ifort /exe:aermod.exe %LINK_FLAGS% mod_struc.obj mod_inp.obj mod_out.obj bpip_gauss..obj 

del *.obj
del *.mod

