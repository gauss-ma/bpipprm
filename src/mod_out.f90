module writeBPIP
        use struc
contains
subroutine writeOUT(oT,sT,title,outputFileName)
        implicit none

        integer ::i
        character(78),  intent(in) :: title
        type(outTable), intent(in) :: ot(:)
        type(stkTable), intent(in) :: st(:)
        character(24),intent(in) :: outputFileName
        

        open(2,file=outputFileName,action="WRITE")
        
        WRITE (2,'(1X,A78,/)') TITLE

        !DATE:
        call writeDATE()
        WRITE (2,'(1X,A78,/)') TITLE


        WRITE(2,*) '============================'
        WRITE(2,*) 'BPIP PROCESSING INFORMATION:'
        WRITE(2,*) '============================'
        !
        WRITE(2,"(/3X,'The ',A2,' flag has been set for preparing downwash',' related data',10X)") 'P '
        WRITE(2,"('          for a model run utilizing the PRIME algorithm.',/)")
        WRITE(2,"(3X,'Inputs entered in ',A10,' will be converted to ','meters using ')") "METERS    "
        WRITE(2,"(3X,' a conversion factor of',F10.4,'.  Output will be in meters.',/)") 1.00
        WRITE(2,"(3X,'UTMP is set to ',A4,'.  The input is assumed to be in',' a local')") "UTMN"
        WRITE(2,"(3x,' X-Y coordinate system as opposed to a UTM',' coordinate system.')")
        WRITE(2,"(3x,' True North is in the positive Y',' direction.',/)")
        WRITE(2,"(3X,'Plant north is set to',F7.2,' degrees with respect to',' True North.  ',//)") 0.00
        !
        WRITE (2,'(1X,A78,///)') TITLE

        !STACK RESULTS
        WRITE(2,"(16X,'PRELIMINARY* GEP STACK HEIGHT RESULTS TABLE')")
        WRITE(2,"(13X,'            (Output Units: meters)',/)")
        WRITE(2,"(8X,'                    Stack-Building            Preliminary*')")
        WRITE(2,"(8X,' Stack    Stack     Base Elevation    GEP**   GEP Stack')")
        WRITE(2,"(8X,' Name     Height    Differences       EQN1    Height Value',//)")
        do i=1,size(sT,1),1
                WRITE(2,'(8X, A8, 4(F8.2,5X))') st(i)%stkName, st(i)%stkHeight, st(i)%BaseElevDiff, st(i)%GEPEQN1, st(i)%GEPSHV
                !                                STKN(S),        SH(S),           DIF,                 GEP(S),       PV
        end do
        WRITE(2,"(/,'   * Results are based on Determinants 1 & 2 on pages 1',' & 2 of the GEP')   ")
        WRITE(2,"( '     Technical Support Document.  Determinant',' 3 may be investigated for')  ")
        WRITE(2,"( '     additional stack height cred','it.  Final values result after')          ")
        WRITE(2,"( '     Determinant 3 has been ta','ken into consideration.')                    ")
        WRITE(2,"( '  ** Results were derived from Equation 1 on page 6 of GEP Tech','nical')     ")
        WRITE(2,"( '     Support Document.  Values have been adjusted for a','ny stack-building') ")
        WRITE(2,"( '     base elevation differences.',/)                                          ")
        WRITE(2,"( '     Note:  Criteria for determining stack heights for modeling',' emission') ")
        WRITE(2,"( '     limitations for a source can be found in Table 3.1 of the')              ")
        WRITE(2,"( '     GEP Technical Support Document.')                                        ")
        WRITE(2,"(/,/,/,/)")

        !DATE (AGAIN)
        call writeDATE()
        WRITE (2,'(//,1X,A78,/)') TITLE
        WRITE (2, *) ' BPIP output is in meters'
        
        !MAIN OUTPUT:
        do i=1,size(oT,1),1
                write(2,'(/)')
                !print*,names(i)
                !!WDIR        
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,1:6  ,1)
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,7:12 ,1)
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,13:18,1)
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,19:24,1)
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,25:30,1)
                !write(2,'(5X,"SO WDIR     ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(i,31:36,1)
                !HGT
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(1:6  ,2)
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(7:12 ,2)
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(13:18,2)
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(19:24,2)
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(25:30,2)
                write(2,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(31:36,2)
                !WID                                                         
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(1:6  ,3)
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(7:12 ,3)
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(13:18,3)
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(19:24,3)
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(25:30,3)
                write(2,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(31:36,3)
                !LEN                                                         
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(1:6  ,4)
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(7:12 ,4)
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(13:18,4)
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(19:24,4)
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(25:30,4)
                write(2,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(31:36,4)
                !XBADJ12                                                     
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(1:6  ,5)
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(7:12 ,5)
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(13:18,5)
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(19:24,5)
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(25:30,5)
                write(2,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(31:36,5)
                !YBADJ12                                                     
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(1:6  ,6)
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(7:12 ,6)
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(13:18,6)
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(19:24,6)
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(25:30,6)
                write(2,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(31:36,6)

        end do

        close(2)!cierro bpip.out
end subroutine

subroutine writeDATE()
        implicit none
        integer :: date_time(8)
        integer :: iyr,imon,iday,ihr,imin,isec
        character(len=12) :: real_clock(3)

        CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), DATE_TIME)
        IYR = DATE_TIME(1); IMON = DATE_TIME(2); IDAY = DATE_TIME(3)
        IHR = DATE_TIME(5); IMIN = DATE_TIME(6); ISEC = DATE_TIME(7)

        !header:
         WRITE (2,'(30X,"BPIP (Dated: 04274)")')
         WRITE (2,'(1X, "DATE : ",I2,"/",I2,"/",I4)') IMON, IDAY, IYR
         WRITE (2,'(1X, "TIME : ",I2,":",I2,":",I2)') IHR, IMIN, ISEC
 end subroutine

end module
