module writeBPIP
contains
subroutine writeOUT(tabla,names)
        implicit none

        integer ::i
        character(8) :: names(:)
        real, intent(in) :: tabla(:,:,:) !tabla final  buildhgt,buildwid,buildlen,xbadj,ybadj

        do i=1,size(tabla,1),1
                print*,names(i)
                !WDIR        
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,1:6  ,1)
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,7:12 ,1)
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,13:18,1)
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,19:24,1)
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,25:30,1)
                write(*,'(5X,"SO WDIR     ",a8,6(f8.2))') names(i),tabla(i,31:36,1)
                write(*,"(/)")
                !HGT
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,1:6  ,2)
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,7:12 ,2)
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,13:18,2)
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,19:24,2)
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,25:30,2)
                write(*,'(5X,"SO BUILDHGT ",a8,6(f8.2))') names(i),tabla(i,31:36,2)
                write(*,"(/)")
                !WID
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,1:6  ,3)
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,7:12 ,3)
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,13:18,3)
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,19:24,3)
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,25:30,3)
                write(*,'(5X,"SO BUILDWID ",a8,6(f8.2))') names(i),tabla(i,31:36,3)
                write(*,"(/)")
                !LEN
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,1:6  ,4)
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,7:12 ,4)
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,13:18,4)
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,19:24,4)
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,25:30,4)
                write(*,'(5X,"SO BUILDLEN ",a8,6(f8.2))') names(i),tabla(i,31:36,4)
                write(*,"(/)")
                !XBADJ
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,1:6  ,5)
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,7:12 ,5)
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,13:18,5)
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,19:24,5)
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,25:30,5)
                write(*,'(5X,"SO XBADJ    ",a8,6(f8.2))') names(i),tabla(i,31:36,5)
                write(*,"(/)")
                !YBADJ
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,1:6  ,6)
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,7:12 ,6)
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,13:18,6)
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,19:24,6)
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,25:30,6)
                write(*,'(5X,"SO YBADJ    ",a8,6(f8.2))') names(i),tabla(i,31:36,6)
                write(*,"(/)")
        
        end do
end subroutine
end module


!Asi lo tenia Don. EPA
!C             OF ALL TIERS, PRINT/SAVE WHICH HAS MOST EFFECT BY STACK AND WD
!        IF (SWT .EQ. 0 .OR. SWT .EQ. 2) THEN
!            WRITE (12, 461)
!            WRITE (12, 462) IMON, IDAY, IYR
!            WRITE (12, 463) IHR, IMIN, ISEC
!            WRITE (12, 297)
!            WRITE (12, 1) TITLE
!            WRITE (12, *) ' BPIP output is in meters'
!
!          DO 510 S = 1, NS
!              L = NDIR / 6
!               WRITE(12,297)
!            DO I = 1, 6
!              J = (I-1) * 6 + 1
!              K = I * 6
!              WRITE (12,293) STKN(S), (MXPBH(S,D) , D = J,K)
!            END DO
!            DO I = 1, 6
!              J = (I-1) * 6 + 1
!              K = I * 6
!              WRITE (12,296) STKN(S), (MXPBW(S,D) , D = J,K)
!            END DO
!CVRT
!            IF (SWT .EQ. 0) THEN
!              DO I = 1, 6
!                J = (I-1) * 6 + 1
!                K = I * 6
!                WRITE (12,299) STKN(S), (MXPBL(S,D) , D = J,K)
!              END DO
!
!              DO I = 1, 6
!                J = (I-1) * 6 + 1
!                K = I * 6
!                WRITE (12,290) STKN(S), (MPADX(S,D) , D = J,K)
!              END DO
!              DO I = 1, 6
!                J = (I-1) * 6 + 1
!                K = I * 6
!                WRITE (12,291) STKN(S), (MPADY(S,D) , D = J,K)
!              END DO
!            END IF
!
!290   FORMAT(5X,'SO XBADJ    ', A8, 6F8.2)
!291   FORMAT(5X,'SO YBADJ    ', A8, 6F8.2)
!292   FORMAT(3(/1X, 8F6.2))
!293   FORMAT(5X,'SO BUILDHGT ', A8, 6F8.2)
!296   FORMAT(5X,'SO BUILDWID ', A8, 6F8.2)
!297   FORMAT(/)
!299   FORMAT(5X,'SO BUILDLEN ', A8, 6F8.2)

