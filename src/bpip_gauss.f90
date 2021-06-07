program bpip_gauss
!MODULES-------------------------------------------------------------
use struc            !objetos
use readBPIP         !lectura de datos de entrada
use writeBPIP        !escritura de salidas

!VARIABLES-----------------------------------------------------------
implicit none
integer i,j,k,l,imax !indices
double precision, parameter :: pi=3.141593

!global
character(24) :: inputFileName="bpip.inp"
TYPE(build), allocatable :: B(:) !array de edificios
TYPE(stack), allocatable :: S(:) !array de stacks
!double precision :: R(2,2)      !matriz de rotacion
double precision :: wdir        !direccion del viento

!output vars:
real :: bh,bl,bw,xbadj,ybadj    !building hgt len width xbadj ybadj
real, allocatable :: TABLA(:,:,:) !(36,6)!tabla final  c/stack y wdir: stackname,buildhgt,buildwid,buildlen,xbadj,ybadj,wdir

!INPUT---------------------------------------------------------------
call readINP(inputFileName,B,S)   !read file y guardar datos en B y S

allocate(tabla(size(S),36,6))     !allocatar tabla de salida

!MAIN----------------------------------------------------------------
DO i=1,size(S),1        !for each stack
DO l=1,36,1             !for each wdir (cada 10 grados)
        wdir=(360-l*10)*2*pi/360.0      !wdir en radianes
        call proyectar(B,S(i),wdir) !proyecta y setea todos los valores de tiers y stacks
        
        DO j=1,size(B),1

            DO k=1,size(B(j)%t),1
                
                !check si hay otra structura cercana (de mismo tier) y combinarlos:
                !


                !check si tier afecta Stack, sino todo a cero
                !print*, "min-DIST: ",minDist( S(i)%xy2, B(j)%T(k)%xy2 )
                if( minDist( S(i)%xy2, B(j)%T(k)%xy2 ) .GT. 5.0 * B(j)%T(k)%L ) then
                   B(j)%T(k)%wid  =0.0
                   B(j)%T(k)%len  =0.0
                   B(j)%T(k)%hgt  =0.0
                   B(j)%T(k)%xbadj=0.0
                   B(j)%T(k)%ybadj=0.0
                end if


            END DO!tiers
                
            !Elijo el tier de máximo GEP Stack Height
            imax=maxloc(B(j)%T(:)%gsh,1)
 
            bw=B(j)%T(imax)%wid
            bl=B(j)%T(imax)%len  !proyected length
            bh=B(j)%T(imax)%hgt
            xbadj=B(j)%T(imax)%xbadj
            ybadj=B(j)%T(imax)%ybadj

        END DO!buildings

        tabla(i,l,1:6)=(/sngl(wdir),bh,bw,bl,xbadj,ybadj/) !copiar en tabla:

END DO!wdir
END DO!stacks

!OUTPUT:-------------------------------------------------------------
call writeOUT(tabla,S(:)%nombre)


contains


!!Calc dist build stack (proyectados)
!function DISTLIN()
!endfunction

!Calculo de GEP Stack Height, xbadj, ybadj
subroutine calcGepStackHeight(Bz0,t,s)!xy,xy2,h,
        implicit none
        real :: Bz0
        type(tier),intent(inout) :: T
        type(stack), intent(inout) :: S

        !calculo widht height L, gepStackHeight para este tier
        T%gsh= Bz0 + t%hgt - s%z0 + 1.5*t%L        !Equation 1 (GEP, page 6)

        T%xbadj= T%ymin - s%xy2(2)                 !XBADJ = YMIN(C) - YPSTK                 !ymin - ystack
        T%ybadj= s%xy2(1) - (T%xmin + T%wid*0.5)   !YBADJ = XPSTK - (XMIN(C) + TW * 0.5)    !xstack- xmin + tw*0.5

end subroutine


!Calculo de nuevas coordenadas
subroutine proyectar(B,S,wdir)
        implicit none
        type(build),intent(inout)   :: B(:)
        type(stack), intent(inout)  :: S
        double precision,intent(in) :: wdir
        double precision :: R(2,2)
        integer :: i,j,k
        
        !matriz de rotación:       
        R(1,1)=dcos(wdir);  R(2,2)=R(1,1)
        R(1,2)=dsin(wdir);  R(2,1)=-1*R(1,2)

        !stack
        S%xy2=sngl(matmul(R,dble(S%xy)))  !proyected coords stack

        !build
        do i=1,size(B),1
                do j=1,size(B(i)%t),1
                        !do k=1,size(B(i)%t(j)%xy,1),1
                        !        B(i)%t(j)%xy2(k,1:2)=sngl(matmul(R,dble(B(i)%t(j)%xy(k,1:2))))  !proyected coords tier
                        !enddo
                        B(i)%t(j)%xy2(:,1:2)=transpose( sngl( matmul( R, dble(transpose( B(i)%t(j)%xy(:,1:2) ) ) ) ) ) !proyected coords tier

                        call setTierValues(B(i)%T(j))
                        call calcGepStackHeight(B(i)%z0,B(i)%T(j),S)
                enddo
        enddo
end subroutine

subroutine setTierValues(T)
        implicit none
        type(tier),intent(inout)   :: T

        T%xmin=minval( T%xy2(:,1) )
        T%xmax=maxval( T%xy2(:,1) )
        T%ymin=minval( T%xy2(:,2) )
        T%ymax=maxval( T%xy2(:,2) )
                                            
        T%wid=T%xmax - T%xmin
        T%len=T%ymax - T%ymin
        !T%hgt=

        T%L=min(T%wid, T%hgt)
end subroutine


function minDist(xy,xxyy) result(dist)
        implicit none
        real,intent(in) :: xy(2)
        real,intent(in) :: xxyy(:,:)
        real :: dist

        dist=sqrt( minval( (xxyy(:,1)-xy(1))**2 + (xxyy(:,2)- xy(2))**2 ) )
endfunction


!  SUBROUTINE CNRLIN (XI,YI,X1, Y1, X2, Y2, BET, DIST, XKP, YKP)
!       calculate corner perpendicular to side distance,
!        intercept point, and determine if intercept on or between
!        corners.

!  SUBROUTINE DISLIN (X1, Y1, X2, Y2, L5, IBET, XSP, YSP)
!     calculate if stack directly downwind of a side and on or
!      within 5L of side.
!  L is used in the subroutine DISLIN, to determine if a stack
!  is at or within 5L directly downwind of one of the sides of each
!  structure. If the stack is at or within 5L directly downwind of
!  the structure, the flag IBET is set.

!  SUBROUTINE MXBWH(MB, MT, MBT, MSK, MD, DPADX, DPADY, DFLG,
!*                 DE, DHWE, DPBH, DPBL, DPBW, BELEV, SB,
!*                 GEP, GEPBH, GEPBW, GEPIN, MPADX, MPADY,
!*                 MHWE, MXPBH, MXPBL, MXPBW,
!*                 MI, MJ, TNUM2, TLIST2, MTNUM, MTLIST,
!*                 D, I, S, C, TW, HTA, WS, TL1, IG, BL,
!+                 XBADJ,YBADJ)

! SUBROUTINE GPC (MB, MBT, MT, MSK, BELEV,
!*                SB, GEP,GEPBH,GEPBW,GEPIN, TNUM2, TLIST2,
!*                GTNUM, GTLIST, GDIRS, MI, MJ,
!*                D, I, C, S, TW, WS, HTA, CH, IG)
!    calculate GEP values

end program
