program bpip_gauss
implicit none
!OBJETCS/TYPES ------------------------------------------------------
type tier
  double precision,allocatable :: xy(:,:)   !coords
  real :: h                                 !height  (leido del .inp) (NO CAMBIA)
  !projected:
  double precision,allocatable :: xy2(:,:)  !coords (cambia para c/wdir)
  real :: xmin,xmax,ymin,ymax  !boundaries          (cambia para c/wdir)
  real :: hgt,wid,len          !height,width,length (cambia para c/wdir)
  real :: L                    !L: min{ wid , hgt } (cambia para c/wdir)
  real :: gsh                  !gep stack height    (cambia para c/stack y wdir)
  real :: xbadj,ybadj          !gep stack height    (cambia para c/stack y wdir)
  real, allocatable :: minDist(:,:)
endtype

type building
  character(8)           :: nombre   !name
  real                   :: z0       !base height
  type(tier),allocatable :: t(:)     !tier
endtype

type stack
  character(8)     :: nombre
  real             :: z0, h
  real             :: xy(2), xy2(2)   !original
  !double precision :: xy(2), xy2(2)
  logical          :: affected_by_tier=.false.
  logical          :: onRoof=.false.
  !character(8)     :: whichRoof='' !es más seguro por indice que por nombre
  integer          :: whichRoof(2)=[0,0]
endtype

type outTable  !output table
  character(8) ::stkName
  real         :: tabla(36,6)  !36 windirs, 6vars: wdir,hgt,wid,len,xbadj,ybadj
end type

type stkTable  !stack's table
  character(8) :: stkName
  real         :: stkHeight, BaseElevDiff, GEPEQN1=0.0, GEPSHV
end type

!VARIABLES-----------------------------------------------------------
integer :: i,j,k,d

!global params:
double precision, parameter :: pi=3.141593 !3.141592653589793_8
double precision, parameter :: deg2rad=pi/180.0
character(24),    parameter :: inputFileName="BPIP.INP"
character(24),    parameter :: outputFileName="bpip.out"

!work variables:
TYPE(building), allocatable :: B(:) !array de edificios
TYPE(stack), allocatable    :: S(:) !array de stacks
!type(tier)                  :: fT   !focal Tier
character(78)    :: title
double precision :: wdir         !direccion del viento
logical          :: SIZ,ROOF,any_tier_affects_stack_at_this_wdir
real             :: refGSH,refWID !max GSH and WID encountred (for given stack and wdir)
integer          :: bmax,tmax     !indices where building and tier of max GSH are stored on "build" objct array

!out tables:
type(outTable), allocatable :: oTable(:)  !output table
type(stkTable), allocatable :: sTable(:)  !stack  table

!INPUT---------------------------------------------------------------
call readINP(inputFileName,B,S,title)   !read file & store data in B & S
allocate(oTable(size(S)))   !allocatar tabla de salida
allocate(sTable(size(S)))   !allocatar tabla de stacks

!MAIN----------------------------------------------------------------
call check_which_stack_over_roof(S,B)
!call check_which_tiers_to_merge(B)
!call calc_mindist_between_tiers(B)

DO d=1,36   !for each wdir (c/10 deg)
     print '("      Wind flow passing", I4," degree direction.")', d*10 
     wdir=d*10.0*deg2rad  !wdir [rad]
     !call rotateCoords(B,S,wdir)       !roto coordenadas según wdir. "BIEN"
     call rotateCoords(B,S,sngl(wdir))  !roto coordenadas según wdir. "MAL" (original)

     DO i=1,size(S)                                 !for each stack
         refGSH=0.0; refWID=0.0; tmax=0; bmax=0
         any_tier_affects_stack_at_this_wdir=.false.
         S(i)%affected_by_tier=(S(i)%affected_by_tier .OR. .false.)

         DO j=1,size(B)                             !for each building
            DO k=1,size(B(j)%T)                     !for each tier

               SIZ  = isInsideSIZ(S(i), B(j)%T(k))
               ROOF = S(i)%whichRoof(1) ==j .and.  S(i)%whichRoof(2) ==k
               if ( SIZ .OR. ROOF ) then 
                  any_tier_affects_stack_at_this_wdir=.true. 
                  S(i)%affected_by_tier=.true.

                  !call mergeCloseTiers(B(j)%T(k),B,k,B(j)%nombre)   !si hay otra structura cercana (de == tier) y combinarlos
                  call calcGepStackHeight(B(j)%z0, B(j)%T(k), S(i)) !calc: GSH, XBADJ, YBADJ
                  
                  if ( hasMaxGSH(B(j)%T(k), refWID, refGSH) ) then
                     refGSH=B(j)%T(k)%gsh; refWID=B(j)%T(k)%wid
                     bmax=j      !building max (index)
                     tmax=k      !tier     max (index)
                  end if
               endif
            END DO!tiers
         END DO!buildings
         
         if ( any_tier_affects_stack_at_this_wdir ) then
            call addToOutTable(oTable(i),S(i),B(bmax)%T(tmax),d,wdir) 
            if ( refGSH > sTable(i)%GEPEQN1 ) then
               call addToStkTable(sTable(i),S(i),B(bmax)%T(tmax),B(bmax)%z0)
            endif
         else
            call no_tier_affects_this_stack_at_this_wdir(S(i),d,wdir,oTable(i),sTable(i))
         endif
     END DO!stacks
END DO!wdir

!OUTPUT:-------------------------------------------------------------
call writeOUT(oTable,sTable,title,outputFileName)
WRITE(*,'(/,A,/)') ' END OF BPIP RUN.'

contains
!GEOMETRY   ***************************************************************************************
subroutine rotateCoords(B,S,wdir) !Calculo de nuevas coordenadas
    implicit none
    type(building),intent(inout) :: B(:)
    type(stack), intent(inout)   :: S(:)
    real,             intent(in) :: wdir  !original
    !double precision,intent(in) :: wdir
    double precision :: R(2,2)    !usan distinta precision para stacks y buildings :/
    real             :: Rs(2,2)   !usan distinta precision para stacks y buildings :/
    integer :: i,j,k
    !matriz de rotación:       
    !R(1,1)=dcos(wdir); R(1,2)=-dsin(wdir)
    R(1,1)=dcos(dble(wdir)); R(1,2)=-dsin(dble(wdir)) !original
    R(2,1)=-R(1,2)   ; R(2,2)= R(1,1) 
    Rs(1,1)=cos(wdir)  ; Rs(1,2)=-sin(wdir); Rs(2,1)=-Rs(1,2)   ; Rs(2,2)= Rs(1,1)    !original

    !stack
    do i=1,size(S)
       !S(i)%xy2=matmul(R,S(i)%xy)   !proyected coords stack
       S(i)%xy2=matmul(Rs,S(i)%xy)   !original
    enddo
    !buildings
    do i=1,size(B)
       do j=1,size(B(i)%T)
           do k=1,size(B(i)%T(j)%xy(:,1))
              B(i)%T(j)%xy2(k,:)=matmul(R,B(i)%T(j)%xy(k,:))   !proyected coords tier corners
           enddo
           call calcTierProyectedValues(B(i)%T(j) )             ! calc tier: xmin,xmax,ymin,ymax
       enddo
    enddo
end subroutine

real function DisLin2Point (X0, Y0, X1, Y1, XP, YP) 
!LOGICAL FUNCTION DisLin2Point (X0, Y0, X1, Y1, XP, YP, maxL)
   !REPLACE of original DISLIN procedure
   !computes min distance between side/line defined by (v0,v1), and point (vp).
   implicit none
   double precision, intent(in)    :: x0,x1,y0,y1,xp,yp
   double precision,dimension(2)   :: v,p,d                !side (v), point (p), and v-p vector (d)
   double precision                :: mod2v,v_dot_p,dist   !length of v, dot_product(v,p)             
   !double precision, intent(in)    :: maxL                 !max distance admisible (for exmple: L or 5L) 

   v=[x1-x0, y1-y0] !v1-v0                                 ! take vectors to same origin
   p=[xp-x0, yp-y0] !vp-v0                                 ! take vectors to same origin
   mod2v = sqrt(dot_product(v,v))
   !proj_of_p_on_v=dot_product(v,p)/mod2v !proyeccion de p en v
   v_dot_p=dot_product(v,p) 

   if ( v_dot_p > mod2v**2 .or. v_dot_p < 0.0 .or. mod2v == 0.0 ) then
       d=v-p                               !diference vector
       dist=sqrt(dot_product(d,d))         !distance between the "tip" of the vectors
   else 
       dist=abs(v(1)*p(2)-v(2)*p(1))/mod2v !distance parametric line to point
   end if
   DISLIN2POINT = dist
   !DISLIN2POINT = dist <= maxL

end function

!STRUCTURE INFLUENCE ZONE (SIZ) *******************************************************************
logical function isInsideSIZ(S,T)       result(SIZ)  
    ! check if Stack is over influence of a tier (Structure Influence Zone, SIZ)
    implicit none
    type(stack), intent(in) :: S
    type(tier), intent(in)  :: T
    double precision:: x_stack,y_stack
    integer ::i,j,n
    x_stack=dble(S%xy2(1) )
    y_stack=dble(S%xy2(2) )
    
    SIZ=           x_stack .GE. (T%xmin - 0.5*T%L)  
    SIZ=SIZ .AND. (x_stack .LE. (T%xmax + 0.5*T%L) )
    SIZ=SIZ .AND. (y_stack .GE. (T%ymin - 2.0*T%L) )
    !SIZ=SIZ .AND. (y_stack .LE. (T%ymax + 5.0*T%L) )  !old
    if ( SIZ ) then
       n=size(T%xy2(:,1))
       do i=1,n
           j=mod(i,n)+1
           SIZ = DISLIN2POINT(T%xy2(i,1), T%xy2(i,2), T%xy2(j,1),T%xy2(j,2), x_stack,y_stack) <= T%L*5.0
           if (SIZ) exit
       enddo
    endif
end function

!STACK IS OVER ROOF? ******************************************************************************
logical function isInsideTier(S,T)    result(inTier) 
    implicit none
    type(stack), intent(in) :: S
    type(tier), intent(in)  :: T
    double precision:: angle,signo=1.0
    double precision:: v1(2), v2(2), p(2)
    integer :: i,j,n!,k
    p=dble(S%xy)
    angle=0.0
    n=size(T%xy(:,1))
    do i=1,n
        j=mod(i,n)+1
        v1=T%xy(i,:)-p
        v2=T%xy(j,:)-p
        signo = dsign(signo,v1(1)*v2(2)-v1(2)*v2(1))      !sign of 3rd-component v1 x v2 (cross prod)
        angle = angle + signo * dacos( dot_product(v1,v2) / sqrt(dot_product(v1,v1)*dot_product(v2,v2))) 
    enddo
    inTier=(ABS(2*pi - ABS(angle)) .LT. 1e-4) 
end function

subroutine check_which_stack_over_roof(S,B)
    implicit none
    type(stack) ,intent(inout) :: S(:)
    type(building),intent(in)  :: B(:)        
    integer                    :: i,j,k

    print*, "Detect if a stack is on top of a roof"
    DO i=1,size(S)                                  !for each stack
       DO j=1,size(B)                               !for each building
           DO k=1,size(B(j)%T)                      !for each tier
             if ( isInsideTier(S(i), B(j)%T(k)) ) then
                print '(A10,"=>",A10)',S(i)%nombre, B(j)%nombre
                S(i)%onRoof=.true.
                S(i)%whichRoof=[j,k] 
                exit
             endif
          enddo
       enddo
    enddo
end subroutine


!BPIP PARAM CALCULATIONS **************************************************************************

subroutine calcTierProyectedValues(T) !Calculo de XMIN XMAX YMIN YMAX
    implicit none
    type(tier),intent(inout)   :: T
    T%xmin=sngl(minval(T%xy2(:,1)) ); T%xmax=sngl(maxval(T%xy2(:,1)) )
    T%ymin=sngl(minval(T%xy2(:,2)) ); T%ymax=sngl(maxval(T%xy2(:,2)) )

    T%wid=T%xmax - T%xmin
    T%len=T%ymax - T%ymin
    T%hgt=T%h
    T%L  =min(T%wid, T%hgt)
end subroutine

subroutine calcGepStackHeight(Bz0,T,S) !Calculo de GEP Stack Height, XBADJ, YBADJ
    implicit none
    real :: Bz0
    type(tier),intent(inout)   :: T
    type(stack), intent(inout) :: S
    real                       :: x_stack, y_stack
    x_stack= S%xy2(1)
    y_stack= S%xy2(2)
    T%gsh  = Bz0 + T%hgt - s%z0 + 1.5*T%L      !Equation 1 (GEP, page 6)
    T%ybadj= x_stack - (T%xmin + T%wid * 0.5)  !YBADJ = XPSTK - (XMIN(C) + TW * 0.5)
    T%xbadj= T%ymin - y_stack                  !XBADJ = YMIN(C) - YPSTK             
end subroutine

logical function hasMaxGSH(T, WID, GSH) !true si tier de máximo GSH. (si hay dos con == GSH) quedarme con el de menor width
        implicit none
        type(tier) :: T
        real       :: WID, GSH
        hasMaxGSH=GSH < T%gsh .OR. ( GSH == T%gsh .AND. WID >= T%wid) 
end function

!MERGE TIERS **************************************************************************************

!subroutine calc_mindist_between_tiers(B)
!Aca habria que ver como implementar alguna matriz o estructura que guarde las distancias minimas entre edificios para usar luego.
!   implicit none
!   type(building),intent(inout) :: B(:)        
!   integer :: i,j,ii,jj,k
!    
!
!
!    n=size(B)
!    m=
!    do i=1,,1                            !on each building
!          do j=1,size(B(i)%T)                   !on each tier
!
!allocate(B(i)%T(j)%minDist,
!    if ( B(i)%nombre /= nombre ) then
!       if ( size(B(i)%T) >= ntier ) then
!             
!             do k=1,size(T%xy2(:,1))
!                !dist = sngl( minDist( T%xy2, B(i)%T(ntier)%xy2 ) )
!                !dist=DISLIN2POINT(T%xy2(i,1), T%xy2(i,2), T%xy2(j,1),T%xy2(j,2), x_stack,y_stack)
!                dist = CNRLIN(T%xy2(i,1), T%xy2(i,2), T%xy2(j,1),T%xy2(j,2), x_stack,y_stack    )   !original
!                minL=min(t%L, B(i)%T(ntier)%L) 
!                if ( dist < minL ) then
!                   call CombineTiers(T, B(i)%T(ntier) )
!                endif
!             enddo
!          enddo
!       endif
!    endif
!    enddo
!end subroutine

real function minDist(xy1,xy2)
     double precision :: xy1(:,:), xy2(:,:)
     integer :: i,j,k,n,m
     m=size(xy1(:,1))
     n=size(xy2(:,1))
     do i=1,n
     do j=1,m
        k=mod(j,m)+1
        !minDist=min(minDist,cnrlin( xy1(j,1),xy1(j,2),  xy1(k,1),xy1(k,2), xy2(i,1), xy2(i,2)) )
        minDist=min(minDist,dislin2point( xy1(j,1),xy1(j,2),  xy1(k,1),xy1(k,2), xy2(i,1), xy2(i,2)) )
        !print*,"DIST=",mindist, dislin2point( xy1(j,1),xy1(j,2),  xy1(k,1),xy1(k,2), xy2(i,1), xy2(i,2)) 
     enddo
     enddo
end function


!Buscar tiers cercanos y combinar
subroutine mergeCloseTiers(T,B,ntier,nombre)
    implicit none
    type(tier)    ,intent(inout) :: T
    type(building),intent(inout) :: B(:)        
    integer,intent(in) :: ntier
    character(8),intent(in) :: nombre
    real :: dist, minL
    integer :: i,j,k,p

    !Cuando tenga una matriz con distancias entre edificios podria usar:
    !do i=1,size(B),1                   !on each building
    !  do j=1,size(B(i)%T)              !on each tier
    !      if (T%minDist(i,j) < T%L )
    !         combineTiers(T,B(i)%T(j))
    !      endif
    !  enddo
    !enddo
    
    do i=1,size(B),1                            !on each building
    if ( B(i)%nombre /= nombre ) then
       if ( size(B(i)%T) >= ntier ) then
          do j=1,size(B(i)%T)                   !on each tier
             
             dist = sngl( minDist( T%xy2, B(i)%T(ntier)%xy2 ) )
             minL=min(t%L, B(i)%T(ntier)%L) 
             if ( dist < minL ) then
                call CombineTiers(T, B(i)%T(ntier) )
             endif

             !do k=1,size(T%xy2(:,1))
             !   np=size(B(i)%T(j)%xy2(:,1))
             !   do p=1,np
             !      q=mod(p,np)+1
             !      dist = CNRLIN(B(i)%T(j)%xy2(p,1),B(i)%T(j)%xy2(p,1), //
             !                    B(i)%T(j)%xy2(q,1),B(i)%T(j)%xy2(q,1), //
             !                    T%xy2(k,1),        T%xy2(k,2) )   !original
             !      !dist=DISLIN2POINT(T%xy2(i,1), T%xy2(i,2), T%xy2(j,1),T%xy2(j,2), x_stack,y_stack)
             !  end do
             !  
             !enddo
          enddo
       endif
    endif
    enddo
end subroutine

!Combinar dos tiers 
subroutine CombineTiers(t1,t2) 
    implicit none
    type(tier),intent(inout) :: T1
    type(tier),intent(in   ) :: T2
    T1%xmin=min(T1%xmin, T2%xmin) 
    T1%xmax=max(T1%xmax, T2%xmax)
    T1%ymin=min(T1%ymin, T2%ymin)
    T1%ymax=max(T1%ymax, T2%ymax)

    !T1%wid=T1%xmax - T1%xmin
    !T1%len=T1%ymax - T1%ymin
    !!T1%hgt=T%h
    !T1%L  =min(T1%wid, T1%hgt)

end subroutine

!ADD TO TABLES: ***********************************************************************************
subroutine addToOutTable(table,S,T,d,wdir)
    implicit none
    type(outTable), intent(inout) :: table
    type(tier), intent(in) :: T
    type(stack),intent(in) :: S
    double precision, intent(in) :: wdir
    integer, intent(in) ::d
    table%stkName=S%nombre
    table%tabla(d,1:6)=[ sngl(wdir),t%hgt,t%wid,t%len,t%xbadj,t%ybadj ] !copy on table
end subroutine

subroutine addToStkTable(sTab,S,T,belev)
    implicit none
    type(stkTable), intent(inout) :: sTab
    type(tier), intent(in) :: T
    type(stack),intent(in) :: S
    real,intent(in) :: belev
    sTab%stkName     = S%nombre
    sTab%stkHeight   = S%h
    sTab%BaseElevDiff= S%z0-belev
    sTab%GEPEQN1     = T%gsh ! max(, gep)
    sTab%GEPSHV      = max(65.0, T%gsh)
end subroutine

subroutine no_tier_affects_this_stack_at_this_wdir(S,d,wdir,oTab,sTab)
    implicit none
    type(outTable), intent(inout) :: oTab
    type(stkTable), intent(inout) :: sTab
    type(stack),intent(in) :: S
    double precision, intent(in) :: wdir
    integer, intent(in) ::d
    oTab%stkName=S%nombre
    oTab%tabla(d,1:6)=[ sngl(wdir),0.0,0.0,0.0,0.0,0.0 ] !copy on table
    if ( .not. S%affected_by_tier) then
       print '("     No tier affects this stack. (",A10,")")',S%nombre
       sTab%stkName=S%nombre
       sTab%stkHeight=S%h
       sTab%BaseElevDiff=-99.99 !S%z0-belev
       sTab%GEPEQN1=0.0         !-99.99 !T%gsh ! max(, gep)
       sTab%GEPSHV=65.0         !max(65.0, T%gsh)
    end if
end subroutine

!INPUT:  ******************************************************************************************
subroutine readINP(inp_file,B,S,title) 
        implicit none
        character(78),intent(inout) :: title
        character(24),intent(in) :: inp_file
        TYPE(building),allocatable, intent(inout) :: B(:)
        TYPE(stack),allocatable, intent(inout) :: S(:)
        integer :: nb,nt,nn,ns !# builds,# tiers # nodes,# stacks
        integer :: i,j,k

        WRITE(*,'(/,A,/)') ' READING INPUT DATA FROM FILE.'
        open(1,file=inp_file,action="READ")
        
        read(1,*) title !titulo
        read(1,*) !options 'P'
        read(1,*) !units        & factor of correction
        read(1,*) !coord system & initial angle
        
        !BUILDINGS:
        read(1,*) nb    !# buildings
        allocate(B(nb)) 
        do i=1,nb,1 !read buildings and tiers
           read(1,*) B(i)%nombre,nt,B(i)%z0       !name ntiers z0
        
           allocate(B(i)%T(nt))
           !print*,"BUILDING: ",B(i)%nombre,B(i)%z0,nt              !debug
           do j=1,nt,1
              read(1,*) nn, B(i)%T(j)%h !hgt  !nn hgt
              allocate(B(i)%T(j)%xy(nn,2))
              allocate(B(i)%T(j)%xy2(nn,2))
              do k=1,nn,1
                 read(1,*) B(i)%T(j)%xy(k,1), B(i)%T(j)%xy(k,2)     !x y
              enddo
           enddo
        end do
        !STACKS:       
        read(1,*)ns  !#stacks
        allocate(S(ns)) 
        do i=1,ns,1 !read stacks
                read(1,*) S(i)%nombre, S(i)%z0, S(i)%h, S(i)%xy(1), S(i)%xy(2)    !name z0 h x y
                !print*,"STACK: ", S(i)%nombre, S(i)%z0,S(i)%h,S(i)%xy(1,:)          !debug
        end do

        close(1)
        WRITE(*,'(/,A,/)') ' END OF READING INPUT DATA FROM FILE.'

end subroutine

!OUTPUT: ******************************************************************************************
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
        WRITE (2,'(1X,A78,///)') TITLE
        !STACK RESULTS
        WRITE(2,"(16X,'PRELIMINARY* GEP STACK HEIGHT RESULTS TABLE')")
        WRITE(2,"(13X,'            (Output Units: meters)',/)")
        WRITE(2,"(8X,'                    Stack-Building            Preliminary*')")
        WRITE(2,"(8X,' Stack    Stack     Base Elevation    GEP**   GEP Stack')")
        WRITE(2,"(8X,' Name     Height    Differences       EQN1    Height Value',//)")
        do i=1,size(sT,1),1
           
           if ( st(i)%BaseElevDiff .EQ. -99.99 ) then
              !                                                   STKN(S),       SH(S),         GEP(S),    PV
              WRITE(2,'(8X, A8, F8.2, 10X, "N/A",5X,3(F8.2,5X))') st(i)%stkName, st(i)%stkHeight, st(i)%GEPEQN1, st(i)%GEPSHV 
           else
              !                                STKN(S),        SH(S),           DIF,                 GEP(S),       PV
              WRITE(2,'(8X, A8, 4(F8.2,5X))') st(i)%stkName, st(i)%stkHeight, st(i)%BaseElevDiff, st(i)%GEPEQN1, st(i)%GEPSHV
           end if
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
     WRITE (2,'(30X,"BPIP (Dated: 24241 )")')
     WRITE (2,'(1X, "DATE : ",I2,"/",I2,"/",I4)') IMON, IDAY, IYR
     WRITE (2,'(1X, "TIME : ",I2,":",I2,":",I2)') IHR, IMIN, ISEC
 end subroutine


!!! FUNCIONES HEREDADAS ********************************************************
!REAL FUNCTION CNRLIN (X1,Y1, X2,Y2, XKP,YKP)!, XI,YI)!, BET,DIST)
!   ! CNRLIN - SUBROUTINE TO CALCULATE THE DISTANCE BETWEEN A TIER CORNER AND
!   !           THE SIDE OF ANOTHER TIER.  SUBROUTINE ALSO CALCULATES WHETHER
!   !           OR NOT A PERPENDICULAR LINE DRAWN FROM THE CORNER TO THE SIDE
!   !           INTERCEPTS THE SIDE BETWEEN THE TWO CORNERS OF THE TIER.
!   !calculate corner perpendicular to side distance,
!   !intercept point, and determine if intercept on or between corners.
!   implicit none
!   double precision, intent(in) :: X1, Y1, X2,Y2, XKP,YKP,
!   double precision :: A1, A2, BET, DIST, SM, XI,YI
!
!   IF ((X1 .NE. X2) .AND. (Y1 .NE. Y2)) THEN
!      SM = (Y2 - Y1) / (X2 - X1)                                   !slope
!      XI = (YKP + XKP / SM - Y1 + X1 * SM) / (SM + 1.0 / SM)       !intercept-X
!      YI = Y1 + (XI - X1) * SM                                     !intercept-Y
!   ELSE
!     IF ((Y2 .EQ. Y1)) THEN
!        XI = XKP
!        YI = Y1
!     ELSE
!        XI = X1
!        YI = YKP
!     END IF
!   END IF
!   DIST = SQRT((YI - YKP) ** 2 + (XI - XKP) ** 2)
!
!!degug:
!print*,"DIST:",DIST,dislin2point(X1,Y1, X2,Y2, XKP,YKP)
!!
!    !
!    ! Is the intercept point between the two corners of the other structure ?
!    A1 = (X1 - XI) ** 2 + (Y1 - YI) ** 2 + (X2 - XI) ** 2 + (Y2 - YI) ** 2
!    A2 = (X1 - X2) ** 2 + (Y1 - Y2) ** 2
!    BET = (A2 - A1)
!
!    CNRLIN=sngl(DIST)
!END FUNCTION

end program
!
!LOGICAL FUNCTION DISLIN (X1, Y1, X2, Y2, XSP, YSP, L5)
!   !calculates distance between a side and stack
!   !and checks it against 5l
!   implicit none
!   double precision, intent(in) :: x1,x2,y1,y2,xsp,ysp
!   double precision, intent(in) :: L5
!   double precision             :: minX,maxX,yI,d1,d2,dist
!!   real :: tmp
!   DISLIN=.false.
!
!   minX = MIN (X1, X2)
!   maxX = MAX (X1, X2)
!
!   IF ( (XSP .LT. minX) .OR. (XSP .GT. maxX) )  RETURN
!
!   IF ( Y1 .EQ. Y2 ) THEN
!      DIST = YSP - Y1
!      IF ((DIST .GE. 0.0) .AND. (DIST .LE. L5)) THEN
!         DISLIN=.true.
!      END IF
!   END IF
!
!   IF (X1 .EQ. X2) THEN
!      IF (XSP .EQ. X1) THEN
!         D1 = YSP - Y1
!         IF ((D1 .LE. L5) .AND. (D1 .GE. 0.0)) THEN
!           DISLIN=.true.
!         END IF
!         D2 = YSP - Y2
!         IF ((D2 .LE. L5) .AND. (D2 .GE. 0.0)) THEN
!           DISLIN=.true.
!         END IF
!      END IF
!
!   ELSE
!      YI = Y2 + (XSP - X2) * (Y1 - Y2) / (X1 - X2)
!      DIST = YSP - YI
!      IF ((DIST .GE. 0.0) .AND. (DIST .LE. L5)) THEN
!         DISLIN=.true.
!      END IF
!   END IF
!!!DEBUG:
!!tmp=distLinePoint([x1,y1],[x2,y2],[xsp,ysp])
!!if (DIST-tmp > 0.5 ) print*,"DIST:",DIST,tmp
!!!ENDDEGUB
!END FUNCTION


!TRASH:
!real function minDist(xy1,xy2) result(dist)  !Calcula distancia mínima entre dos conjuntos de puntos
!    implicit none
!    double precision,intent(in) :: xy1(:,:), xy2(:,:)
!    double precision :: minvalue
!    integer :: i
!    dist=1.0e+12  
!    do i=1,size(xy1,1)
!       minvalue=minval( ( xy2(:,1)-xy1(i,1) )**2 + ( xy2(:,2) - xy1(i,2) )**2 ) 
!       dist=sngl(sqrt(min(Dist**2,minvalue )) )
!    end do
!end function
!
