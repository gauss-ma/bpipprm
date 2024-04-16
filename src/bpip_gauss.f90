program bpip_gauss
implicit none
!OBJETCS/TYPES ------------------------------------------------------------------------------------
type tier
  integer :: id
  double precision,allocatable :: xy(:,:)  !coords              (leido de inp NO CAMBIA)
  real :: h,z0                             !height,base hgt     (leido de inp NO CAMBIA)
  !projected:
  real,allocatable :: xy2(:,:)             !coords              (cambia para c/wdir)
  real :: xmin,xmax,ymin,ymax              !boundaries          (cambia para c/wdir)
  real :: hgt=0,wid=0,len=0.0              !height,width,length (cambia para c/wdir)
  real :: L=0.0                            !L: min{ wid , hgt } (cambia para c/wdir)
  real :: gsh=0.0, xbadj=0.0, ybadj=0.0    !gep values          (cambia para c/stack y wdir)
endtype

type building
  character(8)           :: nombre         !name
  real                   :: z0             !base height
  type(tier),allocatable :: t(:)           !tier
endtype

type stack
  character(8)     :: nombre
  real             :: z0, h
  real             :: xy(2), xy2(2)        !stack coords
  integer          :: whichRoof=0          !id to the roof where this stack is placed
endtype

type outTable  !output table
  character(8) :: stkName
  real         :: tabla(36,6)=0.0           !36 windirs, 6vars: wdir,hgt,wid,len,xbadj,ybadj
endtype

type stkTable  !stacks table
  character(8) :: stkName
  real         :: stkHeight, BaseElevDiff=-99.99, GEPEQN1=0.0, GEPSHV=65.0
endtype

!VARIABLES-----------------------------------------------------------------------------------------
integer :: i,j,k,d

!global params:
double precision, parameter :: pi=3.141593 !3.141592653589793_8
double precision, parameter :: deg2rad=pi/180.0
character(24),    parameter :: inputFileName="BPIP.INP"
character(24),    parameter :: outputFileName="bpip.out"

!work variables:
TYPE(building), allocatable :: B(:)  !array de edificios
TYPE(stack), allocatable    :: S(:)  !array de stacks
type(tier)                  :: fT,mT !"focal" tier, max. GSH Tier
character(78)    :: title
double precision :: wdir             !direccion del viento
logical          :: SIZ,ROOF
real             :: maxGSH,refWID    !max GSH and WID encountred (for given stack and wdir)
integer          :: mxtrs            !max tier number found in input file
! vars used for combined tiers
type(tier)          :: cT,T1,T2      ! "combined", "sub-group" and "merge-candidate" Tier
integer,allocatable :: TLIST(:,:)    !list of combinable (w/focal) tiers indices
real,allocatable    :: DISTMN(:,:)   !distances between all tiers
integer             :: l1,l2,TNUM         
!out tables:
type(outTable), allocatable :: oTable(:)  !output table
type(stkTable), allocatable :: sTable(:)  !stack  table

!INPUT--------------------------------------------------------------------------------------------
call readINP(inputFileName,B,S,title,mxtrs)   !read file & store data in B & S

allocate(oTable(size(S)))                                   !allocatar tabla de salida
allocate(sTable(size(S)))                                   !allocatar tabla de stacks
allocate(DISTMN(size(B)*mxtrs,size(B)*mxtrs)); DISTMN=0.0 
allocate( TLIST(size(B)*mxtrs,2))            ; TLIST=0    
!MAIN----------------------------------------------------------------------------------------------

do i=1,size(B);do j=1,size(B(i)%T)                                     !indexing tiers
  B(i)%T(j)%id=(i-1) * mxtrs + j                                       !give each tier an absolute ID
enddo; enddo;

call check_which_stack_over_roof(S,B)                                  !check which stacks are placed over a roof.
call calc_mindist_between_tiers(B,DISTMN)                              !calc min distance between structures.

DO d=1,36                                                              !for each wdir (c/10 deg)
   print '("      Wind flow passing", I4," degree direction.")', d*10 
   wdir=d*10.0*deg2rad                                                 !get wdir [rad]

   call rotateCoords(B,S,sngl(wdir))                                   !rotate coordinates

   DO i=1,size(S)                                                      !for each stack
     maxGSH=0.0; refWID=0.0

     !Single tiers structs   ----------------------------------------------------------------------
     DO j=1,size(B)                                                        !for each building
        DO k=1,size(B(j)%T)                                                !for each tier

           fT=B(j)%T(k)                                                    !create "focal" tier

           SIZ  = isInsideSIZ(S(i), fT)                                    !check if stack on SIZ
           ROOF = S(i)%whichRoof == fT%id                                  !check if stack over ROOF
           if ( SIZ .or. ROOF ) then                                       !check if any condition meet

              call calcGepValues(fT,S(i))                                  !calc: GSH, XBADJ, YBADJ

              if (fT%gsh>maxGSH .or. (fT%gsh==maxGSH .and. refWID>=fT%wid)) then !check if this tier has
                                                                                 ! > GHS than previous ones
                 maxGSH=fT%gsh                                             !set new ref GSH value
                 refWID=fT%wid                                             !store its WID
                 mT=fT                                                     !set max GSH Tier as this fT
              end if
           endif
        END DO!tiers
     END DO!buildings

     !"COMBINED TIERS" ----------------------------------------------------------------------------
     DO j=1,size(B)                                                           !
        DO k=1,size(B(j)%T)                                                   !for each focal tier

           fT=B(j)%T(k)                                                       !create "focal" tier

           call getListOfCombinableTiers(fT,B,DISTMN,TLIST,TNUM) !list combinable tiers and distances

           !if (.false.   ) then
           if ( TNUM > 0 ) then
             do l1=1,TNUM                                                     !on each combinable tier
                T1=B(TLIST(l1,1))%T(TLIST(l1,2))                              !create "subgrup" tier "T1"

                if ( T1%hgt < fT%hgt .or. T1%id == fT%id  ) then              !common hgt should be < focal tier hgt
                   cT=fT                                                      !create common ("combined") tier "cT"
                   cT%hgt = T1%hgt                                            !use T1_hgt as "common subgroup hgt"
                   cT%L   = min(cT%hgt, cT%wid)                               !update L
!debug if(d==11)   print*,"candidate",B(j)%nombre, B(TLIST(l1,1))%nombre,distmn(T1%id,fT%id),ct%L ,T1%L

                   do l2=1,TNUM                                               !on each combinable candidate tier
                      T2=B(TLIST(l2,1))%T(TLIST(l2,2))                        !create candidate tier "T2"
                      if ( T2%hgt >= cT%hgt .and. T2%id /= cT%id ) then       !if T2_hgt >= common hgt, and is not cT

                         if ( DISTMN(fT%Id,T2%id) < max(T2%L, cT%L)) then     !if dist T2-fT is less than the maxL
                            call combineTiers(cT,T2)                          !combine common Tier w/ T2 ! (update boundaries)
!debug if(d==11)   print*,"combine! ",B(j)%nombre, B(TLIST(l2,1))%nombre,distmn(T2%id,fT%id),ct%L,T2%L
                        endif
                      endif  
                   enddo                                                       !
                   !hasta acá va OK.                                           !Once all tiers has been combined w/tC
                   cT%wid = cT%xmax - cT%xmin                                  !update wid
                   cT%len = cT%ymax - cT%ymin                                  !update len
                   cT%L   = min(cT%hgt, cT%wid)                                !update L

                   !cT%xy2(1,1)  = cT%xmin; cT%xy2(1,2)  = cT%ymin             !define Gap Filling Structure (GFS)
                   !cT%xy2(2,1)  = cT%xmin; cT%xy2(2,2)  = cT%ymax             !over xy2 coordinates
                   !cT%xy2(3,1)  = cT%xmax; cT%xy2(3,2)  = cT%ymax
                   !cT%xy2(4,1)  = cT%xmax; cT%xy2(4,2)  = cT%ymin
                   !cT%xy2(5:,1) = cT%xmin; cT%xy2(5:,2) = cT%ymin
                                                                               !and do the same as w/single tiers
                   SIZ  = isInsideSIZ(S(i), cT)
                   ROOF = S(i)%whichRoof == cT%id  ! (!) esto creo que hay que cambiarlo!
                   if ( SIZ .OR. ROOF ) then 
                                                                                       
                      call calcGepValues(cT,S(i))     !calc: GSH, XBADJ, YBADJ
                                                                                       
                      if ( maxGSH < cT%gsh .OR. ( maxGSH == cT%gsh .AND. refWID >= cT%wid) ) then 
                         maxGSH=cT%gsh
                         refWID=cT%wid
                         mT=cT
                      end if
                   endif
                 endif
             enddo
           endif
        END DO!tiers
     END DO!buildings
     
     !Add results to Out-Table/Stk-Table
     oTable(i)%stkName  = S(i)%nombre
     sTable(i)%stkName  = S(i)%nombre
     sTable(i)%stkHeight= S(i)%h
     if ( maxGSH /= 0.0 ) then    !if ( any_tier_affects_stack_at_this_wdir ) then
        oTable(i)%tabla(d,1:6)=[ sngl(wdir),mT%hgt,mT%wid,mT%len,mT%xbadj,mT%ybadj ] 
        if ( maxGSH > sTable(i)%GEPEQN1 ) then
           sTable(i)%BaseElevDiff = S(i)%z0 - mT%z0
           sTable(i)%GEPEQN1      = mT%gsh 
           sTable(i)%GEPSHV       = max(65.0, mT%gsh)
        endif
     !else   !no tier affects this stack at this wdir
     !   print '("           No tier affects this stack at this wdir. (",A10,")")',S(i)%nombre
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
    double precision :: D(2,2)    !usan distinta precision para stacks y buildings :/
    real             :: R(2,2)   !usan distinta precision para stacks y buildings :/
    integer :: i,j,k
    !matriz de rotación:       
    D(1,1)=dcos(dble(wdir)); D(1,2)=-dsin(dble(wdir)) !original
    D(2,1)=-D(1,2)   ; D(2,2)= D(1,1) 
    R(1,1)=cos(wdir) ; R(1,2)=-sin(wdir);
    R(2,1)=-R(1,2)   ; R(2,2)= R(1,1)

    !stack
    do i=1,size(S)
       S(i)%xy2=matmul(R,S(i)%xy)                                   !projected stack coordinates
    enddo
    !buildings
    do i=1,size(B)
       do j=1,size(B(i)%T)
           do k=1,size(B(i)%T(j)%xy(:,1))
              B(i)%T(j)%xy2(k,:)=sngl(matmul(D,B(i)%T(j)%xy(k,:)) ) !projected tiers coordinates
           enddo
           call calcTierProyectedValues(B(i)%T(j) )                 ! calc tier: xmin,xmax,ymin,ymax
       enddo
    enddo
end subroutine

real function DisLin2Point (X0, Y0, X1, Y1, XP, YP)    !REPLACE of original "DISLIN" procedure
   !computes min distance between side/line defined by (v0,v1), and point (vp).
   implicit none
   real, intent(in)    :: x0,x1,y0,y1,xp,yp
   real, dimension(2)  :: v,p,d                !side (v), point (p), and v-p vector (d)
   real                :: mod2v,v_dot_p,dist   !length of v, dot_product(v,p)             

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
end function

real function minDist(xy1,xy2)    !min distance between two polygons
     real           ,intent(in) :: xy1(:,:), xy2(:,:)
     integer :: i,j,k,n,m
     m=size(xy1(:,1))
     n=size(xy2(:,1))
     minDist=1e20
     do i=1,n
     do j=1,m
        k=mod(j,m)+1
        minDist=min(minDist,dislin2point( xy1(j,1),xy1(j,2),  xy1(k,1),xy1(k,2), xy2(i,1), xy2(i,2)) )  !new! (faster)
     enddo
     enddo
end function

!STRUCTURE INFLUENCE ZONE (SIZ) *******************************************************************
logical function isInsideSIZ(S,T)       result(SIZ)  
    ! check if Stack is over influence of a tier (Structure Influence Zone, SIZ)
    implicit none
    type(stack), intent(in) :: S
    type(tier), intent(in)  :: T
    real            :: x_stack,y_stack
    integer ::i,j,n
    x_stack=S%xy2(1) 
    y_stack=S%xy2(2) 
    
    SIZ=           x_stack .GE. (T%xmin - 0.5*T%L)  
    SIZ=SIZ .AND. (x_stack .LE. (T%xmax + 0.5*T%L) )
    SIZ=SIZ .AND. (y_stack .GE. (T%ymin - 2.0*T%L) )
    !SIZ=SIZ .AND. (y_stack .LE. (T%ymax + 5.0*T%L) )  !old
    if ( SIZ ) then
       n=size(T%xy2(:,1))
       do i=1,n
           j=mod(i,n)+1
           SIZ = dislin2point(T%xy2(i,1), T%xy2(i,2), T%xy2(j,1),T%xy2(j,2), x_stack,y_stack) <= T%L*5.0
           if (SIZ) return!exit
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
                S(i)%whichRoof = (I-1) * MXTRS + J ![j,k] 
                return!exit
             endif
          enddo
       enddo
    enddo
end subroutine

!BPIP PARAM CALCULATIONS **************************************************************************
subroutine calcTierProyectedValues(T) !Calculo de XMIN XMAX YMIN YMAX,WID,HGT,LEN,L
    implicit none
    type(tier),intent(inout)   :: T
    T%xmin= minval(T%xy2(:,1)) ; T%xmax=maxval(T%xy2(:,1)) 
    T%ymin= minval(T%xy2(:,2)) ; T%ymax=maxval(T%xy2(:,2)) 
    T%wid = T%xmax - T%xmin
    T%len = T%ymax - T%ymin
    T%hgt = T%h
    T%L   = min(T%wid, T%hgt)
end subroutine

subroutine calcGepValues(T,S) !Calculo de GEP Stack Height, XBADJ, YBADJ
    implicit none
    type(tier),intent(inout)   :: T
    type(stack), intent(inout) :: S
    real                       :: x_stack, y_stack
    x_stack = S%xy2(1)
    y_stack = S%xy2(2)
    T%gsh   = T%z0 + T%hgt - S%z0 + 1.5*T%L     !Equation 1 (GEP, page 6)
    T%ybadj = x_stack - (T%xmin + T%wid * 0.5)  !YBADJ = XPSTK - (XMIN(C) + TW * 0.5)
    T%xbadj = T%ymin - y_stack                  !XBADJ = YMIN(C) - YPSTK             
end subroutine

!MERGE TIERS **************************************************************************************
subroutine calc_mindist_between_tiers(B,Matrix)
   implicit none
   type(building),intent(in)    :: B(:)        
   real,          intent(inout) :: Matrix(:,:)
   integer :: i,j,n,ii,jj,id1,id2
   
   print*, "Calculate min distance between buildings"
   n=size(B)
   do i=1,n                         !on each building
   do j=1,size(B(i)%T)              !on each tier
      id1= (i-1)*mxtrs + j 
      do ii=1,n                      !on each other building
      do jj=1,size(B(ii)%T)          !on each other tier
          id2= (ii-1)*mxtrs + jj
          if ( id1 > id2 ) then
             Matrix(id1,id2) = mindist(sngl(B(i)%T(j)%xy), sngl(B(ii)%T(jj)%xy)) !store min distance between structures
             Matrix(id2,id1) = Matrix(id1,id2)                                   !(symetry)
          end if
        end do 
        end do 
   end do 
   end do 
   print '(6(F8.3))',Matrix !debug
end subroutine

subroutine combineTiers(t1,t2) 
    implicit none
    type(tier),intent(inout) :: T1
    type(tier),intent(in   ) :: T2
    T1%xmin=min(T1%xmin, T2%xmin) 
    T1%xmax=max(T1%xmax, T2%xmax)
    T1%ymin=min(T1%ymin, T2%ymin)
    T1%ymax=max(T1%ymax, T2%ymax)
end subroutine

subroutine getListOfCombinableTiers(T,B,DISTMN,list,nC)
    implicit none
    type(tier)    ,intent(in)    :: T             !"focal" Tier
    type(building),intent(in)    :: B(:)        
    real          ,intent(in)    :: DISTMN(:,:)
    integer       ,intent(inout) :: list(:,:)       
    integer       ,intent(inout) :: nc              
    real                         :: dist, maxL
    integer                      :: i,j
     nc=0
     do i=1,size(B)                                  !on each building
          do j=1,size(B(i)%T)                        !on each tier    
             
             dist = DISTMN(T%id, B(i)%T(j)%id)-1e-2        !get dist
             maxL = max(T%L, B(i)%T(j)%L)            !"If the GREATER of each pair of Ls is greater than the minimum distance
                                                     ! between the two tiers, then the two tiers are considered to be combinable."
             if ( dist < maxL ) then                 !if dist less than maxL tiers are "COMBINABLE"
                nc=nc+1
                list(nc,1:2) = [i,j]
             endif
             if (T%id ==1 .and. i==2) print*,"c1,c2:",dist,maxL,dist<maxL
          enddo
     enddo
end subroutine

!INPUT:  ******************************************************************************************
subroutine readINP(inp_file,B,S,title,mxtrs) 
        implicit none
        character(78),intent(inout) :: title
        integer      ,intent(inout) :: mxtrs
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
        read(1,*) nb                              !# buildings
        allocate(B(nb)) 
        do i=1,nb,1                               !read buildings and tiers
           read(1,*) B(i)%nombre,nt,B(i)%z0       !name ntiers z0
           mxtrs=max(mxtrs,nt)  
           allocate(B(i)%T(nt))
           !print*,"BUILDING: ",B(i)%nombre,B(i)%z0,nt              !debug
           do j=1,nt,1
              B(i)%T(j)%z0=B(i)%z0  
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
                !print*,"STACK: ", S(i)%nombre, S(i)%z0,S(i)%h,S(i)%xy(1,:)       !debug
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
        WRITE(2,"('          for a model run utilizing the PRIME algorithm.',/)"                 )
        WRITE(2,"(3X,'Inputs entered in ',A10,' will be converted to ','meters using ')"         ) "METERS    "
        WRITE(2,"(3X,' a conversion factor of',F10.4,'.  Output will be in meters.',/)"          ) 1.00
        WRITE(2,"(3X,'UTMP is set to ',A4,'.  The input is assumed to be in',' a local')"        ) "UTMN"
        WRITE(2,"(3x,' X-Y coordinate system as opposed to a UTM',' coordinate system.')"        )
        WRITE(2,"(3x,' True North is in the positive Y',' direction.',/)")
        WRITE(2,"(3X,'Plant north is set to',F7.2,' degrees with respect to',' True North.  ',//)") 0.00
        WRITE (2,'(1X,A78,///)') TITLE
        !STACK RESULTS
        WRITE(2,"(16X,'PRELIMINARY* GEP STACK HEIGHT RESULTS TABLE')")
        WRITE(2,"(13X,'            (Output Units: meters)',/)")
        WRITE(2,"(8X,'                    Stack-Building            Preliminary*')")
        WRITE(2,"(8X,' Stack    Stack     Base Elevation    GEP**   GEP Stack')")
        WRITE(2,"(8X,' Name     Height    Differences       EQN1    Height Value',//)")
        do i=1,size(sT,1)
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

!subroutine writeSUM()
!end subroutine
end program
