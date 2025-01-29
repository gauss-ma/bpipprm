!  **************************************************************************
!  *                       (!) NON OFFICIAL CODE (!)                        *
!  *                                                                        *
!  *         BUILDING PROFILE INPUT PROGRAM for PRIME(DATED 25252)          *
!  *                                                                        *
!  *            *** SEE BPIPPRM MODEL CHANGE BULLETIN MCB#1 ***             *
!  *                                                                        *
!  *     ON THE SUPPORT CENTER FOR REGULATORY AIR MODELS BULLETIN BOARD     *
!  *                                                                        *
!  *                      http://www.epa.gov/scram                          *
!  *                                                                        *
!  **************************************************************************
!
!        Programmed by: Ramiro A. Espada  
!                       Lakes Environmental Software
!                       email: espada(at)agro.uba.ar
!
!  **************************************************************************
!
!        Written to: FORTRAN 90 Standards (Free Format)
!
!  **************************************************************************
 
      program bpipv25252
      IMPLICIT NONE
      !OBJETCS/TYPES ----------------------------------------------------------------
      type tier
        integer :: id
        double precision,allocatable :: xy(:,:)                         !x-y coordinates (from INP)
        real :: h,z0                                                    !height,base hgt (from INP)
        !projected (wdir-dependent) variables
        real, allocatable :: xy2(:,:)                                   !rotated coords          (changes on each wdir)
        real :: xmin,xmax,ymin,ymax                                     !boundaries              (changes on each wdir)
        real :: hgt=0,wid=0,len=0.0                                     !height,width,length     (changes on each wdir)
        real :: L=0.0                                                   !L: min{ wid , hgt }     (changes on each wdir)
        real :: gsh=0.0, xbadj=0.0, ybadj=0.0                           !GEP Stack Height Values (changes on each stack & wdir)
      endtype
      
      type building
        character(8)            :: bldName                              !name                (from INP)
        real                    :: z0                                   !base elevation (z0) (from INP)
        type(tier), allocatable :: t(:)                                 !tiers                    
      endtype
      
      type stack
        integer          :: id
        character(8)     :: stkName                                     !stack name
        real             :: z0, h                                       !base elev. (z0) and stack height
        double precision :: xy(2)                                       !stack coords (from INP)
        real             :: xy2(2)                                      !stack coords (rotated)

        logical          :: isOnRoof =.false.                           !boolean flag that states if stack is atop of a building
        integer          :: whichRoof=0                                 !id to the roof where this stack is placed
      endtype
      
      type outTable  !output table
        character(8) :: stkName
        real         :: tabla(36,6)=0.0                                 !36 windirs, 6vars: wdir,hgt,wid,len,xbadj,ybadj
      endtype
      
      type stkTable  !stacks table
        character(8) :: stkName
        real         :: stkHeight, BaseElevDiff=-99.99, GEPEQN1=0.0, GEPSHV=65.0
        !for summary report:
        real         :: gepbw,gepbh!,gepbl                              !GEP building wid, height and length
        integer      :: gtnum=0                                         !# of tiers affecting this stack
        real         :: gtdir=-99                                       !direcition at wich max gep ocurrs
        integer      :: tlist(10)=0                                     !list of tiers affecting stack
      endtype

      type smryData  !data for summary report
        type(tier)   :: sT(36),cT(36)               !max. GSH single & combined tier for each direction
        integer      :: Bid(36)=0,Tid(36)=0         !building and tier id that affects single tiers max gep on each wdir
        integer      :: gtlist(36,10)=0             !list of tiers contributing to maxGEP for combined tiers
        integer      :: gtnum(36)=0                 !# of combined-tiers affecting stack for each wdir.
      endtype

      !VARIABLES----------------------------------------------------------------------
      !global params:
      double precision, parameter :: pi=3.141593, pi2=6.28318530718_8 !3.141592653589793_8
      double precision, parameter :: deg2rad=pi/180.0, rad2deg=1.0/deg2rad
      character(24)               :: INPFILE="bpip.inp"
      character(24)               :: OUTFILE="bpip.out"
      character(24)               :: SUMFILE="bpip.sum"
      !setup values
      character (len=2)  :: SWTN = 'P '                                 !procesing for: "P":PRIME,"NP": no PRIME
                                                                        !"ST:"ISCST or "LT": ISCLT algorithms
      character (len=10) :: UNITS="METERS    ",UTMP  ="UTMN"            !UNITS !UTMP
      real :: FCONV = 1.00                                              !FCONV (from UNIT to meters)
      real :: PNORTH=0.00                                               !PNORTH (angle offset)
      !work variables:
      character(78)               :: TITLE                              !title of the run
      TYPE(building), allocatable :: B(:)                               !array of buildings
      TYPE(stack), allocatable    :: S(:)                               !array of stacks
      type(tier)                  :: fT,mT                              !"current" or "focal" Tier and "max. GSH" tier
      type(stack)                 :: Si                                 !"current" Stack
      double precision            :: wdir                               !wind direction
      logical                     :: SIZ,GSH                            !struc. influence zone boolean flag
      real                        :: maxGSH,refWID                      !max GSH and WID encountred (for given stack and wdir)
      integer                     :: mxtrs                              !max tier number found in input file
      real, allocatable           :: DISTMN(:,:),DISTMNS(:,:)           !distance Matrices: tier-tier, tier-stack 
      !vars used for tier combination
      type(tier)                  :: cT,T1,T2                           !"combined", "sub-group" and "merge-candidate" Tier
      integer,allocatable         :: TLIST(:,:)!TLIST2(:,:)             !list of combinable (w/focal) tiers indices
      integer,allocatable         :: TLIST2(:),mtlist(:)               !list of combined tiers ids
      integer                     :: TNUM,TNUM2,mtnum                  !TNUM= # of combinable tiers, TNUM2=# of actual combined tiers.
      real                        :: MNTDIST !min_tier_dist             !min distance between stack and combined tiers
      !out tables:
      type(outTable), allocatable  :: oTable(:)                          !output table
      type(smryData), allocatable  :: summary(:)                         !output table
      type(stkTable), allocatable  :: sTable(:)                          !stack  table
      !indices:
      integer :: i,j,k,d,i1,i2,nargs!,dd
      
      !get input, out and sum file paths from stdin
      nargs = command_argument_count()
      if(nargs>=1) CALL get_command_argument(1, INPFILE)
      if(nargs>=2) CALL get_command_argument(2, OUTFILE)
      if(nargs>=3) CALL get_command_argument(3, SUMFILE)
      
      !INPUT--------------------------------------------------------------------------
      call readINP(INPFILE,B,S,title,mxtrs)                       !read file & store data in B & S
      
      !Memory allocations and initializations:
      allocate(sTable(size(S)))                                       !allocate stack   table
      allocate(oTable(size(S)))                                       !allocate output  table
      allocate(summary(size(S)))                                      !allocate summary table
      
      allocate(DISTMNS(size(B)*mxtrs,size(S)))     ; DISTMNS=0.0 
      allocate(DISTMN(size(B)*mxtrs,size(B)*mxtrs));  DISTMN=0.0 
      allocate( TLIST(size(B)*mxtrs, 2 ))          ;   TLIST=0    
      allocate(TLIST2(size(B)*mxtrs))              ;  TLIST2=0    
      allocate(mtlist(size(B)*mxtrs))              ; mtlist=0    
      
      !MAIN---------------------------------------------------------------------------
      
      !Calculate things with "rotational invariance":
      call check_which_stack_over_roof(S,B)                             !check which stacks are placed over a roof.
      call calc_dist_stacks_tiers(S,B,DISTMNS)                          !calc min distance between stacks and tiers.
      call calc_dist_tiers(B,DISTMN)                                    !calc min distance between structures.
      
      !if ( mod(int(pnorth),360) /= 0 ) call rotate_coordinates(B,S,sngl(-pnorth*deg2rad)) 
      
      DO d=1,36                                                         !for each wdir (c/10 deg)
         print '(6X,"Wind flow passing", I4," degree direction.")',d*10 
         wdir=d*10.0*deg2rad - pnorth*deg2rad                           !get wdir [rad] 
         call rotate_coordinates(B,S,sngl(wdir))                        !rotate coordinates (T%xy S%xy -> T%xy2 S%xy2)
 

         DO i=1,size(S)                                                 !for each stack
           maxGSH=0.0; refWID=0.0                                       !reset maxGSH and refWID values

           Si=S(i)                                                      !init "current" stack

           !Single tiers structs   ---------------------------------------------------
           DO j=1,size(B)                                               !for each building
              DO k=1,size(B(j)%T)                                       !for each tier
      
                 fT=B(j)%T(k)                                           !create "focal" tier

                 !check if stack is inside SIZ
                 SIZ=           Si%xy2(1) .GE. (fT%xmin - 0.5*fT%L)  
                 SIZ=SIZ .AND. (Si%xy2(1) .LE. (fT%xmax + 0.5*fT%L) )
                 SIZ=SIZ .AND. (Si%xy2(2) .GE. (fT%ymin - 2.0*fT%L) )
                 SIZ=SIZ .AND. DISTMNS(fT%id,i) <= 5*fT%L
                 !SIZ=SIZ .AND. (Si%xy2(2) .LE. (fT%ymax + 5.0*fT%L) )      !(old)
                 if ( SIZ ) then 
                    fT%gsh   = fT%z0 + fT%hgt - Si%z0 + 1.5*fT%L        !GHS   [ Equation 1 (GEP, page 6) ]
                    fT%ybadj = Si%xy2(1) - (fT%xmin + fT%wid * 0.5)     !YBADJ = XPSTK - (XMIN(C) + TW*0.5)
                    fT%xbadj = fT%ymin - Si%xy2(2)                      !XBADJ = YMIN(C) - YPSTK             

                    GSH=ft%gsh > maxGSH
                    GSH=GSH.or.(fT%gsh == maxGSH .and. fT%wid < refWid)
                    if ( GSH ) then                                     !check if this tier has > GHS than previous ones
                       maxGSH=fT%gsh                                    !set new refGSH value
                       refWID=fT%wid                                    !store its WID
                       mT=fT                                            !set fT as the max GSH Tier
                       MTLIST=mt%id
                       MTNUM=1
                    end if

                    !!!@@(summary report): greatest GEP single tier.
                    if ( fT%gsh > summary(i)%sT(d)%gsh ) then
                       summary(i)%sT(d) = fT                            !save maxGSH tier on summary table
                       summary(i)%BId(d)= j                             !save maxGSH tier on summary table
                       summary(i)%TId(d)= k                             !save maxGSH tier on summary table
                    endif
                    !!!@@(summary report)
                 endif

              END DO!tiers
           END DO!buildings     

           !Combined tiers ----------------------------------------------------------
           DO j=1,size(B)                                               !for each building
              DO k=1,size(B(j)%T)                                       !for each tier
                fT=B(j)%T(k)                                           !create "focal" tier
                call list_combinable_tiers(fT, B, DISTMN,TLIST,TNUM)   !list combinable tiers and distances
                !look on combinable tiers (listed on TLIST) for "subgroups"
                if ( TNUM > 0 ) then
                   do i1=1,TNUM                                         !on each combinable tier
                      T1=B( TLIST(i1,1))%T( TLIST(i1,2) )               !create "subgrup" tier "T1"
           
                      if ( T1%hgt < fT%hgt .or. T1%id == fT%id ) then   !subgroup's common hgt should be < focal tier hgt
                         cT     = fT                                    !init. common ("combined") tier "cT"
                         cT%hgt = T1%hgt                                !use T1_hgt as "common subgroup height"
                         cT%L   = min(cT%hgt, cT%wid)                   !update L
           
                         TNUM2=0                                        !reset counter of combined tiers (tnum2) to 0
                         TLIST2=0                                       !reset list of combined tiers (tlist2) to 0
                         do i2=1,TNUM                                   !on each combinable candidate tier
                           T2=B(TLIST(i2,1))%T(TLIST(i2,2))             !create candidate tier "T2"
                           T2%L=min(cT%hgt, T2%wid)                     !use T2-L=min(T1hgt,T2wid) (line 1213 of original code)
                           if (T2%hgt>=cT%hgt .and. T2%id/= cT%id) then !combine condition 1) T2_hgt >= common hgt, and is not fT
                           if (DISTMN(cT%Id,T2%id)<max(T2%L,cT%L)) then !combine condition 2) dist T2-fT is less than the maxL
           
                              !"combine" tiers (cT,T2)                  !combine common Tier w/ T2 ! (update boundaries)
                              cT%xmin=min(cT%xmin, T2%xmin) 
                              cT%xmax=max(cT%xmax, T2%xmax)
                              cT%ymin=min(cT%ymin, T2%ymin)
                              cT%ymax=max(cT%ymax, T2%ymax)
           
                              TNUM2=TNUM2+1                             !increment counter of combined tiers (tnum2)
                              TLIST2(TNUM2)=T2%id                       !ad T2 id to list of combined tiers (TLIST2)
                           endif
                           endif  
                         enddo                                          !
                           !Once all tiers has been combined w/cT
                       if ( TNUM2 > 0 ) then                            !if there was at least 1 combined tier
           
                        cT%wid = cT%xmax - cT%xmin                      !update width
                        cT%len = cT%ymax - cT%ymin                      !update length
                        cT%L   = min(cT%hgt, cT%wid)                    !update "L"

                        !!
                        !!Here I need to define the Gap Filling Structure (GFS) polygon
                        !!and use it to determine if SIZ-L5 distance stack-GFS is satisfied
                        !!(could be solved using a "CONVEX HULL" algorithm, for example: Graham Scan Algorithm )
                        !!
                       MNTDIST=minval(DISTMNS([fT%id,TLIST2(:TNUM2)],i))!min dist from stack to all combined tiers

                        !check if stack is inside combined SIZ
                        SIZ=          Si%xy2(1) .GE. cT%xmin-0.5*cT%L   !this is how SIZ is defined for comb Tiers
                        SIZ=SIZ .AND. Si%xy2(1) .LE. cT%xmax+0.5*cT%L   !Struc Influence Zone (SIZ)
                        SIZ=SIZ .AND. Si%xy2(2) .GE. cT%ymin-2.0*cT%L   
                        SIZ=SIZ .AND. MNTDIST <= cT%L*5.0               !any of tiers combined is closer than cT%L*5.0 from stack
                        if (SIZ) then                                   !check if is in SIZ

                          cT%gsh  =cT%z0 + cT%hgt - Si%z0 + 1.5*cT%L     !GSH    [ Equation 1 (GEP, page 6) ]
                          cT%ybadj=Si%xy2(1) - (cT%xmin + cT%wid*0.5)    !YBADJ = XPSTK - (XMIN(C) + TW*0.5 )
                          cT%xbadj=cT%ymin - Si%xy2(2)                   !XBADJ = YMIN(C) - YPSTK

                          GSH=ct%gsh > maxGSH
                          GSH=GSH .or. ( cT%gsh == maxGSH .and. cT%wid < refWid )
                          if ( GSH ) then                                !check if comb tier has > GHS than previous ones
                             maxGSH=cT%gsh                               !set new ref GSH value
                             refWID=cT%wid                               !store its WID
                             mT=cT                                       !set cT as the max GSH Tier                   
                             MTLIST=tlist2                               !max gsh tiers affecting stack list
                             MTNUM =tnum2                                !num of tiers contributing to max gsh
                          end if                                         

                          !!!@@(summary report): greatest GEP combined tiers
                          if ( cT%gsh >= summary(i)%cT(d)%gsh ) then
                             summary(i)%cT(d)=cT                           !save maxGSH tier on summary table
                             summary(i)%gtlist(d,:)=[ft%id,tlist2(1:9)]    !
                             summary(i)%gtnum(d)=tnum2+1
                          endif
                          !!!@@(summary report)

                        endif !SIZ
                       endif!TNUM2>0
                      endif!T1%hgt < fT%hgt
                   enddo!T1
                endif!TNUM>0
              END DO!tiers
           END DO!buildings

           !Add results to Out-Table/Stk-Table
           oTable(i)%stkName  = Si%stkName
           sTable(i)%stkName  = Si%stkName
           sTable(i)%stkHeight= Si%h
           if ( maxGSH /= 0.0 ) then               !if ( any_tier_affects_stack_at_this_wdir ) then
              oTable(i)%tabla(d,1:6)=[ sngl(wdir),mT%hgt,mT%wid, mT%len,mT%xbadj,mT%ybadj ] 
           
              if ( maxGSH > sTable(i)%GEPEQN1 ) then
                 sTable(i)%BaseElevDiff = Si%z0 - mT%z0
                 sTable(i)%GEPEQN1      = mT%gsh 
                 sTable(i)%GEPSHV       = max(65.0, mT%gsh)
           
                 !!@@(summary report)
                 sTable(i)%gepbh = mT%hgt            !for sum report
                 sTable(i)%gepbw = mT%wid            !for sum report (!) este parece estar MAL
                 stable(i)%gtdir =sngl(wdir*rad2deg) !for sum report
                 sTable(i)%gtnum = mtnum             !for sum report
                 sTable(i)%tlist = mtlist            !for sum report
                 !!@@(summary report)
              endif

           else   !no tier affects this stack at this wdir
                 continue
           endif
           
           !if ( summary(i)%gtnum(d) < 2 ) print*,"no combined tiers affects stack",si%stkName,"(",si%id,") for this wdir."

         END DO!stacks
      END DO!wdir
      
      close(14)

      !OUTPUT:-----------------------------------------------------------------------
      call writeOUT(oTable,sTable,title,OUTFILE)
      call writeSUM(sumFile,S,B,summary,sTable,TITLE) !Stack summary

      print   '(/," END OF BPIP RUN." ,/)' 

      contains
      !GEOMETRY   ********************************************************************
      subroutine rotate_coordinates(B,S,wdir) 
      !Calc "new" (rotated) coordinates: xy --> xy2
          implicit none
          type(building),intent(inout) :: B(:)
          type(stack), intent(inout)   :: S(:)
          real,             intent(in) :: wdir  !original
          double precision :: D(2,2)                                   !usan distinta precision para stacks y buildings :/
          real             :: R(2,2)                                   !usan distinta precision para stacks y buildings :/
          integer :: i,j,k
          !matriz de rotación:       
          D(1,1)=dcos(dble(wdir)); D(1,2)=-dsin(dble(wdir)) !original
          D(2,1)=-D(1,2)   ; D(2,2)= D(1,1) 
          R(1,1)=cos(wdir) ; R(1,2)=-sin(wdir);
          R(2,1)=-R(1,2)   ; R(2,2)= R(1,1)
      
          !stacks
          do i=1,size(S)
             !S(i)%xy2=sngl(matmul(D,S(i)%xy))                                 !rotated stack coordinates
             S(i)%xy2=     matmul(R,sngl(S(i)%xy))                             !rotated stack coordinates
          enddo
          !buildings
          do i=1,size(B)
             do j=1,size(B(i)%T)
                do k=1,size(B(i)%T(j)%xy(:,1))
                  B(i)%T(j)%xy2(k,:)=sngl(matmul(D,B(i)%T(j)%xy(k,:))) !rotated tiers coordinates
                enddo
                call calc_tier_projected_values(B(i)%T(j) )            !calc tier: xmin,xmax,ymin,ymax,wid,len,hgt,L
             enddo
          enddo
      end subroutine
      
      real function dist_line_point (X0, Y0, X1, Y1, XP, YP)            !REPLACE of original "DISLIN" procedure
         !computes min distance between side/line defined by (v0,v1), and point (vp).
         !idea: use meaning of dot prod and cross product to get min dist from two vectors with same origin (v0)
         implicit none
         real, intent(in)    :: x0,x1,y0,y1,xp,yp
         real, dimension(2)  :: v,p,d                                   !side (v), point (p), and v-p vector (d)
         real                :: mod2v,v_dot_p,dist                      !length of v, dot_product(v,p)             
         v=[x1-x0, y1-y0] !v1-v0                                        ! take vectors to same origin
         p=[xp-x0, yp-y0] !vp-v0                                        ! take vectors to same origin
         mod2v = sqrt(dot_product(v,v))
         !proj_of_p_on_v=dot_product(v,p)/mod2v !proyeccion de p en v
         v_dot_p=dot_product(v,p) 
      
         if ( v_dot_p>mod2v**2 .or. v_dot_p<0.0 .or. mod2v==0.0 ) then
            d=v-p                                                       !diference vector
            dist=sqrt(dot_product(d,d))                                 !distance between the "tip" of the vectors
         else 
            dist=abs(v(1)*p(2)-v(2)*p(1))/mod2v                         !distance parametric line to point:  |v x p| / |v|
         end if
         dist_line_point = dist
      end function
      
      real function min_dist(xy1,xy2)            
         !min distance between two polygons, or poly (xy1) to point (xy2)
         real           ,intent(in) :: xy1(:,:), xy2(:,:)
         integer :: i,j,k,n,m
         m=size(xy1(:,1))
         n=size(xy2(:,1))
         min_dist=1e20
         do i=1,n
           do j=1,m
             k=mod(j,m)+1
             min_dist=min(min_dist,dist_line_point(xy1(j,1),xy1(j,2), xy1(k,1),xy1(k,2), xy2(i,1),xy2(i,2)) )  !new! (faster)
           enddo
         enddo
      end function
      
     logical function point_is_in_poly(point,poly)
         !idea: if point inside poly, then sum of angles from p to consecutive corners (sides) should be == 2*pi
         implicit none
         double precision, intent(in) :: point(2)
         double precision, intent(in) :: poly(:,:)
         double precision :: angle_sum=0.0,angle=0.0,signo=1.0
         double precision :: v1_dot_v1,v2_dot_v2
         double precision :: v1(2), v2(2)!, p(2)
         integer :: i,j,n!,k

         point_is_in_poly=.false.

         n=size(poly(:,1))
         do i=1,n
             j=mod(i,n)+1
             v1=poly(i,:) - point
             v2=poly(j,:) - point
             v1_dot_v1=dot_product(v1,v1)
             v2_dot_v2=dot_product(v2,v2) 
             if ( v1_dot_v1 == 0 .or. v2_dot_v2 == 0 ) then !this would means that v1 or v2 == point
              !point_is_in_poly=.false.! or .true.? (something to discuss)
              angle_sum=-999.9
              exit
             else
              signo=dsign(signo,v1(1)*v2(2)-v1(2)*v2(1))                 !sign of 3rd-component v1 x v2 (cross prod)
              angle=dacos(dot_product(v1,v2)/sqrt(v1_dot_v1*v2_dot_v2))
              angle_sum = angle_sum + signo * angle
             end if
         enddo
         point_is_in_poly=(abs( pi2 - abs(angle_sum) ) .lt. 1e-4) 
     end function point_is_in_poly


      !BPIP PARAM CALCULATIONS **************************************************************************     
      subroutine calc_tier_projected_values(T) 
          !calculate rotated/projected tier values: XMIN,XMAX,YMIN,YMAX,WID,HGT,LEN,L
          implicit none
          type(tier),intent(inout)   :: T
          T%xmin= minval(T%xy2(:,1)); T%xmax=maxval(T%xy2(:,1)) 
          T%ymin= minval(T%xy2(:,2)); T%ymax=maxval(T%xy2(:,2)) 
          T%wid = T%xmax - T%xmin
          T%len = T%ymax - T%ymin
          T%hgt = T%h
          T%L   = min(T%wid, T%hgt)
      end subroutine

      !STACK IS OVER ROOF? ******************************************************************************
      subroutine check_which_stack_over_roof(S,B)
          implicit none
          type(stack)   ,intent(inout) :: S(:)
          type(building),intent(in)    :: B(:)        
          integer                      :: i,j,k
          double precision              :: point(2)
          double precision, allocatable :: poly(:,:)
          double precision :: xmin,xmax,ymin,ymax
          print '("Detect if a stack is on top of a roof..")'

          do j=1,size(B)                                           !for each building
          !do k=1,size(B(j)%T)
          k=1
          if ( allocated(poly) )  deallocate(poly);
          allocate(poly(size(B(j)%T(k)%xy(:,1)),2))
          poly=B(j)%T(k)%xy(:,:)
          xmin=minval(poly(:,1)); xmax=maxval(poly(:,1)) 
          ymin=minval(poly(:,2)); ymax=maxval(poly(:,2)) 

          do i=1,size(S)                                           !for each stack
             if ( .not. S(i)%isOnRoof ) then !A stack should't be atop of >1 building at a time
                point=S(i)%xy(:)

                !Quick-check: (point whithin min-max of poly)
                if ( point(1) >= xmin .and. point(1) <= xmax )then
                if ( point(2) >= ymin .and. point(2) <= ymax )then

                   !if ( point_is_in_poly(point,poly) ) then
                   !   print '(A10,"=>",A10)',S(i)%stkName, B(j)%bldName
                      S(i)%isOnRoof  = .true.
                      S(i)%whichRoof = j !int((B(j)%T(k)%id+1) /mxtrs)
                   !else
                   !   continue
                   !end if
                end if
                end if

             end if
          end do
          !end do
          end do
      end subroutine
      
      !DISTANCES CALCULATIONS ***************************************************************************
      subroutine calc_dist_stacks_tiers(S,B,Matrix)
         implicit none
         type(stack),intent(in)       :: S(:)        
         type(building),intent(in)    :: B(:)        
         real,          intent(inout) :: Matrix(:,:)
         integer :: i,j,k,idt
         print '("Calculate min distance between tiers and stacks..")'
         do i=1,size(S)                                                 !on each stack
         do j=1,size(B)                                                 !on each building
         do k=1,size(B(j)%T)                                            !on each tier
            idt= B(j)%T(k)%id                                           !get tier id 
            if ( S(i)%whichRoof == B(j)%T(k)%id ) then
               Matrix(idt,i) = 0.0
            else
               Matrix(idt,i) = min_dist(sngl(B(j)%T(k)%xy), reshape(sngl(S(i)%xy),[1,2]))  !store min distance between structures
               !Matrix(idt,i) = min_dist(sngl(B(j)%T(k)%xy), reshape(S(i)%xy,[1,2]))  !store min distance between structures
            end if
         end do 
         end do 
         end do 
         !print '(25(F9.4))',Matrix !debug
      end subroutine
      
      subroutine calc_dist_tiers(B,Matrix)
         implicit none
         type(building),intent(in)    :: B(:)        
         real,          intent(inout) :: Matrix(:,:)
         integer :: i,j,n,ii,jj,id1,id2
         print '("Calculate min distance between buildings..")'
         n=size(B)
         do i=1,n                                                       !on each building
         do j=1,size(B(i)%T)                                            !on each tier
           id1= (i-1)*mxtrs + j 
           do ii=1,n                                                    !on each other building
           do jj=1,size(B(ii)%T)                                        !on each other tier
               id2= (ii-1)*mxtrs + jj
               if ( id1 > id2 ) then
                  Matrix(id1,id2) = min_dist(sngl(B(i)%T(j)%xy), sngl(B(ii)%T(jj)%xy)) !store min distance between structures
                  Matrix(id2,id1) = Matrix(id1,id2)                     !(symetry)
               end if
           end do 
           end do 
         end do 
         end do 
      end subroutine
      
      !MERGE TIERS **************************************************************************************
      subroutine list_combinable_tiers(T,B,DISTMN,TLIST,TNUM)
          implicit none
          type(tier)    ,intent(in)    :: T                             !"focal" Tier
          type(building),intent(in)    :: B(:)        
          real          ,intent(in)    :: DISTMN(:,:)
          integer       ,intent(inout) :: TLIST(:,:)                    !list of indices of combinable tiers
          integer       ,intent(inout) :: TNUM                          !# of combinable tiers
          real                         :: dist, maxL
          integer                      :: i,j
      
          tnum=0
          tlist=0
          do i=1,size(B)                                               !on each building
            do j=1,size(B(i)%T)                                        !on each tier    
              dist = DISTMN(T%id, B(i)%T(j)%id)                        !get dist
              maxL =    max(T%L , B(i)%T(j)%L )                        !"If the GREATER of each pair of Ls is greater than the minimum distance
              if ( dist < maxL ) then                                  !if dist less than maxL tiers are "COMBINABLE"
                 tnum=tnum+1
                 tLIST(tnum,:)=[i,j]
              endif
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
      
         print '(/," READING INPUT DATA FROM FILE.",/)' 
         OPEN(11, file=inp_file,action="READ",status='OLD')
           !HEADER: 
           read(11,*) TITLE                                              !titulo
           read(11,*) SWTN                                               !run options 'P', 'NP', 'ST', 'LT'
           read(11,*) UNITS, FCONV                                       !units        & factor of correction
           read(11,*) UTMP , PNORTH                                      !coord system & initial angle
      
           if ( SWTN(1:1) == 'S' .or. SWTN(1:1) =='L' ) then
           stop "This version of BPIP doesn't support 'SL'nor 'LS'"
           end if
           !BUILDINGS:
           read(11,*) nb                                                 !# buildings
           allocate(B(nb)) 
         
           do i=1,nb,1                                                  !read buildings and tiers
              read(11,*) B(i)%bldName,nt,B(i)%z0                          !name ntiers z0
              B(i)%z0=B(i)%z0*FCONV                                     !convert hgt units to meters
              mxtrs=max(mxtrs,nt)  
              allocate(B(i)%T(nt))
              do j=1,nt,1
                 B(i)%T(j)%z0=B(i)%z0                                   !asign tier same base elev tha building
                 read(11,*) nn, B(i)%T(j)%h !hgt                         !nn hgt
                 B(i)%T(j)%h=B(i)%T(j)%h*FCONV                          !convert hgt units to meters
                 allocate(B(i)%T(j)%xy(nn,2))
                 allocate(B(i)%T(j)%xy2(nn,2))
                 do k=1,nn,1
                    read(11,*) B(i)%T(j)%xy(k,1), B(i)%T(j)%xy(k,2) !x y
                 enddo
              enddo
           end do
           !STACKS:       
           read(11,*)ns                                                  !#stacks
           allocate(S(ns)) 
           do i=1,ns                                                    !read stacks
             read(11,*) S(i)%stkName,S(i)%z0,S(i)%h,S(i)%xy(1),S(i)%xy(2) !name z0 h x y
             S(i)%id=i                                                  !stacks are indexed by order of aparence on input file
             S(i)%z0=S(i)%z0*FCONV
           end do
         CLOSE(11)
         print '(/," END OF READING INPUT DATA FROM FILE.",/)'

         !----------------------------------------------------------------------
         ! INDEXING
         !Building Indexing:
         do i=1,size(B);do j=1,size(B(i)%T)                             !indexing tiers
           B(i)%T(j)%id=(i-1) * mxtrs + j                               !give each tier an absolute ID
         enddo; enddo;
         !Stacks Indexing: stacks are indexed while reading input file by the order of ocurrence
         !----------------------------------------------------------------------
      end subroutine
      !OUTPUT: ******************************************************************************************
      subroutine writeOUT(oT,sT,title,outputFileName)
          implicit none
          integer ::i,l
          character(78),  intent(in) :: title
          type(outTable), intent(in) :: ot(:)
          type(stkTable), intent(in) :: st(:)
          character(24),intent(in) :: outputFileName
      
          open(12,file=outputFileName,action="WRITE")
              WRITE(12,'(1X,A78,/)') TITLE
              !DATE:
              call writeDATE(12)
         WRITE(12,'(1X,A78,/)') TITLE
         WRITE(12,*) '============================'
         WRITE(12,*) 'BPIP PROCESSING INFORMATION:'
         WRITE(12,*) '============================'
         !
         WRITE(12,"(/3X,'The ',A2,' flag has been set for preparing downwash',' related data',10X)")  SWTN   !'P '
         WRITE(12,"(10X,'for a model run utilizing the PRIME algorithm.',/)"                       ) 
         WRITE(12,"(3X,'Inputs entered in ',A10,' will be converted to ','meters using ')"         )  UNITS  !"METERS    "
         WRITE(12,"(3X,' a conversion factor of',F10.4,'.  Output will be in meters.',/)"          )  FCONV  !1.00
         WRITE(12,"(3X,'UTMP is set to ',A4,'.  The input is assumed to be in',' a local')"        )  UTMP   !"UTMN"
         WRITE(12,"(3x,' X-Y coordinate system as opposed to a UTM',' coordinate system.')"        ) 
         WRITE(12,"(3x,' True North is in the positive Y',' direction.',/)")                                 ! 
         WRITE(12,"(3X,'Plant north is set to',F7.2,' degrees with respect to',' True North.  ')"  ) PNORTH ! 0.00
         WRITE(12,"(//1X,A78,///)") TITLE
         !STACK RESULTS
         WRITE(12,"(16X,'PRELIMINARY* GEP STACK HEIGHT RESULTS TABLE')")
         WRITE(12,"(13X,'            (Output Units: meters)',/)")
         WRITE(12,"(8X,'                    Stack-Building            Preliminary*')")
         WRITE(12,"(8X,' Stack    Stack     Base Elevation    GEP**   GEP Stack')")
         WRITE(12,"(8X,' Name     Height    Differences       EQN1    Height Value',//)")
         do i=1,size(sT,1)
            if ( st(i)%BaseElevDiff .EQ. -99.99 ) then
               !                                                   STKN(S),       SH(S),         GEP(S),    PV
               WRITE(12,'(8X, A8, F8.2, 10X, "N/A",5X,3(F8.2,5X))') st(i)%stkName, st(i)%stkHeight, st(i)%GEPEQN1, st(i)%GEPSHV 
            else
               !                                STKN(S),        SH(S),           DIF,                 GEP(S),       PV
               WRITE(12,'(8X, A8, 4(F8.2,5X))') st(i)%stkName,st(i)%stkHeight,st(i)%BaseElevDiff,st(i)%GEPEQN1, st(i)%GEPSHV
            end if
         end do
         WRITE(12,"(/,'   * Results are based on Determinants 1 & 2 on pages 1',' & 2 of the GEP')   ")
         WRITE(12,"(  '     Technical Support Document.  Determinant',' 3 may be investigated for')  ")
         WRITE(12,"(  '     additional stack height cred','it.  Final values result after')          ")
         WRITE(12,"(  '     Determinant 3 has been ta','ken into consideration.')                    ")
         WRITE(12,"(  '  ** Results were derived from Equation 1 on page 6 of GEP Tech','nical')     ")
         WRITE(12,"(  '     Support Document.  Values have been adjusted for a','ny stack-building') ")
         WRITE(12,"(  '     base elevation differences.',/)                                          ")
         WRITE(12,"(  '     Note:  Criteria for determining stack heights for modeling',' emission') ")
         WRITE(12,"(  '     limitations for a source can be found in Table 3.1 of the')              ")
         WRITE(12,"(  '     GEP Technical Support Document.')                                        ")
         WRITE(12,"(/,/,/,/)")
         !DATE (AGAIN)
         call writeDATE(12)
         WRITE(12,'(//,1X,A78,/)') TITLE
         WRITE(12, *) ' BPIP output is in meters'
      
         !MAIN OUTPUT:
         do i=1,size(oT,1),1
           write(12,'(/)')
           !HGT
           do l=1,36,6
              write(12,'(5X,"SO BUILDHGT ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(l:l+5,2)
           enddo
           !WID
           do l=1,36,6
              write(12,'(5X,"SO BUILDWID ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(l:l+5,3)
           enddo
           if ( SWTN(1:1) == 'P' .or. SWTN(1:1) == 'p' ) then
           !LEN
           do l=1,36,6
              write(12,'(5X,"SO BUILDLEN ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(l:l+5,4)
           enddo
           !XBADJ
           do l=1,36,6
              write(12,'(5X,"SO XBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(l:l+5,5)
           enddo
           !YBADJ
           do l=1,36,6
              write(12,'(5X,"SO YBADJ    ",a8,6(f8.2))') oT(i)%stKName,oT(i)%tabla(l:l+5,6)
           enddo
           endif
         end do
         close(12)!cierro bpip.out
      end subroutine
      
      subroutine writeDATE(iounit)
          implicit none
          integer, intent(in) :: iounit      
          integer :: date_time(8)
          integer :: iyr,imon,iday,ihr,imin,isec
          character(len=12) :: clock(3)
          CALL DATE_AND_TIME (CLOCK(1), CLOCK(2), CLOCK(3), DATE_TIME)
          IYR = DATE_TIME(1); IMON = DATE_TIME(2); IDAY = DATE_TIME(3)
          IHR = DATE_TIME(5); IMIN = DATE_TIME(6); ISEC = DATE_TIME(7)
          !header:
          WRITE (iounit,'(30X,"BPIP (Dated: 24244)")')
          WRITE (iounit,'(1X, "DATE : ",I2,"/",I2,"/",I4)')IMON,IDAY,IYR
          WRITE (iounit,'(1X, "TIME : ",I2,":",I2,":",I2)')IHR,IMIN,ISEC
      end subroutine
      
      !SUM File  ===================================================
      subroutine writeSUM(fileName,S,B,smT,sT,title)
         implicit none
         character(24),intent(in) :: FileName           !bpip.sum
         character(78) , intent(in) :: title            !run title
         type(stack)   , intent(in) :: S(:)             !all stacks
         type(building), intent(in) :: B(:)             !all buildings/tiers
         type(stkTable), intent(in) :: St(:)            !stacks table
         type(smryData), intent(in) :: smT(:)       !summary data
         type(tier)                 :: fT
         integer :: d,i,j,k
      
         print '("Writing summary report..")'

         open(14,file=FileName, action="WRITE")
            !-----------
            !Part 1: Report of Input Data.

            WRITE(14,'(1X,A78,/)') TITLE
            !DATE:
            call writeDATE(14)
            WRITE(14,'(1X,A78,/)') TITLE
            WRITE(14,*) "============================"
            WRITE(14,*) "BPIP PROCESSING INFORMATION:"
            WRITE(14,*) "============================"
            !Global options:
            WRITE(14,"(/3X,'The ',A2,' flag has been set for preparing downwash',' related data',10X)") SWTN
            WRITE(14,"(10X,'for a model run utilizing the PRIME algorithm.',/)"                       )
            WRITE(14,"(3X,'Inputs entered in ',A10,' will be converted to ','meters using ')"         ) UNITS
            WRITE(14,"(3X,' a conversion factor of',F10.4,'.  Output will be in meters.',/)"          ) FCONV
            WRITE(14,"(3X,'UTMP is set to ',A4,'.  The input is assumed to be in',' a local')"        ) UTMP
            WRITE(14,"(3x,' X-Y coordinate system as opposed to a UTM',' coordinate system.')"        )
            WRITE(14,"(3x,' True North is in the positive Y',' direction.',/)")
            WRITE(14,'(3X,"Plant north is set to",F7.2," degrees with respect to"," True North.  ",////)') PNORTH
            !
            WRITE(14,*) "=============="
            WRITE(14,*) "INPUT SUMMARY:"
            WRITE(14,*) "=============="
            WRITE(14, '(//,1X,"Number of buildings to be processed :",I4)') size(B) ! NB
            do i=1,size(B)
             WRITE(14,'(//1X,A8," has",I2," tier(s) with a base elevation of",F8.2," ",A10)') B(i)%bldName,size(B(i)%T),B(i)%z0,UNITS
               !TABLE:
               WRITE(14,'(" BUILDING  TIER  BLDG-TIER  TIER   NO. OF      CORNER   COORDINATES")')
               WRITE(14,'("   NAME   NUMBER   NUMBER  HEIGHT  CORNERS        X           Y"/)')
               do j=1,size(B(i)%T)
                 WRITE(14,'(1X,A8,I5,5X,I4,4X,F6.2,I6)') B(i)%bldName, j, B(i)%T(j)%id, B(i)%T(j)%h, size(B(i)%T(j)%xy(:,1))
                 do k=1, size(B(i)%T(j)%xy(:,1))
                    WRITE(14,'(42X,2F12.2, 1X,"meters")') B(i)%T(j)%xy(k,1),B(i)%T(j)%xy(k,2)
                    if (mod(int(pnorth),360) /= 0.0 ) then
                        WRITE(14,'(41X,"[",2F12.2,"] meters")') B(i)%T(j)%xy2(k,1),B(i)%T(j)%xy2(k,2)
                    endif
                 enddo
               enddo
            enddo
            !Stacks table:
            WRITE(14,'(/,1X,"Number of stacks to be processed :",I4,/)') size(S)
            WRITE(14, '("                    STACK            STACK   COORDINATES")')
            WRITE(14, '("  STACK NAME     BASE  HEIGHT          X           Y"/)')
            do i=1,size(S) 
               WRITE(14,'(2X, A8,3X, 2F8.2, 1X, A10)') S(i)%stkName, S(i)%z0, S(i)%h, UNITS
               WRITE(14,'(31X,2F12.2, " meters")') S(i)%xy(1),S(i)%xy(2)
               if (mod(int(pnorth),360) /= 0.0 ) then
                   WRITE(14,'(30X,"[",2F12.2,"] meters")') S(i)%xy2(1),S(i)%xy2(2)
               endif
            enddo
            !stack on roof table
            if (ANY( S(:)%isOnRoof) ) then
                WRITE(14,'(/,/," The following lists the stacks that have been identified",/," as being atop the noted building-tiers.",/)')
                WRITE(14,'("          STACK            BUILDING         TIER",/,9X," NAME      NO.    NAME        NO.  NO.")')
                do i=1,size(S) 
                  if ( s(i)%isOnRoof ) then
                     WRITE(14,'(10X, A8, I4, 5X, A8, 2(1X, I5))') S(i)%stkName, S(i)%id, B(S(i)%whichRoof)%bldName,S(i)%whichRoof, 1
                  end if
                enddo
            else
               WRITE(14,*) " "
               WRITE(14,*) "   No stacks have been detected as being atop"
            endif

           !-----------
           !Part 2: Overall MAX GEP stack summary (only needs stack table)
           WRITE(14,'(//,21X,"Overall GEP Summary Table",/,26X,"(Units: meters)")') 
           do i=1,size(St)
              WRITE(14,'(//," StkNo:", I3,"  Stk Name:", A8," Stk Ht:",F7.2," Prelim. GEP Stk.Ht:",F8.2)') i,st(i)%stkName,st(i)%stkHeight,st(i)%GEPSHV
              WRITE(14,'(12X,"GEP:  BH:",F7.2,"  PBW:",F8.2, 11X,"  *Eqn1 Ht:",F8.2)') st(i)%GEPBH,st(i)%GEPBW,st(i)%GEPEQN1                                                  !1022
              if ( st(i)%gtnum > 0) then       
                WRITE(14,'(10X,"*adjusted for a Stack-Building elevation difference"," of",F8.2)') st(i)%BaseElevDiff      !1025
                WRITE(14,'("  No. of Tiers affecting Stk:", I3,"  Direction occurred:", F8.2)') st(i)%gtnum,st(i)%gtdir    !1023
                WRITE(14,'("   Bldg-Tier nos. contributing to GEP:", 10I4)') [ (st(i)%tlist(j), j=1, min(10,st(i)%gtnum)) ]!1024
              else
                 WRITE(14,*) "     No tiers affect this stack."
              endif
           enddo
           !@-----------
           !@Part 3a: Single tier stack summary (need summary table (single tier only)+ stack table + buildings)
               WRITE(14,'(////,21X,"Summary By Direction Table",/,26X,"(Units:  meters)",/)')
               !-- SINGLE TIERS:      
               WRITE(14,'(" Dominate stand alone tiers:",/)')
               do d=1,36
               WRITE(14,'(/,1X,"Drtcn: ", F6.2/)') real(d*10)                                                         !604
               do i=1,size(S)
                  fT=smt(i)%sT(d) !single focal-max tier for this stack "i" and direction "d".
                  WRITE(14,'(" StkNo:", I3,"  Stk Name:", A8, 23X,"   Stack Ht:", F8.2)') S(i)%id,S(i)%stkName,S(i)%h    !2022
                  WRITE(14,'(17X,"GEP:  BH:",F7.2,"  PBW:",F7.2,"   *Equation 1 Ht:", F8.2)') st(i)%gepbh,st(i)%gepbw,st(i)%GEPEQN1  !2027

                  if ( smt(i)%BId(d) /= 0 ) then
                     WRITE(14,'(5X,"Single tier MAX:  BH:",F7.2,"  PBW:",F7.2,"  PBL:",F7.2,"  *Wake Effect Ht:", F8.2)')fT%hgt,fT%wid,fT%len,ft%gsh
                     WRITE(14,'(5X,"Relative Coordinates of Projected Width Mid-point: XADJ: ",F7.2,"  YADJ: ",F7.2/5X)')fT%xbadj,fT%ybadj 
                     WRITE(14,'(10X,"*adjusted for a Stack-Building elevation difference of",F8.2)') st(i)%BaseElevDiff  !S(i)%z0-fT%z0
                     WRITE(14,'(15X," BldNo:", I3,"  Bld Name:", A8, "  TierNo:", I3)') smt(i)%bId(d),B(smt(i)%bId(d))%bldName,smt(i)%tId(d)
                  else
                    WRITE(14,*) '    No single tier affects this stack for this direction.'!,st(i)%gtnum
                  endif
               enddo
               enddo
           !@-----------
           !@Part 3b: Combined tier stack summary (need summary table (combined tier only)+ stack table + buildings)
               !-- COMBINED TIERS:    
               WRITE(14,'(//," Dominant combined buildings:")')
               do d=1,36
               WRITE(14,'(/1X,"Drtcn: ", F6.2/)') real(d*10)                                                      !604
               do i=1,size(S)

               fT=smt(i)%cT(d) !combined focal-max tier for this stack "i" and direction "d".
               WRITE(14,'(" StkNo:", I3,"  Stk Name:", A8, 23X,"   Stack Ht:", F8.2)') S(i)%id,S(i)%stkName,S(i)%h !2022
               WRITE(14,'(17X,"GEP:  BH:",F7.2,"  PBW:",F7.2,"   *Equation 1 Ht:", F8.2)') st(i)%gepbh,st(i)%gepbw,st(i)%GEPEQN1  !2027
               
               if ( smt(i)%gtnum(d) /= 0) then !st(i)%gtnum /= -9 .and. 
                  WRITE(14,'(3X,"Combined tier MAX:  BH:",F7.2,"  PBW:",F7.2,"  PBL:",F7.2,"  *WE Ht:", F8.2)') ft%hgt,ft%wid,ft%len,ft%gsh
                  WRITE(14,'(5X,"Relative Coordinates of Projected Width Mid-point: XADJ: ",F7.2,"  YADJ: ",F7.2/5X)')ft%xbadj,ft%ybadj !2026
                  WRITE(14,'(10X,"*adjusted for a Stack-Building elevation difference of",F8.2)') st(i)%BaseElevDiff  !S(i)%z0-fT%z0
                  WRITE(14,'("  No. of Tiers affecting Stk:", I3)') smt(i)%gtnum(d)
                  WRITE(14,'("   Bldg-Tier nos. contributing to MAX:", 10I4)') smt(i)%gtlist(d,1:min(10,smt(i)%gtnum(d)))

               else
                  WRITE(14,*) '     No combined tiers affect this stack for this direction.'
               endif
               enddo
               enddo
            close(14)
      end subroutine

      end program
