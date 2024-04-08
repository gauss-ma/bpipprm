program bpip_gauss
!MODULES-------------------------------------------------------------
use struc            !objetos
use readBPIP         !lectura de datos de entrada
use writeBPIP        !escritura de salidas

!VARIABLES-----------------------------------------------------------
implicit none
integer :: i,j,k,d           !indices(para: stacks,buildings,tiers,wdirs)

!global
double precision, parameter :: pi=3.141593 !3.141592653589793_8
double precision, parameter :: deg2rad=pi/180.0
character(24),    parameter :: inputFileName="BPIP.INP"
character(24),    parameter :: outputFileName="bpip.out"

!objects
TYPE(build), allocatable :: B(:) !array de edificios
TYPE(stack), allocatable :: S(:) !array de stacks

double precision :: wdir         !direccion del viento

logical :: L5,SIZ,inTier,tier_affects_stack
!output vars:
character(78) :: title
real    :: refGSH,refWID
integer :: bmax,tmax         !index of max GSH for building and tier respectivamente

!out tables:
type(outTable), allocatable :: oTable(:)  !output table
type(stkTable), allocatable :: sTable(:)  !stack  table

!INPUT---------------------------------------------------------------
call readINP(inputFileName,B,S,title)   !read file y guardar datos en B y S

allocate(oTable(size(S)))   !allocatar tabla de salida
allocate(sTable(size(S)))   !allocatar tabla de stacks
!MAIN----------------------------------------------------------------

!print*, "Detect if a stack is on top of a roof"
!DO i=1,size(S)                                  !for each stack
!   DO j=1,size(B)                                !for each building
!   print '("(",A10,"-",A10,")")',S(i)%nombre,B(j)%nombre
!
!       DO k=1,size(B(j)%T)                       !for each tier
!         if ( isInsideTier(S(i), B(j)%T(k)) ) then
!            print '(A10,"=>",A10)',S(i)%nombre, B(j)%nombre
!         endif
!      enddo
!   enddo
!enddo

DO d=1,36   !for each wdir (c/10 grados)
     print '("      Wind flow passing", I4," degree direction.")', d*10 
     wdir=d*10.0*deg2rad    !wdir [rad]

     call rotar(B,S,wdir)  !roto coordenadas según wdir

     DO i=1,size(S)                                    !for each stack
         tier_affects_stack=.false.
         S(i)%affected_by_tier=S(i)%affected_by_tier .OR. .false.
         refGSH=0.0; refWID=0.0; tmax=0; bmax=0

         DO j=1,size(B)                                !for each building
             DO k=1,size(B(j)%T)                       !for each tier
          
                !check si tier afecta stack (adentro de 5L y SIZ)
                SIZ=isInsideSIZ(S(i), B(j)%T(k)) 
                L5 =isInsideL5( S(i), B(j)%T(k))
                inTier=( SIZ .and. S(i)%xy2(2) .LE. B(j)%T(k)%ymax ) 
                if ( ( L5 .AND. SIZ) .OR. inTier ) then 
                    S(i)%affected_by_tier=.true.
                    tier_affects_stack=.true.
                    call mergeCloseTiers(B(j)%T(k),B,k,B(j)%nombre)   !si hay otra structura cercana (de == tier) y combinarlos
                    call calcGepStackHeight(B(j)%z0, B(j)%T(k), S(i)) !calc: GSH, XBADJ, YBADJ
                   !guardar el tier de máximo GSH, (si hay dos con == GSH) quedarme con el de narrower width
                   if ( refGSH < B(j)%T(k)%gsh .OR. ( refGSH == B(j)%T(k)%gsh .AND. refWID >= B(j)%T(k)%wid) ) then
                        refGSH=B(j)%T(k)%gsh   
                        refWID=B(j)%T(k)%wid
                        bmax=j      !building max (index)
                        tmax=k      !tier     max (index)
                   end if
                endif
             END DO!tiers
         END DO!buildings
         
         if (tier_affects_stack) then
            call addToOutTable(oTable(i),S(i),B(bmax)%T(tmax),d,wdir) !Agregar a Tablas
            if ( refGSH > sTable(i)%GEPEQN1 ) then
               call addToStkTable(sTable(i),S(i),B(bmax)%T(tmax), B(bmax)%z0)
            endif
         else
            call no_tiers_affect_this_stack(S(i),d,wdir, oTable(i),sTable(i))
         endif
     END DO!stacks
END DO!wdir

!OUTPUT:-------------------------------------------------------------
call writeOUT(oTable,sTable,title,outputFileName)

contains


subroutine setTierProyectedValues(T) !Calculo de XMIN XMAX YMIN YMAX
        implicit none
        type(tier),intent(inout)   :: T
        T%xmin=sngl(minval(T%xy2(:,1)) ); T%xmax=sngl(maxval(T%xy2(:,1)) )
        T%ymin=sngl(minval(T%xy2(:,2)) ); T%ymax=sngl(maxval(T%xy2(:,2)) )

        T%wid=T%xmax - T%xmin
        T%len=T%ymax - T%ymin
        T%hgt=T%h
        T%L  =min(T%wid, T%hgt)
end subroutine

subroutine rotar(B,S,wdir) !Calculo de nuevas coordenadas
        implicit none
        type(build),intent(inout)   :: B(:)
        type(stack), intent(inout)  :: S(:)
        double precision,intent(in) :: wdir
        double precision :: R(2,2)
        integer :: i,j,k
        !matriz de rotación:       
        R(1,1)=dcos(wdir); R(1,2)=-dsin(wdir)
        R(2,1)=-R(1,2)   ; R(2,2)= R(1,1) 
        !stack
        do i=1,size(S)
           S(i)%xy2=matmul(R,S(i)%xy)                    !proyected coords stack
        enddo
        !buildings
        do i=1,size(B)
           do j=1,size(B(i)%T)
               do k=1,size(B(i)%T(j)%xy(:,1))
                  B(i)%T(j)%xy2(k,:)=matmul(R,B(i)%T(j)%xy(k,:))   !proyected coords tier corners
               enddo
               call setTierProyectedValues(B(i)%T(j) )       ! calc tier: xmin,xmax,ymin,ymax
           enddo
        enddo
end subroutine

!Calcula distancia mínima entre dos conjuntos de puntos
function minDist(xy1,xy2) result(Dist)
        implicit none
        double precision,intent(in) :: xy1(:,:), xy2(:,:)
        double precision :: Dist
        integer :: i
        Dist=1.0e+12  
        do i=1,size(xy1,1)
           Dist=sqrt(min(Dist**2, minval( ( xy2(:,1) - xy1(i,1) )**2 + ( xy2(:,2) - xy1(i,2) )**2 ) ) )
        end do
end function

function isInsideL5(S,T)   result(L5)
        implicit none
        logical :: L5
        type(stack), intent(in) :: S
        type(tier), intent(in)  :: T
        double precision        :: dist
        !double precision:: x_stack,y_stack
        !  Distancia Stack Tier debe ser menor a 5L
        dist=sqrt( minval( ( T%xy2(:,1) - S%xy2(1) )**2 + ( T%xy2(:,2) - S%xy2(2) )**2 ) ) 
        L5=( dist .LE. 5.0*T%L )
        !x_stack=S%xy2(1)!1,
        !y_stack=S%xy2(2)!1,
        !L5=          y_stack .LE. (T%ymax + 5.0*T%L) 
        !L5=L5 .AND. (y_stack .GE. T%ymin )
        !L5=L5 .AND. (x_stack .LE. T%xmax ) .AND. ( x_stack .GE. T%xmin ) 
end function

function isInsideSIZ(S,T)       result(SIZ)
        implicit none
        type(stack), intent(in) :: S
        type(tier), intent(in)  :: T
        logical :: SIZ
        double precision:: x_stack,y_stack
        x_stack=S%xy2(1)
        y_stack=S%xy2(2)
        !  Stack está dentro del SIZ de la estructura?  SIZ= [xmin-0.5*L , xmax + 0.5*L] , [ymin -2L , ymax+5L ]
        SIZ=           x_stack .GE. (T%xmin - 0.5*T%L)    
        SIZ=SIZ .AND. (x_stack .LE. (T%xmax + 0.5*T%L) )  
        SIZ=SIZ .AND. (y_stack .GE. (T%ymin - 2.0*T%L) )  
end function

function isInsideTier(S,T)    result(inTier) 
        implicit none
        type(stack), intent(in) :: S
        type(tier), intent(in)  :: T
        logical :: inTier
        double precision:: angle
        double precision:: v1(2), v2(2), p(2)
        integer :: i,j,n!,k

        p=S%xy
        angle=0.0
        n=size(T%xy(:,1))
        do i=1,n
            j=mod(i,n)+1
            v1=T%xy(i,:)-p
            v2=T%xy(j,:)-p
            angle = angle + dacos( dot_product(v1,v2) / sqrt( dot_product(v1,v1) * dot_product(v2,v2)) ) 
        enddo
        !NO FUNCIONA BIEN!
        inTier=(ABS(2*pi - ABS(angle)) .LT. 0.1) 
        !NO FUNCIONA BIEN!
end function

!!set to null al tier parameters
!subroutine nullifyTier(t)
!        implicit none
!        type(tier),intent(inout) :: T
!        T%hgt  =0.00
!        T%wid  =0.00
!        T%len  =0.00
!        T%L    =0.00
!        T%xbadj=0.00
!        T%ybadj=0.00
!        T%gsh  =0.00
!end subroutine

!Calculo de GEP Stack Height, XBADJ, YBADJ
subroutine calcGepStackHeight(Bz0,t,s)
        implicit none
        real :: Bz0
        type(tier),intent(inout) :: T
        type(stack), intent(inout) :: S
        T%gsh  = Bz0 + T%hgt - s%z0 + 1.5*T%L             !Equation 1 (GEP, page 6)
        T%xbadj= T%ymin - sngl(s%xy2(2))                  !XBADJ = YMIN(C) - YPSTK                 !ymin - ystack
        T%ybadj= sngl(s%xy2(1)) - (T%xmin + T%wid*0.5)    !YBADJ = XPSTK - (XMIN(C) + TW * 0.5)    !xstack- xmin + tw*0.5
end subroutine

!Buscar tiers cercanos y combinar
subroutine mergeCloseTiers(t,B,ntier,nombre)
        implicit none
        type(tier) ,intent(inout) :: T
        type(build),intent(inout) :: B(:)        
        integer,intent(in) :: ntier
        character(8),intent(in) :: nombre
        real :: dist, minL
        integer :: i
        do i=1,size(B),1
        if ( B(i)%nombre /= nombre ) then
           if ( size(B(i)%t) >= ntier ) then
              !dist = sngl( minDist( T%xy, B(i)%T(ntier)%xy ) )
              dist = sngl( minDist( T%xy2, B(i)%T(ntier)%xy2 ) )
              minL=min(t%L, B(i)%T(ntier)%L) 
              if ( dist < minL ) then
                 call CombineTiers(t, B(i)%t(ntier) )
              endif
           endif
        endif
        enddo
end subroutine

!Combinar dos tiers 
subroutine CombineTiers(t1,t2) 
        implicit none
        type(tier),intent(inout) :: T1
        type(tier),intent(inout) :: T2
        T1%xmin=min(T1%xmin, T2%xmin) 
        T1%xmax=max(T1%xmax, T2%xmax)
        T1%ymin=min(T1%ymin, T2%ymin)
        T1%ymax=max(T1%ymax, T2%ymax)
end subroutine

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

subroutine addToStkTable(sT,S,T,belev)
        implicit none
        type(stkTable), intent(inout) :: sT
        type(tier), intent(in) :: T
        type(stack),intent(in) :: S
        real,intent(in) :: belev
        st%stkName=S%nombre
        st%stkHeight=S%h
        st%BaseElevDiff=S%z0-belev
        st%GEPEQN1=T%gsh ! max(, gep)
        st%GEPSHV=max(65.0, T%gsh)
end subroutine

subroutine no_tiers_affect_this_stack(S,d,wdir,oTab,sTab)
        implicit none
        type(outTable), intent(inout) :: oTab
        type(stkTable), intent(inout) :: sTab
        type(stack),intent(in) :: S
        double precision, intent(in) :: wdir
        integer, intent(in) ::d
        oTab%stkName=S%nombre
        oTab%tabla(d,1:6)=[ sngl(wdir),0.0,0.0,0.0,0.0,0.0 ] !copy on table
        if ( .not. S%affected_by_tier) then
           WRITE(*,*) '     No tiers affect this stack.'
           sTab%stkName=S%nombre
           sTab%stkHeight=S%h
           sTab%BaseElevDiff=S%z0 !-belev
           sTab%GEPEQN1=-99.99 !T%gsh ! max(, gep)
           sTab%GEPSHV=65.0 !max(65.0, T%gsh)
        end if
end subroutine
end program
