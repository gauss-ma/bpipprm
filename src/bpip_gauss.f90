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
character(24),    parameter :: inputFileName="bpip.inp"
character(24),    parameter :: outputFileName="bpip.out"

!objects
TYPE(build), allocatable :: B(:) !array de edificios
TYPE(stack), allocatable :: S(:) !array de stacks

double precision :: wdir         !direccion del viento

logical :: L5,SIZ
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
DO d=1,36   !for each wdir (c/10 grados)
     print '("      Wind flow passing", I4," degree direction.")', d*10 
     wdir=d*10.0*deg2rad    !wdir [rad]

     call rotar(B,S,wdir)  !roto coordenadas según wdir

     DO i=1,size(S)                                    !for each stack
         refGSH=0.0; refWID=0.0; tmax=0; bmax=0
         DO j=1,size(B)                                !for each building
             DO k=1,size(B(j)%T)                       !for each tier

                !call setTierValues( B(j)%T(k) )        ! calc tier: wid len hgt L
                   
                !check si tier afecta stack (adentro de 5L y SIZ)
                SIZ=isInsideSIZ(S(i), B(j)%T(k)) 
                L5 =isInsideL5( S(i), B(j)%T(k))
                if ( L5 .AND. SIZ ) then 
                    call mergeCloseTiers(B(j)%T(k),B,k,B(j)%nombre)   !si hay otra structura cercana (de == tier) y combinarlos
                    call calcGepStackHeight(B(j)%z0, B(j)%T(k), S(i)) !calc: GSH, XBADJ, YBADJ
                else 
                   call nullifyTier(B(j)%T(k)) !set all tier's params = 0.0, (so it is not considered)
                endif
                
                !guardar el tier de máximo GSH, (si hay dos con == GSH) quedarme con el de narrower width
                if ( refGSH < B(j)%T(k)%gsh .OR. ( refGSH == B(j)%T(k)%gsh .AND. refWID >= B(j)%T(k)%wid) ) then
                     refGSH=B(j)%T(k)%gsh   
                     refWID=B(j)%T(k)%wid
                     bmax=j      !building max (index)
                     tmax=k      !tier     max (index)
                end if
             END DO!tiers
         END DO!buildings
         
         call addToOutTable(oTable(i),S(i),B(bmax)%T(tmax),d,wdir) !Agregar a Tablas
         if ( refGSH > sTable(i)%GEPEQN1 ) then
            call addToStkTable(sTable(i),S(i),B(bmax)%T(tmax), B(bmax)%z0)
         endif
     END DO!stacks
END DO!wdir

!OUTPUT:-------------------------------------------------------------
call writeOUT(oTable,sTable,title,outputFileName)

contains


subroutine setTierProyectedBoundary(T) !Calculo de XMIN XMAX YMIN YMAX
        implicit none
        type(tier),intent(inout)   :: T
        T%xmin=sngl(minval(T%xy2(:,1)) ); T%xmax=sngl(maxval(T%xy2(:,1)) )
        T%ymin=sngl(minval(T%xy2(:,2)) ); T%ymax=sngl(maxval(T%xy2(:,2)) )

        T%wid=T%xmax - T%xmin
        T%len=T%ymax - T%ymin
        T%hgt=T%h
        T%L  =min(T%wid, T%hgt)
end subroutine

!subroutine setTierValues(T)            !Calculo de WID LEN HGT L 
!        implicit none
!        type(tier),intent(inout)   :: T
!        T%wid=T%xmax - T%xmin
!        T%len=T%ymax - T%ymin
!        T%hgt=T%h
!        T%L  =min(T%wid, T%hgt)
!end subroutine

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
               call setTierProyectedBoundary(B(i)%T(j) )       ! calc tier: xmin,xmax,ymin,ymax
               !call setTierValues(           B(i)%T(j) )       ! calc tier: wid len hgt L
           enddo
        enddo
end subroutine

!Calcula distancia mínima entre dos conjuntos de puntos
function minDist(xy1,xy2) result(Dist)
        implicit none
        double precision,intent(in) :: xy1(:,:), xy2(:,:)
        double precision :: Dist
        integer :: i
        Dist=1000000.0
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
        !x_stack=S%xy2(1,1)
        !y_stack=S%xy2(1,2)
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

!set to null al tier parameters
subroutine nullifyTier(t)
        implicit none
        type(tier),intent(inout) :: T
        T%hgt  =0.00
        T%wid  =0.00
        T%len  =0.00
        T%L    =0.00
        T%xbadj=0.00
        T%ybadj=0.00
        T%gsh  =0.00
end subroutine

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

end program
