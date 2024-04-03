program bpip_gauss
!MODULES-------------------------------------------------------------
use struc            !objetos
use readBPIP         !lectura de datos de entrada
use writeBPIP        !escritura de salidas

!VARIABLES-----------------------------------------------------------
implicit none
integer :: i,j,k,d           !indices(para: stacks,buildings,tiers,wdirs)

!global
double precision, parameter :: pi=3.141592653589793_8
double precision, parameter :: deg2rad=2*pi/360.0
character(24),    parameter :: inputFileName="bpip.inp"
character(24),    parameter :: outputFileName="bpip.out"

!objects
TYPE(build), allocatable :: B(:) !array de edificios
TYPE(stack), allocatable :: S(:) !array de stacks

double precision :: wdir         !direccion del viento

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
DO d=1,36,1 !for each wdir (c/10 grados)
     print '("      Wind flow passing", I4," degree direction.")', d*10 
     wdir=(360-d*10)*deg2rad   !wdir [rad]

     call rotar(B,S,wdir)  !rotar coordenadas según wdir

     DO i=1,size(S),1 !for each stack
         refGSH=0.0; refWID=0.0; tmax=0; bmax=0
         DO j=1,size(B),1 !for each building
             DO k=1,size(B(j)%t),1 !for each tier
                    
                call setTierProyectedBoundary(B(j)%T(k))  ! xmin,xmax,ymin,ymax
                call setTierValues(B(j)%T(k))             ! tier: wid len hgt L

                !check si tier afecta stack (adentro de 5L y SIZ)
                if ( isInsideL5(S(i),B(j)%T(k)) .AND. isInsideSIZ(S(i), B(j)%T(k)) ) then
                   !check si hay otra structura cercana (de mismo tier) y combinarlos:    
                    call mergeCloseTiers(B(j)%T(k),B,k,B(j)%nombre)
                    call calcGepStackHeight(B(j)%z0, B(j)%T(k), S(i)) !calc: GSH, XBADJ, YBADJ
                else 
                   !if ( .not. isInsideSIZ(S(i), B(j)%T(k)) ) then
                   !   print*,"Stack not in SIZ!", S(i)%nombre,B(j)%nombre
                   !endif
                   !set all tier's params = 0.0
                   call nullifyTier(B(j)%T(k))
                endif
                
                !guardar el tier de máximo GSH, (si hay dos con == GSH) quedarme con el de narrower width
                if ( refGSH < B(j)%T(k)%gsh .OR. ( refGSH == B(j)%T(k)%gsh .AND. refWID >= B(j)%T(k)%wid) ) then
                     refGSH=B(j)%T(k)%gsh       
                     refWID=B(j)%T(k)%wid
                     bmax=j      !index building max
                     tmax=k      !index tier     max
                end if
     
             END DO!tiers
         END DO!buildings
          
         !Agregar a Tablas       
         call addToOutTable(oTable(i),S(i),B(bmax)%T(tmax),d,wdir)

         if (refGSH > sTable(i)%GEPEQN1 ) then
            call addToStkTable(sTable(i),S(i),B(bmax)%T(tmax), B(bmax)%z0)
         endif
     END DO!stacks
END DO!wdir

!OUTPUT:-------------------------------------------------------------
call writeOUT(oTable,sTable,title,outputFileName)

contains
        
!Calcula distancia mínima entre dos conjuntos de puntos
function minDist(xy1,xy2) result(Dist)
        implicit none
        double precision,intent(in) :: xy1(:,:), xy2(:,:)
        double precision :: Dist
        integer :: i
        Dist=1000000.0
        do i=1,size(xy1,1),1
           Dist=min(Dist, sqrt( minval( ( xy2(:,1) - xy1(i,1) )**2 + ( xy2(:,2) - xy1(i,2) )**2 ) ) )
        end do
end function

function isInsideL5(S,T)   result(L5)
        implicit none
        logical :: L5
        type(stack), intent(in) :: S
        type(tier), intent(in)  :: T
        !  Distancia Stack Tier debe ser menor a 5L
        L5=( minDIST( S%xy, T%xy ) .LE. 5.0*T%L )
end function
function isInsideSIZ(S,T)       result(SIZ)
        implicit none
        type(stack), intent(in) :: S
        type(tier), intent(in)  :: T
        logical :: SIZ
        double precision:: x_stack,y_stack
        x_stack=S%xy2(1,1)
        y_stack=S%xy2(1,2)
        !  Stack está dentro del SIZ de la estructura?  SIZ= [xmin-0.5*L , xmax + 0.5*L] , [ymin -2L , ymax+5L ]
        SIZ=(           x_stack .GE. (T%xmin - 0.5*T%L)  )       !*sign(T%xmin,1.0)                                   
        SIZ=(SIZ .AND. (x_stack .LE. (T%xmax + 0.5*T%L) ))       !*sign(T%xmax,1.0)
        SIZ=(SIZ .AND. (y_stack .GE. (T%ymin - 2.0*T%L) ))       !*sign(T%ymin,1.0)
        !print*,"Stack xy:",S%xy2(1,:)
        !print*,"Tier min, max, L:",T%xmin, T%xmax,T%ymin,T%ymax,T%L
        !print*,"SIZ=",SIZ
endfunction

!Calculo de nuevas coordenadas
subroutine rotar(B,S,wdir)
        implicit none
        type(build),intent(inout)   :: B(:)
        type(stack), intent(inout)  :: S(:)
        double precision,intent(in) :: wdir
        double precision :: R(2,2)
        integer :: i,j!,k
        !matriz de rotación:       
        R(1,1)=dcos(wdir); R(2,2)= R(1,1) 
        R(1,2)=dsin(wdir); R(2,1)=-R(1,2) 
        !stack
        do i=1,size(S),1
           S(i)%xy2(1,:)=matmul(R,S(i)%xy(1,:))  !proyected coords stack
        enddo
        !buildings
        do i=1,size(B),1
           do j=1,size(B(i)%t),1
              B(i)%t(j)%xy2(:,1:2)=transpose(  matmul(R,transpose(B(i)%t(j)%xy(:,1:2)) ) ) !proyected coords tier
              call setTierProyectedBoundary(B(i)%T(j)) !(?) va acá o puede salir del loop?
           enddo
        enddo
end subroutine

!Calculo de XMIN XMAX YMIN YMAX
subroutine setTierProyectedBoundary(T)
        implicit none
        type(tier),intent(inout)   :: T
        T%xmin=sngl(minval(T%xy2(:,1)) ); T%xmax=sngl(maxval(T%xy2(:,1)) )
        T%ymin=sngl(minval(T%xy2(:,2)) ); T%ymax=sngl(maxval(T%xy2(:,2)) )
end subroutine

!Calculo de WID LEN HGT L 
subroutine setTierValues(T)
        implicit none
        type(tier),intent(inout)   :: T
        T%wid=T%xmax - T%xmin
        T%len=T%ymax - T%ymin
        T%hgt=T%h
        T%L=min(T%wid, T%hgt)
end subroutine

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
        !T%xy2=0.0       !coords proyectadas
end subroutine

!Calculo de GEP Stack Height, XBADJ, YBADJ
subroutine calcGepStackHeight(Bz0,t,s)!xy,xy2,h,
        implicit none
        real :: Bz0
        type(tier),intent(inout) :: T
        type(stack), intent(inout) :: S
        !calculo widht height L, gepStackHeight para este tier
        T%gsh  = Bz0 + T%hgt - s%z0 + 1.5*T%L       !Equation 1 (GEP, page 6)
        T%xbadj= T%ymin - sngl(s%xy2(1,2))                !XBADJ = YMIN(C) - YPSTK                 !ymin - ystack
        T%ybadj= sngl(s%xy2(1,1)) - (T%xmin + T%wid*0.5)  !YBADJ = XPSTK - (XMIN(C) + TW * 0.5)    !xstack- xmin + tw*0.5
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
              dist = sngl( minDist( T%xy, B(i)%T(ntier)%xy ) )
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
        type(tier),intent(inout)    :: T2
        T1%xmin=min(T1%xmin, T2%xmin) 
        T1%xmax=max(T1%xmax, T2%xmax)
        T1%ymin=min(T1%ymin, T2%ymin)
        T1%ymax=max(T1%ymax, T2%ymax)
        !T1%hgt =max(T1%hgt , T2%hgt )
end subroutine

subroutine addToOutTable(table,S,T,d,wdir)
        implicit none
        type(outTable), intent(inout) :: table
        type(tier), intent(in) :: T
        type(stack),intent(in) :: S
        double precision, intent(in) :: wdir
        integer, intent(in) ::d
        table%stkName=S%nombre
        table%tabla(d,1:6)=[ sngl(wdir),t%hgt,t%wid,t%len,t%xbadj,t%ybadj ] !copiar en tabla:
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
