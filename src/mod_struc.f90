module struc
!implicit none

type tier
        !integer :: nn   !# nodes
        double precision,allocatable :: xy(:,:)    !coords 
        real :: h                                  !height  (leido del .inp) (NO CAMBIA)
        !proyectados:
        double precision,allocatable :: xy2(:,:)   !coords (cambia para c/wdir)
        real :: xmin,xmax,ymin,ymax!               (cambia para c/wdir)
        real :: hgt,wid,len   !height,width,length (cambia para c/wdir)
        real :: L             !L: min{ wid , hgt } (cambia para c/wdir)
        real :: gsh           !gep stack height    (cambia para c/stack y wdir)
        real :: xbadj,ybadj   !gep stack height    (cambia para c/stack y wdir)
endtype

type build
  !integer :: nt                        !# tiers
  character(8) :: nombre                !name
  real :: z0                            !base height
  type(tier),allocatable :: t(:)        !tier
endtype 

type stack
        character(8) :: nombre
        real :: z0, h
        !real :: xy(1,2),xy2(1,2)
        double precision :: xy(1,2),xy2(1,2)
endtype

!TABLAS:
type outTable
     character(8) ::stkName  
     real :: tabla(36,6)        !36 windirs, 6vars: wdir,hgt,wid,len,xbadj,ybadj
end type

type stkTable
     character(8) ::stkName 
     real :: stkHeight, BaseElevDiff, GEPEQN1=0.0, GEPSHV
end type



end module
