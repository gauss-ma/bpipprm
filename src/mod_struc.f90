module struc
!implicit none

type tier
        !integer :: nn   !# nodes
        real,allocatable :: xy(:,:)
        !proyectados:
        real,allocatable :: xy2(:,:)
        real :: hgt,wid,len   !height,width,length (cambia para c/wdir)
        real :: L             !L: min{ wid , hgt } (cambia para c/wdir)
        real :: xmin,xmax,ymin,ymax!               (cambia para c/wdir)
        real :: gsh           !gep stack height    (cambia para c/stack y wdir)
        real :: xbadj,ybadj   !gep stack height    (cambia para c/stack y wdir)
endtype

type build
  !integer :: nt                        !# tiers
  character(8) :: nombre                !name
  real :: z0                            !base height
  type(tier),allocatable :: t(:)        !tier
endtype 

!'CALD1' 9.0 4.5325 371654.70 6174202.42
type stack
        character(8) :: nombre
        real :: z0, h
        real :: xy(2),xy2(2)
endtype

end module
