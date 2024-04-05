module readBPIP
        use struc
contains

subroutine readINP(inp_file,B,S,title) 
        implicit none
        character(78),intent(inout) :: title
        character(24),intent(in) :: inp_file
        TYPE(build),allocatable, intent(inout) :: B(:)
        TYPE(stack),allocatable, intent(inout) :: S(:)
        integer :: nb,nt,nn,ns !# builds,# tiers # nodes,# stacks
        integer :: i,j,k

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
                !read(1,*) S(i)%nombre, S(i)%z0, S(i)%h, S(i)%xy(1,1), S(i)%xy(1,2)    !name z0 h x y
                read(1,*) S(i)%nombre, S(i)%z0, S(i)%h, S(i)%xy(1), S(i)%xy(2)    !name z0 h x y
                !print*,"STACK: ", S(i)%nombre, S(i)%z0,S(i)%h,S(i)%xy(1,:)          !debug
        end do

        close(1)

endsubroutine
end module
