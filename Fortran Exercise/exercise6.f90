module Poisson_equation
implicit none
integer, parameter :: dp = selected_real_kind(15, 307)
!-----------------
contains
!---------------
function secdderivative(V,h) result(f)
implicit none
real(kind=dp),intent(in) :: V(:,:)
real(kind=dp),intent(in) :: h
integer :: i,j,nx,ny
real(kind=dp),dimension(size(V,dim=1),size(V,dim=2)) :: f
nx=size(V,dim=1)
ny=size(V,dim=2)
f=0
forall(i=2:nx-1,j=2:ny-1)
f(i,j)=(V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1)-4*V(i,j))/(h*h)
end forall
f(1,:)=0
f(nx,:)=0
f(:,1)=0
f(:,ny)=0
end function secdderivative
!---------------
real(dp) function iteration_2DPoisson(u,f,h,alpha)
implicit none
real(kind=dp),intent(in)::f(:,:)
real(kind=dp),intent(in)::h,alpha
real(kind=dp),intent(inout)::u(:,:)
real(kind=dp),dimension(size(f,dim=1),size(f,dim=2))::res
integer::i,j,nx,ny
nx=size(f,dim=1)
ny=size(f,dim=2)
res(1,:)=0
res(nx,:)=0
res(:,1)=0
res(:,ny)=0
u(1,:)=0
u(nx,:)=0
u(:,1)=0
u(:,ny)=0
do i=2,nx-1
   do j=2,ny-1
      res(i,j)=(u(i-1,j)+u(i+1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/h**2-f(i,j)
      u(i,j)=u(i,j)+alpha*res(i,j)*h*h/4
   end do
end do
iteration_2DPoisson=sqrt(sum(res**2)/(nx*ny))
end function iteration_2DPoisson
!---------------
subroutine residue_2DPoisson(u,f,h,res)
implicit none
real(kind=dp),intent(in)::u(:,:),f(:,:)
real(kind=dp),intent(in)::h
real(kind=dp),dimension(size(f,dim=1),size(f,dim=2)),intent(out)::res
real(kind=dp),dimension(size(f,dim=1),size(f,dim=2))::u_d2
integer::i,j,nx,ny
u_d2=secdderivative(u,h)
nx=size(u,dim=1)
ny=size(u,dim=2)
forall(i=2:nx-1,j=2:ny-2)
res(i,j)=u_d2(i,j)-f(i,j)
end forall
res(1,:)=0
res(nx,:)=0
res(:,1)=0
res(:,ny)=0
end subroutine residue_2DPoisson
!---------------
subroutine restrict(fine,coarse)
implicit none
real(kind=dp),intent(in)::fine(:,:)
real(kind=dp),intent(out)::coarse(:,:)
integer::i,j,nx,ny
nx=(size(fine,dim=1)-1.0)/2+1
ny=(size(fine,dim=2)-1.0)/2+1
do i=1,nx
   do j=1,ny
      coarse(i,j)=fine(2*i-1,2*j-1)
   end do
end do
end subroutine restrict
!---------------
subroutine prolongate(coarse,fine)
implicit none
integer::nx,ny,nxc,nyc,i,j
real(kind=dp),intent(in)::coarse(:,:)
real(kind=dp),intent(out)::fine(:,:)
nxc=size(coarse,1)	
nyc=size(coarse,2)
nx=2*nxc-1
ny=2*nyc-1
forall(i=1:nxc,j=1:nyc) 
   fine(2*i-1,2*j-1) = coarse(i,j)
end forall
    do i=2,nx-1,2
      do j=1,ny,2
        fine(i,j) = (fine(i+1,j)+fine(i-1,j))/2.
      end do
	end do 
    do i=1,nx,2
      do j=2,ny-1,2
    	fine(i,j) = (fine(i,j+1)+fine(i,j-1))/2.
      end do
    end do
	do i=2,nx-1,2	! this loop must be the last one
      do j=2,ny-1,2
        fine(i,j) = (fine(i+1,j+1)+fine(i-1,j-1)+fine(i+1,j-1)+fine(i-1,j+1))/4.
      end do
    end do
end subroutine prolongate
!---------------
  recursive function Vcycle_2DPoisson(u_f,rhs,h) result (resV)
    implicit none
    real resV
    real(kind=dp),intent(inout):: u_f(:,:)  ! arguments
    real(kind=dp),intent(in)   :: rhs(:,:),h
    real(kind=dp),allocatable  :: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)    
    integer		:: nx, ny, nxc,nyc, i  ! local variables
    real (kind=dp)       	:: alpha=0.8
    real(kind=dp)			::res_rms

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size
    if(nxc*2-1 /= nx .or. nyc*2-1 /= ny) stop "ERROR: not possible to divide grid by 2"

    if (min(nx,ny)>5) then  ! not the coarsest level

       allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))

       !---------- take 2 iterations on the fine grid--------------
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha) 
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)

       !---------- restrict the residue to the coarse grid --------
       call residue_2DPoisson(u_f,rhs,h,res_f) 
       call restrict(res_f,res_c)
       !print*,res_rms
       !---------- solve for the coarse grid correction -----------
       corr_c = 0.  
       res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2) ! *RECURSIVE CALL*
       !print*,res_rms
       !---- prolongate (interpolate) the correction to the fine grid 
       call prolongate(corr_c,corr_f)
       !---------- correct the fine-grid solution -----------------
       u_f = u_f - corr_f  
       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       deallocate(res_f,corr_f,res_c,corr_c)

    else  

       !----- coarsest level (ny=5): iterate to get 'exact' solution

       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha)
       end do

    end if
    resV = res_rms   ! returns the rms. residue
  end function Vcycle_2DPoisson
end module Poisson_equation
!---------Programme below----------!
program TwoDPoisson_Solver
   use Poisson_equation
   implicit none
integer::nx,ny,i,j,source, iteration
real(kind=dp)::alpha
real(kind=dp)::h,r_rms,rvc_rms,f_rms
real(kind=dp),allocatable::u(:,:),f(:,:)

namelist/inputs/ nx, ny, alpha, source, iteration
open(1,file='parameter_Poisson.txt')
read(1,inputs)
close(1)
print*,'My input parameters are:'
write(*,inputs)

allocate(f(nx,ny),u(nx,ny))
h=1./(ny-1)
u=0
if (source==1) then
f(:,:)=0
f(nx/2+1,ny/2+1)=1/h**2
print*,'This is a spike source.'
end if
if (source==2) then
call random_number(f)
print*,'This is a source with random number.'
end if
   open(3,file='exercise6_1.dat')
    do i=1,ny
        write(3,'(1000(1pe13.5))') f(:,i)
    end do
    close(3)
i=0
j=0
f_rms=sqrt(sum(f*f)/(nx*ny))
rvc_rms=100
if (iteration==1) then
  print*,'This is the iteration method'
do 
  r_rms=iteration_2DPoisson(u,f,h,alpha)
  if (r_rms/f_rms<=10.**(-5)) exit
  i=i+1
end do
print*,'Iteration number is :',i
print*,r_rms
print*,'The u_rms/f_rms is :',r_rms/f_rms
end if
if (iteration==2) then
   print*,'This is the vcycle method'
do
   if (rVc_rms/f_rms<=10.**(-5)) exit
   rVc_rms=Vcycle_2DPoisson(u,f,h)
   j=j+1
end do
print*,'Vcycle number is ',j 
print*,rvc_rms
print*,'The u_rms/f_rms is :',rvc_rms/f_rms
end if
   open(2,file='exercise6.dat')
    do i=1,ny
        write(2,'(1000(1pe13.5))') u(:,i)
    end do
    close(2)
deallocate(u,f)
end program TwoDPoisson_Solver
