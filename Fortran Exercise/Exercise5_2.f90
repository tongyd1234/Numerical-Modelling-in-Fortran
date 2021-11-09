Module advection_diffusion_equation
implicit none
type grid
integer nx,ny
real xsize,ysize 
real B  ! flowstrength
real k  ! diffusivity
real a_dif,a_adv
real total_time
real time
real dx,dy
real dt
real,allocatable,dimension(:,:) :: T
real,allocatable,dimension(:,:,:) :: V
real,allocatable,dimension(:,:) :: S
end type grid

contains


subroutine initialize(A)
type(grid),intent(inout) :: A
integer :: i,j
real :: pi
real :: vxmax,vymax
real :: dt_dif,dt_adv
pi=4.*atan(1.)


a%dx=a%xsize*1./(a%nx-1)
a%dy=a%ysize*1./(a%ny-1)

allocate(a%T(a%ny,a%nx),a%V(2,a%ny,a%nx),a%S(a%ny,a%nx))

a%T=0
a%V=0

forall(i=1:a%ny,j=1:a%nx)
! phi=Bsin(piy/ymax)sin(pix/xmax)
a%S(i,j)=a%B*sin(pi*(i-1)*1./(a%ny-1))*sin(pi*(j-1)*1./(a%nx-1))
end forall

forall(i=2:a%ny-1,j=2:a%ny-1)
! vx=dphi/dy
a%v(1,i,j)=(-a%S(i-1,j)+a%S(i+1,j))/(2*a%dy)
! vy=-dphi/dx
a%v(2,i,j)=-(-a%S(i,j-1)+a%S(i,j+1))/(2*a%dx)
end forall

vxmax=maxval(a%v(1,:,:))
vymax=maxval(a%v(2,:,:))

dt_dif=a%a_dif*(min(a%dx,a%dy))**2/a%k
dt_adv=a%a_adv*min(a%dx/vxmax,a%dy/vymax)
a%dt=min(dt_dif,dt_adv)

a%time=0

end subroutine initialize

!2D Second-Derivative
function D2T(T,dx,dy)
real,intent(in) :: T(:,:)
real,intent(in) :: dx,dy
integer :: i,j,Nx,Ny
real,dimension(size(T,dim=1),size(T,dim=2)) :: D2T


Ny=size(T,dim=1)
Nx=size(T,dim=2)
D2T=0

forall(i=2:Ny-1,j=2:Nx-1)
D2T(i,j)=(T(i-1,j)-2*T(i,j)+T(i+1,j))/(dy**2)+(T(i,j-1)-2*T(i,j)+T(i,j+1))/(dx**2)
end forall
end function D2T

function vdT(T,V,dx,dy)
real,intent(in) :: T(:,:)
real,intent(in) :: V(:,:,:)
real,intent(in) :: dx,dy
real,dimension(size(T,dim=1),size(T,dim=2)):: vdT,vxdT,vydT

integer :: i,j,Nx,Ny
Ny=size(T,dim=1)
Nx=size(T,dim=2)
vdT=0
vxdT=0
vydT=0

!DT_DX
do i=2,Ny-1
do j=2,Nx-1
if (V(1,i,j)>0) then
! vx*dT/dx
vxdT(i,j)=V(1,i,j)*(T(i,j)-T(i,j-1))/dx
else
vxdT(i,j)=V(1,i,j)*(T(i,j+1)-T(i,j))/dx
end if 
! vy*dT/dy
if (V(2,i,j)>0) then
vydT(i,j)=V(2,i,j)*(T(i,j)-T(i-1,j))/dy
else
vydT(i,j)=V(2,i,j)*(T(i+1,j)-T(i,j))/dy
end if 
end do
end do

vdT=vxdT+vydT

end function vdT



end module advection_diffusion_equation

program main
use advection_diffusion_equation
implicit none
type (grid):: G
real :: xsize,ysize
integer :: nx,ny
real :: B
real :: k
real :: a_dif,a_adv
real :: total_time
real :: time
real,allocatable:: LaT(:,:)
real,allocatable:: vNaT(:,:)
integer :: i
namelist /inputs/ xsize,ysize,nx,ny,B,k,&
a_dif,a_adv,total_time,time
open(1,file='inputs.txt')
read(1,inputs)
close(1)
G%xsize=xsize
G%ysize=ysize
G%nx=nx
G%ny=ny
G%B=B
G%k=k
G%a_dif=a_dif
G%a_adv=a_adv
G%total_time=total_time
G%time=time

call initialize(G)
do while(G%time<G%total_time)

!Boundary Conditions
G%T(1,:)=0
G%T(G%ny,:)=1
G%T(:,1)=G%T(:,2)
G%T(:,G%nx)=G%T(:,G%nx-1)

LaT=D2T(G%T,G%dx,G%dy)
vNaT=vdT(G%T,G%V,G%dx,G%dy)
G%T=G%T+G%dt*(G%k*LaT-vNaT)
G%time=G%time+G%dt
end do

open(2,file='exercise5_1.dat')

do i=1,G%ny
write(2,'(1000(1pe12.3))'),G%T(i,:)
end do

close(2)

end program main
