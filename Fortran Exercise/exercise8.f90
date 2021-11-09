module Poisson_equation
implicit none
integer, parameter :: dp = selected_real_kind(15, 307)
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
real(kind=dp),intent(in)::coarse(:,:)
real(kind=dp),intent(out)::fine(:,:)
integer::i,j,nx,ny,nxc,nyc
    nxc = size(coarse,1);	nyc= size(coarse,2)
    nx	= 2*nxc-1; 			ny= 2*nyc-1
	forall(i=1:nxc,j=1:nyc) fine(2*i-1,2*j-1) = coarse(i,j)
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
    real (kind=dp)       	:: alpha=1
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
!---------------------------
!---------------------------
module advdifequation
  implicit none
  integer, parameter :: db = selected_real_kind(15, 307)
  contains
!---------------------
     subroutine firdderivative_vx (a,nx,ny,h,aprime)!to calculate vx
         implicit none
         real(kind=db),dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         real(kind=db),intent(in)::h
         real(kind=db),dimension(nx,ny),intent(out)::aprime
         integer::i,j
         do i=2,nx-1
            do j=1,ny
              aprime(i,j)=(a(i,j+1)-a(i,j-1))/(2*h)
            end do
         end do
         do i=1,nx,nx-1
              aprime(i,:)=0
         end do
     end subroutine firdderivative_vx
!---------------------
     subroutine firdderivative_vy (a,nx,ny,h,aprime)!to calculate vy
         implicit none
         real(kind=db),dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         real(kind=db),intent(in)::h
         real(kind=db),dimension(nx,ny),intent(out)::aprime
         integer::i,j
         do j=2,ny-1
            do i=1,nx
              aprime(i,j)=-(a(i+1,j)-a(i-1,j))/(2*h)
            end do
         end do
         do j=1,ny,ny-1
              aprime(:,j)=0
         end do
     end subroutine firdderivative_vy
!-----------------
     subroutine v_grad_cal(vx,vy,T,h,nx,ny,v_grad)
         implicit none
         integer,intent(in)::nx,ny
         real(kind=db),dimension(nx,ny),intent(in)::vx,vy,T
         real(kind=db),intent(in)::h
         real(kind=db),dimension(nx,ny)::vx_dt,vy_dt
         real(kind=db),dimension(nx,ny),intent(out)::v_grad
         integer::i,j
         do i=2,nx-1
           do j=2,ny-1
              if (vx(i,j)>=0) then
                    vx_dt(i,j)=vx(i,j)*(T(i,j)-T(i-1,j))/h
              end if
              if (vx(i,j)<0) then
                    vx_dt(i,j)=vx(i,j)*(T(i+1,j)-T(i,j))/h
              end if
              if (vy(i,j)>=0) then
                    vy_dt(i,j)=vy(i,j)*(T(i,j)-T(i,j-1))/h
              end if
              if (vy(i,j)<0) then
                    vy_dt(i,j)=vy(i,j)*(T(i,j+1)-T(i,j))/h
              end if
           end do
         end do
         do i=1,nx,nx-1
             v_grad(i,:)=0
         end do
         do j=1,ny,ny-1
             v_grad(:,j)=0
         end do
         v_grad=vx_dt+vy_dt
     end subroutine v_grad_cal
!--------------------------------
     subroutine rhside(T,nx,ny,h,Ray_Num,rhs)
        implicit none
         integer,intent(in)::nx,ny
         real(kind=db),dimension(:,:),intent(in)::T
         real(kind=db),intent(in)::h,Ray_Num
         real(kind=db),dimension(nx,ny),intent(out)::rhs
         integer::i,j
         do j=2,ny-1
            do i=2,nx-1
              rhs(i,j)=-(T(i+1,j)-T(i-1,j))/(2*h)*(Ray_Num)
            end do
         end do
         do i=1,nx,nx-1
              rhs(i,:)=0
         end do
         do j=1,ny,ny-1
              rhs(:,j)=0
         end do
     end subroutine rhside
end module advdifequation
!-----------Main Programme Below---------------!
program convection
   use advdifequation
   use Poisson_equation
   implicit none
   integer::nx,ny,i,j,L,a
   real(kind=dp)::a_dif,a_adv,total_time,h,dt_dif,dt_adv,time,Ray_Num,dt,alpha,Pr
   real(kind=dp)::S_rms,T_rms,W_rms,rhs_rms,t1
   real(kind=dp),allocatable::T(:,:),S(:,:),W(:,:),vx(:,:),vy(:,:),rhs(:,:),v_grad(:,:),w_grad(:,:)

   namelist /inputs/ L,total_time,a_dif,a_adv,nx,ny,Ray_Num,Pr
   open(1,file='parameters_Prandtl.txt',status='old')
   read(1,inputs)
   close(1)
   print*,'My input parameters are:'
   write(*,inputs)
   allocate(T(nx,ny),S(nx,ny),W(nx,ny),vx(nx,ny),vy(nx,ny),rhs(nx,ny),v_grad(nx,ny),w_grad(nx,ny))
! Initialise parameters
   h=1.0/(ny-1)
   dt_dif=a_dif*h*h/max(Pr,1.0)
   alpha=1
   time=0
!initialise T
  do 
     write(*,*) 'Enter 1 for spike T source and 2 for random T source: '
     read(*,*) a
     if ((a>=1).and.(a<=2)) then
        exit
     end if
  end do
  if (a==1) then
   print*,'You are calling the spike T source'
   T(nx/2+1,ny/2+1)=1./h**2
  end if
  if (a==2) then
   print*,'You are calling the random T source'
   call random_number(T)
  end if
!Initialse W
  do 
     write(*,*) 'Enter 1 for spike Omiga source and 2 for random Omiga source: '
     read(*,*) a
     if ((a>=1).and.(a<=2)) then
        exit
     end if
  end do
  if (a==1) then
   print*,'You are calling the spike Omiga source'
   W(nx/2+1,ny/2+1)=1./h**2
  end if
  if (a==2) then
   print*,'You are calling the random Omiga source'
   call random_number(W)
  end if
  print*,'Calculating, please wait.'
   T(:,1)=1
   T(:,ny)=0
   T(1,:)=T(2,:)
   T(nx,:)=T(nx-1,:)
   S=0
   vx=0
   vy=0
   i=0
   W=-W
   W_rms=sqrt(sum(W*W)/(nx*ny))
   S_rms=W_rms/10
!-----Run the do loop -----------------

 do while(time<=total_time)
!Boundary condition for T
   T(:,1)=0
   T(:,ny)=1
   T(1,:)=T(2,:)
   T(nx,:)=T(nx-1,:)
   W_rms=sqrt(sum(W*W)/(nx*ny))
!Boundary condition for W
   W(1,:)=0
   W(:,1)=0
   W(:,ny)=0
   W(nx,:)=0
   S_rms=W_rms/10  !Set s_rms to make sure the loop would run
   do
      if (S_rms/W_rms<=10.**(-3)) exit
      S_rms=Vcycle_2DPoisson(S,W,h)
!      print*,S_rms
   end do
!Boundary condition for W
   S(1,:)=0
   S(:,1)=0
   S(:,ny)=0
   S(nx,:)=0
!Calculate Velocity both in x and y
   call firdderivative_vx (S,nx,ny,h,vx)
   call firdderivative_vy (S,nx,ny,h,vy)
!Determine dt
   dt_adv=a_adv*min(h/maxval(abs(vx)),h/maxval(abs(vy)))
   dt=min(dt_dif,dt_adv)
   call v_grad_cal(vx,vy,T,h,nx,ny,v_grad)
   call v_grad_cal(vx,vy,W,h,nx,ny,w_grad)
   call rhside(T,nx,ny,h,Ray_Num,rhs)
   T=T+dt*(secdderivative(T,h)-v_grad)
   W=W+dt*(Pr*secdderivative(W,h)-w_grad-Pr*rhs)
   time=time+dt
end do
   open(1,file='exercise8_3.dat')
      do j=1,ny
      write(1,'(1000(1pe13.5))') W(:,j)
      end do
   close(1)
   open(2,file='exercise8_1.dat')
      do j=1,ny
      write(2,'(1000(1pe13.5))') S(:,j)
      end do
   close(2)
   open(3,file='exercise8_2.dat')
      do j=1,ny
      write(3,'(1000(1pe13.5))') T(:,j)
      end do
   close(3)
   deallocate(T,S,W,vx,vy,rhs,v_grad,w_grad)
   call cpu_time(t1)
   print*,'It takes',t1,'seconds to calculate.'
end program convection
