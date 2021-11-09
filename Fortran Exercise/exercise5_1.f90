module advdifequation
  implicit none
  contains
     subroutine secdderivative (a,nx,ny,h,aprime)
         implicit none
         double precision,dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         double precision,dimension(nx,ny),intent(out)::aprime
         integer::i,j
         do i=2,nx-1
            do j=2,ny-1
              aprime(i,j)=(a(i+1,j)-2*a(i,j)+a(i-1,j))/(h*h)+(a(i,j+1)-2*a(i,j)+a(i,j-1))/(h*h)
            end do
         end do
         do j=1,ny,ny-1
               aprime(:,j)=0
         end do
         do i=1,nx,nx-1
               aprime(i,:)=0
         end do
     end subroutine secdderivative
!---------------------
     subroutine initialiseS(nx,ny,h,L,B,S)
         implicit none
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         real,intent(in)::L,B
         double precision,dimension(nx,ny),intent(out)::S
         integer::i,j
         real::pi=4*atan(1.0)
         do i=1,nx
             do j=1,ny
                  S(i,j)=B*sin(pi*(i-1)*h/L)*sin(pi*(j-1)*h/1.0)
             end do
        end do
     end subroutine initialiseS
!---------------------
     subroutine initialiseV(vx,vy,nx,ny,V)
         implicit none
        integer,intent(in)::nx,ny
         double precision,dimension(nx,ny),intent(in)::vx,vy
         double precision,dimension(2,nx,ny),intent(out)::V
         integer::i,j
         do i=1,nx
            do j=1,ny
               V(1,i,j)=vx(i,j)
               V(2,i,j)=vy(i,j)
            end do
         end do
     end subroutine initialiseV

!---------------------
     subroutine firdderivative_vx (a,nx,ny,h,aprime)!to calculate vx
         implicit none
         double precision,dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         double precision,dimension(nx,ny),intent(out)::aprime
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
         double precision,dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         double precision,dimension(nx,ny),intent(out)::aprime
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
         double precision,dimension(nx,ny),intent(in)::vx,vy,T
         double precision,intent(in)::h
         double precision,dimension(nx,ny)::vx_dt,vy_dt
         double precision,dimension(nx,ny),intent(out)::v_grad
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
!----------------
     subroutine findmax(v,nx,ny,vmax)
         integer,intent(in)::nx,ny
         double precision,dimension(nx,ny)::v
         double precision,intent(out)::vmax
         integer::i,j
         double precision::v_tem
         do j=1,ny
           do i=1,nx-1
              if(v(i,j)>v(i+1,j)) then
                  v_tem=v(i+1,j)
                  v(i+1,j)=v(i,j)
                  v(i,j)=v_tem
              end if
           end do
         end do
         do j=1,ny-1
             if(v(nx,j)>v(nx,j+1)) then
                  v_tem=v(nx,j+1)
                  v(nx,j+1)=v(nx,j)
                  v(nx,j)=v_tem
            end if
         end do
         vmax=v(nx,ny)
      end subroutine findmax
end module advdifequation
!------------- Programme Below--------------
program secdiffusionequation
   use advdifequation
   implicit none
      type initialise
           integer::nx,ny
           double precision,allocatable::T(:,:),S(:,:),vx(:,:),vy(:,:),V(:,:,:)
           double precision::h
      end type initialise
      type(initialise)::ETH
      integer::i,j,xgrids,ygrids
      real::B
      real::pi=4*atan(1.0)
      real::kappa,L
      real::a_dif, a_adv
      double precision::total_time,delta_t,time
      double precision::vxmax,vymax
      double precision,allocatable::v_grad(:,:),vxm(:,:), vym(:,:),dT2(:,:)
      namelist/inputs/ L,total_time,kappa,a_dif,a_adv,xgrids,ygrids,B
      open(1,file='parameters.txt',status='old')
      read(1,inputs)
      close(1)
      ETH%nx=xgrids
      ETH%ny=ygrids
      allocate(ETH%T(ETH%nx,ETH%ny),ETH%S(ETH%nx,ETH%ny),ETH%vx(ETH%nx,ETH%ny),ETH%vy(ETH%nx,ETH%ny),v_grad(ETH%nx,ETH%ny))
      allocate(vxm(ETH%nx,ETH%ny),vym(ETH%nx,ETH%ny),dT2(ETH%nx,ETH%ny),ETH%V(2,ETH%nx,ETH%ny))
      time=0
      ETH%V(:,:,:)=0
      ETH%T(:,:)=0
      ETH%T(:,1)=1
      ETH%T(:,ETH%ny)=0
      ETH%h=1.0/(ETH%ny-1)
      call initialiseS(ETH%nx,ETH%ny,ETH%h,L,B,ETH%S)
      call firdderivative_vx (ETH%S,ETH%nx,ETH%ny,ETH%h,ETH%vx)
      call firdderivative_vy (ETH%S,ETH%nx,ETH%ny,ETH%h,ETH%vy)
      call initialiseV(ETH%vx,ETH%vy,ETH%nx,ETH%ny,ETH%V)
      vxm=ETH%vx
      vym=ETH%vy
      call findmax(vxm,ETH%nx,ETH%ny,vxmax)
      call findmax(vym,ETH%nx,ETH%ny,vymax)
      delta_t=min(a_dif*ETH%h**2/kappa,a_adv*min(ETH%h/vxmax,ETH%h/vymax))
!Here dx=dy=ETH%h
      do while(time<total_time)
         ETH%T(:,1)=1
         ETH%T(:,ETH%ny)=0
         ETH%T(1,:)=ETH%T(2,:)
         ETH%T(ETH%nx,:)=ETH%T(ETH%nx-1,:)
         call secdderivative(ETH%T,ETH%nx,ETH%ny,ETH%h,dT2)
         call v_grad_cal(ETH%vx,ETH%vy,ETH%T,ETH%h,ETH%nx,ETH%ny,v_grad)
         ETH%T=ETH%T+delta_t*(kappa*dT2-v_grad)
         time=time+delta_t
      end do
      open(3,file='exercise5_2.dat')
        do j=1,ETH%ny
        write(3,'(1000(1pe13.5))') ETH%S(:,j)
        end do
      close(3)
      open(2,file='exercise5_1.dat')
        do j=1,ETH%ny
         write(2,'(1000(1pe13.5))') ETH%T(:,j)
        end do
      close(2)
      deallocate(ETH%T,ETH%S,ETH%vx,ETH%vy,v_grad,vxm,vym,dT2)
end program secdiffusionequation
