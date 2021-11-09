module delta_square
  implicit none
  contains
     function derivative (a,nx,ny,h)
         implicit none
         double precision,dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         real,dimension(nx,ny)::derivative
         integer::i,j
         do i=2,nx-1
            do j=2,ny-1
              derivative(i,j)=(a(i+1,j)-2*a(i,j)+a(i-1,j))/(h*h)+(a(i,j+1)-2*a(i,j)+a(i,j-1))/(h*h)
            end do
         end do
         do j=1,ny,ny-1
            do i=1,nx
               derivative(i,j)=0
            end do
         end do
         do i=1,nx,nx-1
            do j=1,ny
               derivative(i,j)=0
            end do
         end do
     end function derivative
end module delta_square

program secdiffusionequation
   use delta_square
   implicit none
      real::L,total_time,time,kappa,a
      integer::N,M
      integer::i,j
      double precision,allocatable::T(:,:)
      double precision::x,dx,dt,dy,h
      namelist/inputs/ L,total_time,kappa,a,N,M,i,j
      open(1,file='parameters.txt',status='old')
      read(1,inputs)
      close(1)
!      write(*,inputs)
!      write(*,'(a,$)')'Input length of domain（assume equal in both x and y axis）:'
!     read*,L
!      write(*,'(a,$)')'Input number of grid points in x-direction:'
!      read*,N
!      write(*,'(a,$)')'Input number of grid points in y-direction:'
!      read*,M
!      write(*,'(a,$)')'Input total time:'
!      read*,total_time
!      write(*,'(a,$)')'Input thermal diffusivity:'
!      read*,kappa
!      write(*,'(a,$)')'Input time step constant:'
!      read*,a
      h=L/(N-1)
      dx=L/(M-1)
      dt=a*dx**2/kappa
      allocate (T(N,M))
      T(:,:)=0
      T(N/2,M/2)=1.0
!      print*,T
      open(2,file='exercise4_3.dat')
        do j=1,M
         write(2,'(1000(1pe13.5))') T(:,j)
        end do
      close(2)
      do while(time<total_time)
         T=T+dt*kappa*derivative(T,N,M,h)
         time=time+dt
      end do
!      print*,T
      open(2,file='exercise4_31.dat')
        do j=1,M
         write(2,'(1000(1pe13.5))') T(:,j)
        end do
      close(2)
      deallocate(T)
end program secdiffusionequation
