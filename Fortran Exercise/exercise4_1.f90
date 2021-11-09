module dderivative
  implicit none
  contains
     function derivative (a,h)
         implicit none
         double precision,dimension(:),intent(in)::a
         double precision,intent(in)::h
         real,dimension(size(a))::derivative
         integer::i
         do i=2,size(a)-1
            derivative(i)=(a(i+1)-2*a(i)+a(i-1))/(h*h)
         end do
         derivative(1)=0
         derivative(size(a))=0
     end function derivative
end module dderivative

program diffusionequation
   use dderivative
   implicit none
      real::L,total_time,time
      integer::N
      integer::i=1,j=1,k=1
      double precision,allocatable::T(:)
      double precision::x,dx,dt
      write(*,'(a,$)')'Input length of domain:'
      read*,L
      write(*,'(a,$)')'Input number of grid points:'
      read*,N
      write(*,'(a,$)')'Input total time:'
      read*,total_time
      dx=L/(N-1)
      dt=0.1*dx**2
      allocate (T(n))
      T(:)=0
      T(N/2)=1.0
!      print*,T
      open(2,file='exercise3_3.xls')
      do i=1,N
         write(2,*) T(i)
      end do
      close(2)
      do while(time<total_time)
         T=T+dt*(derivative(T,dx))
         time=time+dt
      end do
!      print*,T   !to test if it worked
      open(2,file='exercise3_31.xls')
      do i=1,N
         write(2,*) T(i)
      end do
      close(2)
      deallocate(T)
end program diffusionequation
