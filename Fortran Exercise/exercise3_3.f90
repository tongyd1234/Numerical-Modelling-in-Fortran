module dderivative
  implicit none
  contains
     subroutine derivative (a,np,h,aprime)
         implicit none
         integer,intent(in)::np
         double precision,intent(in)::a(np),h
         double precision,intent(out)::aprime(np)
         integer::i
         do i=1,np-2
            aprime(i)=(a(i+1)-2*a(i)+a(i-1))/(h*h)
         end do
     end subroutine derivative
end module dderivative

program diffusionequation
   use dderivative
   implicit none
      real::L,total_time,time
      integer::N
      integer::i=1,j=1
      double precision,allocatable::T(:),d2Tdx2(:)
      double precision::x,dx,dt
      write(*,'(a,$)')'Input length of domain:'
      read*,L
      write(*,'(a,$)')'Input number of grid points:'
      read*,N
      write(*,'(a,$)')'Input total time:'
      read*,total_time
      dx=L/(N-1)
      dt=0.1*dx**2
      allocate (T(n),d2Tdx2(n))
      T(:)=0
      T(N/2)=1.0
      print*,T
      open(2,file='exercise3_3.xls')
      do i=1,N
         write(2,*) T(i)
      end do
      close(2)
      do while(time<total_time)
      call derivative(T,N,dx,d2Tdx2)
         T=T+dt*d2Tdx2
         time=time+dt
      end do
      print*,d2Tdx2
      print*,T
      open(2,file='exercise3_31.xls')
      do i=1,N
         write(2,*) T(i)
      end do
      close(2)
      deallocate(T,d2Tdx2)
end program diffusionequation
