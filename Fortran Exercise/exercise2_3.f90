program exercise2_3
      implicit none 
      integer::n,i
      double precision,allocatable::y(:),d2ydx2(:)
      double precision::x,dx
      write(*,'(a,$)')'Input number of grid points:'
      read*,n
      allocate (y(n),d2ydx2(n))
      dx=10./(n-1)
      do i=1,n
         x=(i-1)*dx
         y(i)=10*x**2
      end do
!      print*,y
      call derivative (y,n,dx,d2ydx2)
      do i=1,n-2   !Actually we could just work out the 2nd derivative of the first n-2 numbers according to the formation below. aprime(i)=(a(i+2)-2*a(i+1)+a(i))/(h*h). So it would be better if we just output those n-2 numbers. 
         x=(i-1)*dx
         print*, d2ydx2(i),20,20-d2ydx2(i)
      end do
      deallocate(y,d2ydx2)
  contains
      subroutine derivative (a,np,h,aprime)
         integer,intent(in)::np
         double precision,intent(in)::a(np),h
         double precision,intent(out)::aprime(np)
         integer::i
         do i=1,np-2
            aprime(i)=(a(i+2)-2*a(i+1)+a(i))/(h*h)
         end do
      end subroutine derivative
end program exercise2_3
