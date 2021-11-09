module mean_std
     implicit none
     real::SD
     contains
             real function meanArray(a,n,SD)
              implicit none
              integer,intent(in)::n
              real,intent(in)::a(n)
              integer::i
              real::sum=0, sum_squared=0
              real::mean
              real,intent(out)::SD
              do i=1,n
                 sum=sum+a(i)
                 sum_squared=sum_squared+a(i)**2
              end do
              mean=sum/n
              SD=sqrt((sum_squared/n)-mean*mean)
              meanArray=mean
          end function meanArray
    subroutine derivative (a,np,h,aprime)
         implicit none
         integer,intent(in)::np
         double precision,intent(in)::a(np),h
         double precision,intent(out)::aprime(np)
         integer::i
         do i=1,np-2
            aprime(i)=(a(i+2)-2*a(i+1)+a(i))/(h*h)
         end do
     end subroutine derivative
end module mean_std

!program test
!end program test
