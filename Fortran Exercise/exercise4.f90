program exercise4
      implicit none
      integer::n,i
      real::sum=0, sum_squared=0,mean, m
      double precision:: SD 
      print*,"How many numbers do you want to input(n>=1 and n has to be an integer)?"
      read*,n
      if (n<1) then
         print*, "YOU CANNOT INPUT LESS THEN ONE NUMBER!!!"
      else
         do i=1,n
           print*,"Please enter them one by one.(",n-i+1,"number(s) left)"
           read*, m
           sum=sum+m
           sum_squared=sum_squared+m*m
         end do
         mean=sum/n
         SD=sqrt((sum_squared/n)-mean*mean)
         print*,"The mean of those number you input is", mean
         print*,"The Standard Deviation of those number you input is", SD
      end if
end program exercise4
