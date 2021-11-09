program exercise2_2
      implicit none
      integer::i
      real:: SD 
      real,allocatable::f(:)
      integer m(1)
      real,external::meanArray
      print*,"How many numbers do you want to input(n>=1 and n has to be an integer)?"
      read*,m(1)
      if (m(1)<1) then
         print*, "YOU CANNOT INPUT LESS THEN ONE NUMBER!!!"
      else
           allocate(f(m(1)))
           do i=1,m(1)
           print*,"Please enter them one by one.(",m(1)-i+1,"number(s) left)"
           read*, f(i)
         end do
      print*,'The mean of the numbers you input is ',meanArray(f,m(1),SD)
      print*,'The Standard Deviation of the numbers you input is ',SD
      end if
end program exercise2_2
!------------------------
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
