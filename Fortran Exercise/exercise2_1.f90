program exercise2_1
      implicit none
!Define a a string of length 15 characters
      character(len=15)::a
!Declare an integer parameter = 5
      integer, parameter::b=5
!Declare a 1-dimensional array with indices running from -1 to +10
      real::c(-1:10)
!Declare a 4-dimensional allocatable array
      integer, allocatable::d( , , , )
!Convert a real number to the nearest integer
!And Calculate the remainder after a is divided by b
      real::e
      integer ::f,g
      integer ::i=2,j=1
      integer:: sum=0
      real::h(1:100)
      print*,'Input the number you want to convert to the nearest integer.'
      read*,e
      print*,'Input the number you want to get the remainder divided by b.'
      read*,f
      print*,'The nearest integer of ',e,'is ', nint(e)
      print*,'The remainder of ',f,'divided by b is', mod(f,b)
      print*,c
!Add up all even numbers between 12 and 124
      do while (i<=124)
         sum=sum+i
         i=i+2
      end do
      print*,sum
!–  Test each element of array a(1:100) starting from 1; if the element is positive print a message to the screen and leave the loop
      do
         if (h(j)<=0) then
            j=j+1
         else
            print*,'It is a positive number, now EXIT!'
            exit
         end if
      end do
      print*,j
end program exercise2_1
      
