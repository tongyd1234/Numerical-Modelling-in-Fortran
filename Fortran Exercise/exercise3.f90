program exercise3
       implicit none
       real::a
       integer:: i, factorial=1
       print*,"Please enter a positive integer: "
       read*,a
       if(a>0) then
           if(a-floor(a)==0) then
             do i=1, floor(a)
                 factorial=i*factorial
             end do
             print*,"The factorial of ",int(a)," is ",factorial," . "
           else
             print*,"The number is not a positive integer."
           end if
       else
             print*,"The number you entered is not even positive!"
       end if

end program exercise3
                
