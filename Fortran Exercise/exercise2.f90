program exercise2
   implicit none
      real::a,b,c
      real::Amean, Hmean, Gmean
      print*,"Please enter 3 real numbers:"
      read*,a,b,c
      Amean=(a+b+c)/3.0
      Hmean=3.0/(1.0/a+1.0/b+1.0/c)
      Gmean=(a*b*c) **(1.0/3.0)
      print*,"The Arithmetic Mean of the 3 numbers is: ", Amean
      print*,"The Harmonic Mean of the 3 numbers is: ", Hmean
      print*,"The Geometric Mean of the 3 numbers is: ", Gmean
end program exercise2
