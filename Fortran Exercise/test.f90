program test
   implicit none
      real::L,total_time,time,kappa,a
      integer::N,M
      integer::i,j
      namelist/inputs/ L,total_time,kappa,a,N,M,i,j
      open(1,file='parameters.txt',status='old')
      read(1,inputs)
      close(1)
      print*,L
end program test
