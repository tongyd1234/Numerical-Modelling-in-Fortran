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
end module mean_std
!--------------------------------------------------
program exercise3_2
        use mean_std
        implicit none
        integer::n=0, i
        real::b
        real,allocatable::a(:)
        open(1,file='exercise2.dat',status='old')
        do
            read(1,*,iostat=i) b
            if(i<0) exit
            if(i/=0) stop 'Error reading data'
            n=n+1
        end do
        print*,'Found',n,'value(s)'
        allocate (a(n))
        rewind(1)
        do i=1,n
           read(1,*) a(i)
        end do
        print*, a
        print*,'The mean of the numbers you input is ', meanArray(a,n,SD)
        print*,'The Standard Deviation of the numbers you input is ',SD
end program exercise3_2
