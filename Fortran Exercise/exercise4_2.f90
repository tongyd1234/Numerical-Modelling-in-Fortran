module delta_square
  implicit none
  contains
     function derivative (a,nx,ny,h)
         implicit none
         double precision,dimension(:,:),intent(in)::a
         integer,intent(in)::nx,ny
         double precision,intent(in)::h
         real,dimension(nx,ny)::derivative
         integer::i,j
         do i=2,nx-1
            do j=2,ny-1
              derivative(i,j)=(a(i+1,j)-2*a(i,j)+a(i-1,j))/(h*h)+(a(i,j+1)-2*a(i,j)+a(i,j-1))/(h*h)
            end do
         end do
         do j=1,ny,ny-1
            do i=1,nx
               derivative(i,j)=0
            end do
         end do
         do i=1,nx,nx-1
            do j=1,ny
               derivative(i,j)=0
            end do
         end do
     end function derivative
end module delta_square
