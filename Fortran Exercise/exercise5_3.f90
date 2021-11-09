Module advdifequation
implicit none
      type parameters
           integer nx,ny
           real:: xsize,ysize 
           real:: B  ! flowstrength
           real:: kappa  ! diffusivity
           real:: a_dif,a_adv
           real:: total_time
           real:: time
           real:: dx,dy
           real:: dt
           real,allocatable,dimension(:,:) :: T
           real,allocatable,dimension(:,:,:) :: V
           real,allocatable,dimension(:,:) :: S
      end type parameters
contains
