! program myprog
!     x = -3.0
!     y = 2.0* x + 3.0
!     print * , y
! end program

! subroutine basin_diffusion(z, SA, CT, t, SA_bc, CT_bc)
!     intent (in) :: z, t, SA_bc, CT_bc
!     intent (inout) :: SA, CT

! end subroutine  basin_diffusion

module my_fortran_library
    implicit none
  
    ! Public interface of the library
    public :: initialize, compute
  
    ! Library functions
    integer :: initialize_called = 0
  
  contains
  
    ! Function to initialize the library
    subroutine initialize()
      initialize_called = 1
      ! Additional initialization code goes here
    end subroutine initialize
  
    ! Function to compute something
    function compute(a, b) result(result)
      integer, intent(in) :: a, b
      integer :: result
  
      if (initialize_called == 0) then
        print *, "Error: Library not initialized!"
        result = -999
      else
        result = a + b
      end if
    end function compute
  
  end module my_fortran_library    

