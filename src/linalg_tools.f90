module linalg_tools
contains
    function dydx_centre(x, y, i) result(dydx)
    implicit none
    integer, intent(in) :: i
    double precision, dimension(i), intent(in) :: x, y
    double precision :: dydx

    dydx = (y(i+1) - y(i-1)) / (x(i+1) - x(i-1))
    end function dydx_centre

    function dydx_forward(x, y, i) result(dydx)
    implicit none
    integer, intent(in) :: i
    double precision, dimension(i), intent(in) :: x, y
    double precision :: dydx

    dydx = (y(i+1) - y(i)) / (x(i+1) - x(i))
    end function dydx_forward

    function dydx_backward(x, y, i) result(dydx)
    implicit none
    integer, intent(in) :: i
    double precision, dimension(i), intent(in) :: x, y
    double precision :: dydx

    dydx = (y(i) - y(i-1)) / (x(i) - x(i-1))
    end function dydx_backward

    subroutine dydx_array(x, y, dydx, n)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x, y
    double precision, dimension(n), intent(out) :: dydx

    integer :: i

    dydx(1) = dydx_forward(x, y, 1)
    dydx(n) = dydx_backward(x, y, n)

    do i = 2, n - 1
        dydx(i) = dydx_centre(x, y, i)
    end do
    end subroutine dydx_array

    subroutine init_matrix(matrix, m, n)
      implicit none
      integer, intent(in) :: m, n
      double precision, dimension(m, n), intent(inout) :: matrix
      integer :: i, j
  
      do i = 1, m
        do j = 1, n
          matrix(i, j) = 0.d0
        end do
      end do
    end subroutine init_matrix
end module linalg_tools