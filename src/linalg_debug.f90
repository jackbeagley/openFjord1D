module linalg_debug
  implicit none

  public :: print_matrix, print_vector
  character(len=*), parameter :: fmt = '(ES10.3, " ")'

contains
  subroutine print_matrix(A, n, m, title)
    implicit none
    integer, intent(in) :: n, m
    double precision, dimension(n, m), intent(in) :: A
    integer :: i, j
    character(len=*), intent(in) :: title

    write (*, *) title

    do i = 1, n
      do j = 1, m
        write (*, fmt, advance='no') A(i, j)
      end do
      write (*, *) ""
    end do
  end subroutine print_matrix

  subroutine print_vector(b, n, title)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: b
    integer :: i
    character(len=*), intent(in) :: title

    write (*, *) title

    do i = 1, n
      write (*, fmt, advance='no') b(i)
    end do
    write (*, *) ""
  end subroutine print_vector
end module linalg_debug
