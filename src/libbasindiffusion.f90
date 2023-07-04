module basindiffusion
use gsw_mod_toolbox
use matrix_debug

implicit none

  ! Public interface of the library
  public :: compute

  ! Declare a constant for gravity
  double precision, parameter :: g = 9.80665

  ! Set latitude to Doubtful Sound
  double precision, parameter :: lat = -45.5

contains
  subroutine compute(z, a, sa, tis, t, sa_bc, tis_bc, n_z, n_t)
    implicit none
    integer, intent(in) :: n_z, n_t
    double precision, dimension(n_z), intent(inout) :: z, a
    double precision, dimension(n_z, n_t), intent(inout) :: sa, tis
    double precision, dimension(n_t), intent(in) :: t, sa_bc, tis_bc
    integer :: i, j, k

    double precision, dimension(n_z) :: rho, p, ct, n_freq, n_freq2, kappa, dadz, dkappadz

    double precision :: dz, dt

    double precision, dimension(n_z, n_z) :: A_mat
    double precision, dimension(n_z) :: b_vec, x_vec

    double precision, dimension(2) :: A_coefs_z

    ! Output variables for DGESV
    integer, dimension(n_z) :: ipiv
    integer :: info


    dz = (z(n_z) - z(1))/(n_z - 1)
    dt = (t(n_t) - t(1))/(n_t - 1)

    ! Calculate the derivative of area with respect to depth
    call dydx_array(z, a, dadz, n_z)

    ! ##############################
    ! #  Initialise the A matrix   #
    ! ##############################
    do i = 1, n_z
      do j = 1, n_z
        A_mat(i, j) = 0.d0
      end do
    end do

    do i = 1, (n_t-1)
      do j = 1, n_z
        ! Pressure from depth
        p(j) = gsw_p_from_z(z(j), lat)
        ! Conservative temperature from in-situ temperature
        ct(j) = gsw_ct_from_t(sa(j, i), tis(j, i), p(j))
        ! Density from absolute salinity, conservative temperature and pressure
        rho(j) = gsw_rho(sa(j, i), ct(j), p(j))
      end do

      ! Calculate the Brunt-Vaisala (buoyancy) frequency
      call calc_buoyancy_frequency(rho, z, n_freq, n_freq2, n_z)
  
      ! Calculate the diffusivity
      do j = 1, n_z
        kappa(j) = calc_diffusivity(n_freq(j), 2.d-3, 1.6d0)
        ! kappa(j) = 5d-3
      end do

      call print_vector(n_freq, n_z, "n_freq")
      call print_vector(kappa, n_z, "kappa")
  
      ! Calculate the derivative of diffusivity with respect to depth
      call dydx_array(z, kappa, dkappadz, n_z)

      ! Calculate the A matrix
      call compute_A(A_mat, kappa, dkappadz, a, dadz, dz, dt, n_z)

      ! Calculate the RHS vector
      ! Dirichlet boundary condition at the top
      b_vec(1) = sa_bc(i+1)
      ! Neuamnn boundary condition at the bottom
      b_vec(n_z) = 0d0

      ! Set the interior points to the current salinity values
      do j = 2, (n_z-1)
        b_vec(j) = sa(j, i)
      end do

      call dgesv(n_z, 1, A_mat, n_z, ipiv, b_vec, n_z, info)

      do j = 1, n_z
        sa(j, i+1) = b_vec(j)
      end do
    end do
  end subroutine compute

  subroutine compute_A(A_mat, kappa, dkappadz, a, dadz, dz, dt, n_z)
    implicit none
    integer, intent(in) :: n_z
    double precision, dimension(n_z, n_z), intent(out) :: A_mat
    double precision, dimension(n_z), intent(in) :: kappa, dkappadz, a, dadz
    double precision, intent(in) :: dz, dt
    double precision :: alpha, beta
    integer :: i

    call init_A(A_mat, n_z)

    ! Set the top of the water column as a Dirichlet boundary condition
    A_mat(1, 1) = 1.d0

    do i = 2, (n_z-1)
      ! Calculate the FDM coefficients
      alpha = (dadz(i) * kappa(i) / a(i) + dkappadz(i)) / (2 * dz)
      beta = kappa(i) / (dz**2)
    
      ! Fill out the A matrix
      A_mat(i, i-1) = dt * (alpha - beta)
      A_mat(i, i) = 1.d0 + 2.d0 * dt * beta
      A_mat(i, i+1) = -dt * (alpha + beta)
    end do

    ! Set the bottom of the water column as a Neumann boundary condition
    A_mat(n_z, n_z-1) = -1.d0
    A_mat(n_z, n_z) = 1.d0
  end subroutine compute_A

  subroutine init_A(A_mat, n_z) 
    implicit none
    integer, intent(in) :: n_z
    double precision, dimension(n_z, n_z), intent(out) :: A_mat
    integer :: i, j

    do i = 1, n_z
      do j = 1, n_z
        A_mat(i, j) = 0.d0
      end do
    end do
  end subroutine init_A

  function calc_diffusivity(n_freq, a0, q) result(kappa)
    ! Estimate diffusivity, as per Gargett et al. (1984), Young et al. (1988)
    implicit none
    double precision, intent(in) :: n_freq, a0, q
    double precision :: kappa

    ! check if n_freq is nan
    if (isnan(n_freq)) then
      kappa = 0d0
      return
    end if

    kappa = a0 * exp(-q * n_freq)
  end function calc_diffusivity

  subroutine calc_buoyancy_frequency(rho, z, n_freq, n_freq2, n_z)
    implicit none
    integer, intent(in) :: n_z
    double precision, dimension(n_z), intent(in) :: rho, z
    double precision, dimension(n_z), intent(out) :: n_freq, n_freq2
    double precision, dimension(n_z) :: drho_dz

    integer :: i

    ! Calculate the derivative of density with respect to depth
    call dydx_array(z, rho, drho_dz, n_z)

    ! Calculate the Brunt-Vaisala (buoyancy) frequency
    do i = 1, n_z
      n_freq2(i) = -g * drho_dz(i) / rho(i)
      n_freq(i) = sqrt(abs(n_freq2(i)))
    end do
  end subroutine calc_buoyancy_frequency

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

end module basindiffusion
