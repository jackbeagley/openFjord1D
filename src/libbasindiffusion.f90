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
  subroutine compute(z, a, sa, tis, rho, kappa, t, sa_bc, tis_bc, kappa_f_coefs, n_z, n_t)
    implicit none
    integer, intent(in) :: n_z, n_t
    double precision, dimension(n_z), intent(inout) :: z, a
    double precision, dimension(n_z, n_t), intent(inout) :: sa, tis, rho, kappa
    double precision, dimension(n_t), intent(in) :: t, sa_bc, tis_bc
    integer :: i, j

    double precision, dimension(n_z) :: p, n_freq, n_freq2, dadz, dkappadz

    double precision :: dz, dt

    ! double precision, dimension(n_z, n_z) :: A_mat
    double precision, dimension(3+1, n_z) :: A_mat_band
    double precision, dimension(n_z, 2) :: bx_vec

    double precision, dimension(2) :: kappa_f_coefs
    
    ! Output variables for DGESV
    integer, dimension(n_z) :: ipiv
    integer :: info

    dz = (z(n_z) - z(1))/(n_z - 1)
    dt = (t(n_t) - t(1))/(n_t - 1)

    ! Calculate the derivative of area with respect to depth
    call dydx_array(z, a, dadz, n_z)

    do i = 1, (n_t-1)
      do j = 1, n_z
        ! Pressure from depth
        p(j) = gsw_p_from_z(z(j), lat)

        ! Density from absolute salinity, in-situ temperature and pressure
        rho(j, i) = gsw_rho_t_exact(sa(j, i), tis(j, i), p(j))
      end do

      ! Calculate the Brunt-Vaisala (buoyancy) frequency
      call calc_buoyancy_frequency(rho(:, i), z, n_freq, n_freq2, n_z)
  
      ! Calculate the diffusivity
      do j = 1, n_z
        kappa(j, i) = calc_diffusivity(n_freq(j), kappa_f_coefs)
        ! kappa(j) = 5d-3
      end do
  
      ! Calculate the derivative of diffusivity with respect to depth
      call dydx_array(z, kappa(:, i), dkappadz, n_z)

      ! Calculate the A matrix
      call compute_A(A_mat_band, kappa(:, i), dkappadz, a, dadz, dz, dt, n_z)

      ! Calculate the RHS vector
      ! Dirichlet boundary condition at the top
      bx_vec(1, 1) = sa_bc(i+1)
      bx_vec(1, 2) = tis_bc(i+1)
      ! Neuamnn boundary condition at the bottom
      bx_vec(n_z, 1) = 0d0
      bx_vec(n_z, 2) = 0d0

      ! Set the interior points to the current salinity values
      do j = 2, (n_z-1)
        bx_vec(j, 1) = sa(j, i)
        bx_vec(j, 2) = tis(j, i)
      end do

      ! call dgesv(n_z, 2, A_mat, n_z, ipiv, bx_vec, n_z, info)
      call dgbsv(n_z, 1, 1, 2, A_mat_band, 4, ipiv, bx_vec, n_z, info)

      do j = 1, n_z
        sa(j, i+1) = bx_vec(j, 1)
        tis(j, i+1) = bx_vec(j, 2)
      end do
    end do
  end subroutine compute

  subroutine compute_A(A_mat_band, kappa, dkappadz, a, dadz, dz, dt, n_z)
    implicit none
    integer, intent(in) :: n_z
    double precision, dimension(3+1, n_z), intent(out) :: A_mat_band
    double precision, dimension(n_z), intent(in) :: kappa, dkappadz, a, dadz
    double precision, intent(in) :: dz, dt
    double precision :: alpha, beta
    integer :: i

    ! call init_A(A_mat, n_z)
    call init_A_band(A_mat_band, n_z, 3+1)

    ! Set the top of the water column as a Dirichlet boundary condition
    A_mat_band(2+1, 1) = 1.d0

    do i = 2, (n_z-1)
      ! Calculate the FDM coefficients
      alpha = (dadz(i) * kappa(i) / a(i) + dkappadz(i)) / (2 * dz)
      beta = kappa(i) / (dz**2)
    
      ! Fill out the A matrix
      A_mat_band(3+1, i-1) = dt * (alpha - beta)
      A_mat_band(2+1, i) = 1.d0 + 2.d0 * dt * beta
      A_mat_band(1+1, i+1) = -dt * (alpha + beta)
    end do

    ! Set the bottom of the water column as a Neumann boundary condition
    A_mat_band(3+1, n_z-1) = -1.d0
    A_mat_band(2+1, n_z) = 1.d0
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

  subroutine init_A_band(A_mat, n_z, n_band) 
    implicit none
    integer, intent(in) :: n_z, n_band
    double precision, dimension(n_band, n_z), intent(out) :: A_mat
    integer :: i, j

    do i = 1, n_band
      do j = 1, n_z
        A_mat(i, j) = 0.d0
      end do
    end do
  end subroutine init_A_band

  function calc_diffusivity(n_freq, kappa_f_coefs) result(kappa)
    ! Estimate diffusivity, as per Gargett et al. (1984), Young et al. (1988)
    implicit none
    double precision, intent(in) :: n_freq
    double precision, dimension(2), intent(in) :: kappa_f_coefs
    double precision :: kappa

    ! check if n_freq is nan
    if (isnan(n_freq)) then
      kappa = 0d0
      return
    end if

    kappa = kappa_f_coefs(1) * n_freq**(-kappa_f_coefs(2))
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
