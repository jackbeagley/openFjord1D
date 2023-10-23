#define N_DERIVED 3
#define I_KAPPA 1
#define I_N_FREQ 2
#define I_RHO 3

module indices
  ! Export variables for accessing the derived variables
  integer, parameter :: n_derived = N_DERIVED
  integer, parameter :: kappa = I_KAPPA
  integer, parameter :: n_freq = I_N_FREQ
  integer, parameter :: rho = I_RHO
end module indices

module diffusivity_modes
  ! Export variables for different diffusivity modes
  integer, parameter :: constant = 0
  integer, parameter :: buoyancy = 1
end module diffusivity_modes

module computation_modes
  ! Export variables for different computation modes
  integer, parameter :: normal = 0
  integer, parameter :: closed_top = 1
end module

module basin
use gsw_mod_toolbox
use linalg_debug
use linalg_tools
use oceanography_tools
use diffusivity_modes, only: &
  diff_constant => constant, &
  diff_buoyancy => buoyancy
use computation_modes, only: &
  comp_normal => normal, &
  comp_closed_top => closed_top

implicit none
  ! Set latitude to Doubtful Sound
  double precision, parameter :: lat = -45.5

contains
  subroutine compute(z, a, sa, tis, derived, t, sa_bc, tis_bc, computation_mode, diffusivity_mode, kappa_coefs, n_z, n_t)
    ! Compute the salinity and temperature evolution along with other parameters in the
    ! water column as a consequence of the diffusivity and boundary conditions.
    ! Input:
    !   z - depths to compute at
    !   a - area at each depth
    !   sa - absolute salinity at each depth and time (first column is initial condition)
    !   tis - in-situ temperature at each depth and time (first column is initial condition)
    !   derived - derived parameters at each depth and time
    !   t - times to compute at
    !   sa_bc - boundary condition for absolute salinity
    !   tis_bc - boundary condition for in-situ temperature
    !   computation_mode - mode to compute in (0 = normal, 1 = closed top)
    !   diffusivity_mode - mode to compute diffusivity in (0 = constant, 1 = buoyancy varying)
    !   kappa_coefs - k=k_coefs(1)*N**k_coefs(2) OR k=k_coefs(1) for constant diffusivity 
    !   n_z - number of depths
    !   n_t - number of times

    implicit none
    integer, intent(in) :: computation_mode, diffusivity_mode, n_z, n_t
    double precision, dimension(n_z), intent(inout) :: z, a
    double precision, dimension(n_z, n_t), intent(inout) :: sa, tis
    double precision, dimension(n_z, n_t, N_DERIVED), intent(inout) :: derived
    double precision, dimension(n_t), intent(in) :: t, sa_bc, tis_bc
    integer :: i, j

    double precision, dimension(n_z) :: p, n_freq, n_freq2, dadz, dkappadz

    double precision :: dz, dt

    ! double precision, dimension(n_z, n_z) :: A_mat
    double precision, dimension(3+1, n_z) :: A_mat_band
    double precision, dimension(n_z, 2) :: bx_vec

    double precision, dimension(2) :: kappa_coefs
    
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
        derived(j, i, I_RHO) = gsw_rho_t_exact(sa(j, i), tis(j, i), p(j))
      end do

      ! Calculate the Brunt-Vaisala (buoyancy) frequency
      call calc_buoyancy_frequency(derived(:, i, I_RHO), z, n_freq, n_freq2, n_z)
  
      ! Calculate the diffusivity
      do j = 1, n_z
        derived(j, i, I_KAPPA) = calc_diffusivity(n_freq(j), kappa_coefs)
        ! kappa(j) = 5d-3
      end do
  
      ! Calculate the derivative of diffusivity with respect to depth
      call dydx_array(z, derived(:, i, I_KAPPA), dkappadz, n_z)

      ! Calculate the A matrix
      call compute_A(A_mat_band, derived(:, i, I_KAPPA), dkappadz, a, dadz, dz, dt, n_z)

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
    call init_matrix(A_mat_band, n_z, 3+1)

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
end module basin
