module oceanography_tools
  use linalg_tools
  use constants

contains
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

    kappa = kappa_f_coefs(1)*n_freq**(-kappa_f_coefs(2))
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
      n_freq2(i) = -g*drho_dz(i)/rho(i)
      n_freq(i) = sqrt(abs(n_freq2(i)))
    end do
  end subroutine calc_buoyancy_frequency
end module oceanography_tools
