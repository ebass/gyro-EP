!-----------------------------------------------------
! GLLCGM_globals.f90
!
! PURPOSE:
!  GQLCGM variables and arrays
!-----------------------------------------------------

module GQLCGM_globals

  integer :: i_do_GQLCGM !0 no ! 1 yes
  integer :: i_update_gyro_GQLCGM !2.22.21pm
  integer :: n_grid_EP ! set to 51 in GQLCGM_setup
  integer :: i_EP   !1,2...n_grid_EP
  integer :: i_sEP  !1,2   index over number of EP species no more than 2

  !caution EXPRO_n_exp = n_grid_exp = 51 assumed

  !time changing
  real :: time_EP                       !   [sec]
  real :: time_prev_EP                  !   [sec]

  real, dimension(2,51) :: den_EP       !   [10**19/m**3]
  real, dimension(2,51) :: den_prev_EP  !   [10**19/m**3]

  real, dimension(2,51) :: flow_EP      !   [10**19/sec]

  !not changed in time
  real, dimension(2,51) :: den_sd_EP    !   [10**19/m**3]
  real, dimension(2,51) :: temp_sd_EP   !   [keV]

  real, dimension(2,51) :: tau_sd_EP    !   [1/sec]
  real :: flux_gbnorm_EP                !   [10**19/m**3]*[m/sec]  from gb 
  real, dimension(2,51) :: flux_EP      !   [10**19/m**3]*[m/sec]

  real, dimension(2,51) :: D_EP_bkg     !   [m**2/sec]
  real, dimension(2) :: error_EP        !   RMS difference den_EP with den_prev_EP

  real, dimension(51) :: rho_EP         !   [m]
  real, dimension(51) :: r_EP           !   [m]
  real, dimension(51) :: rmaj_EP        !   [m]
  real, dimension(51) :: kappa_EP       !
  real, dimension(51) :: Vprime_EP      !   [m**2]


end module GQLCGM_globals
