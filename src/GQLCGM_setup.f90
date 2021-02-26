!--------------------------------------------------------
! GQLCGM_setup.f90
!
! PURPOSE:
!  Setup intial time=0 values for the GQLCGM_transport callsC
!--------------------------------------------------------

subroutine GQLCGM_setup

  use mpi
  use math_constants
  use gyro_profile_exp
  use gyro_globals
  use GQLCGM_globals

  !----------------------------------------
  implicit none
  !----------------------------------------

  real :: kt
  real , dimension(2) :: mass_EP
  real , dimension(2) :: Z_EP
  real , dimension(2) :: Amass_EP

  call send_line('**[ GQLCGM_setup starting]')

  n_grid_EP = 51 ! set for n_grid_exp = 51
 
  do i_EP = 1,n_grid_EP
   r_EP(i_EP) = rmin_exp(i_EP)            ![m]
   rmaj_EP(i_EP) = rmaj_exp(i_EP)         ![m]
   kappa_EP(i_EP) = kappa_exp(i_EP)       

   den_sd_EP(1,i_EP) = den_exp(2,i_EP)     ![10**19/m**3]
   den_sd_EP(2,i_EP) = den_exp(3,i_EP)     ![10**19/m**3]

   temp_sd_EP(1,i_EP) = tem_exp(2,i_EP)    ![keV]
   temp_sd_EP(2,i_EP) = tem_exp(3,i_EP)    ![keV]
  enddo
  call send_line('**[ GQLCGM_setup after 1st i_EP loop')

  do i_EP = 1,n_grid_EP
      Vprime_EP(i_EP) = (2.*pi*rmaj_EP(i_EP))*(2.*pi*r_EP(i_EP)*kappa_EP(i_EP)) !simple ellipse geo
  enddo
  call send_line('**[ GQLCGM_setup after 2nd i_EP loop')

  !starting EP densities
    den_prev_EP(:,:)=0.97*den_sd_EP(:,:)

! prep for tau_sd_EP
! z_vec(1:10) ions z_vec(0) = -1. (electrons)
! mu_vec(1:10) ions  mu_vec(0) = sqrt(mass(1)/mass_electron)
! mu_i = sqrt(m(1)/m(i))
! mu_vec(1:10)   mu_vec(0) = gyro_mu_electron_in
! real, parameter :: kg_proton = 1.6726e-27   (kg)

! gyro_alloc_big.f90:     allocate(gbflux_i(n_kinetic,n_field,p_moment,n_x))
! gyro_globals.f90:  real, dimension(:,:,:,:), allocatable :: gbflux_i
! gyro_gbflux.f90
 !  if (nonlinear_flag == 0 .and. n_1(in_1) /= 0) then
 !     if (lindiff_method == 2) then
 !        ! [Gamma_a, Q_a, Pi_a, S_a] / Q_i normalization for linear diffusivities.
 !        gbflux_norm = gbflux(1,1,2)
 !         if(i_do_GQLCGM == 1) gbflux_norm =1.0  !GQLCGM 01.22.21
 !        gbflux_i(:,:,:,:)         = gbflux_i(:,:,:,:)/gbflux_norm
 !     endif
 !  endif


   flux_gbnorm_EP = den_norm*csda_norm*a_meters*rhos_norm**2                !   [10**19/m**3]*[m/sec]

  if(i_proc == 0) then
   ! kT in MJ (note the conversion 1.6022e-22 MJ/keV)
   kt = 1.6022e-22*tem_norm

   print *, 2.0*kg_proton,'m_ref (kg)'
   print *, b_unit_norm,'b_unit (Tesla)'
   print *, a_meters,'a (m)'
   print *,  csda_norm,'csD/a (1/s)'
   print *,  csda_norm*a_meters,'csD (m/s)'
   print *,  tem_norm,'Te (keV)'
   print *, den_norm,'ne (10^19/m^3)'
   print *, rhos_norm*a_meters,'rho_sD (m)'
   print *,  csda_norm*(rhos_norm*a_meters)**2,&
       'chi_gBD (m^2/s)'
   print *, 1e19*den_norm*(csda_norm*a_meters)*rhos_norm**2/0.624e22,&
       'Gamma_gBD (0.624e22/m^2/s) = (MW/keV/m^2)'
   print *,  1e19*den_norm*(csda_norm*a_meters)*kt*rhos_norm**2,&
       'Q_gBD (MJ/m^2/s) = (MW/m^2)'
   print *,  1e19*den_norm*a_meters*kt*rhos_norm**2*1e6,&
       'Pi_gBD (J/m^2) = (Nm/m^2)'
   print *,  1e19*den_norm*csda_norm*kt*rhos_norm**2,&
       'S_gBD (MJ/m^3/s) = (MW/m^3)'

  endif

  !from ALPHA code
  !       ln_lambda=17.
  !  tau_ee=1.088E-3*(T_e_rho(i))**1.5/n_e_rho(i)/ln_lambda ! in sec
  !  tau_s = 1836.*tau_ee &                  !!!for alphas
  !     *(M_DT/4.0)/(1.0/2.0)**2             !!!M_DT = atomic number of EP proton with the chargeZ=+1.  and Amass=1.


 Amass_EP(1) = 2.0  !HARDWIRED 
 Amass_EP(2) = 4.0  !HARDWIRED

  mass_EP(:) = Amass_EP(:)*kg_proton     !kg 

   z_EP(1) = z_vec(2)
   z_EP(2) = z_vec(3)

!  do i_sEP =1,2
  do i_sEP =1,i_do_GQLCGM
   do i_EP = 1,n_grid_EP
     tau_sd_EP(i_sEP,i_EP) = 1836.*1.088E-3*(tem_exp(n_spec,i_EP))**1.5/den_exp(n_spec,i_EP)/(17.)*(Amass_EP(i_sEP)/4.0)/(Z_EP(i_sEP)/2.)**2
   enddo
  enddo

  if(i_proc == 0) then
   print *, 'HARDWIRED Amass_EP(:)=', Amass_EP(:)         
   print *, 'Z_EP(:)=',Z_EP(:)
   print *, 'tau_sd_EP(:,25)=',tau_sd_EP(:,25),'sec'
   print *, 'flux_gbnorm_EP =',flux_gbnorm_EP,'[10**19/m**3]*[m/sec]'
  endif

!  D_EP_bkg(:,:) = 0.2 ![m**2/sec]   !replace with Angioni model  or TGLF-EP likely read from external file
!  D_EP_bkg(:,:) = 0.05  !OK BC=0. 10 BC=1.0 11
!  D_EP_bkg(:,:) = 0.04   !
!  D_EP_bkg(:,:) = 0.02
  D_EP_bkg(:,:) = 0.015   !2.18.21

  do i_sEP = 1,i_do_GQLCGM  !GQLCGM 218.21
   den_sd_EP_s(i_sEP,:) = den_s(1+i_do_GQLCGM,:)*den_norm  !2.22.21 den_norm fixed
   den_EP_s(i_sEP,:) = den_s(1+i_do_GQLCGM,:)*den_norm
   dlnndr_EP_s(i_sEP,:)=dlnndr_s(1+i_do_GQLCGM,:)
  enddo

  i_update_gyro_GQLCGM = 0

  call send_line('**[ GQLCGM_setup done]')

   if (debug_flag == 1 .and. i_proc == 0) then
     print '(t2,2(a,i5,3x))','-> step =',step,'data_step =',data_step
     print *,'**[ GQLCGM_setup done]'
  endif

end subroutine GQLCGM_setup
