!--------------------------------------------------------------
! EXPRO_read_driver.f90
!
! PURPOSE:
!  Low-level read control for EXPRO.
!--------------------------------------------------------------

subroutine EXPRO_read_driver

  use EXPRO_globals
  use EXPRO_interface

!  use Alpha_use_input
!  use Alpha_use_output

  implicit none

  integer, parameter :: io=2
  integer :: i
  integer :: ierr

  integer :: EXPRO_n_shot
  integer :: EXPRO_nion_file
  real :: EXPRO_ip_exp
  real :: EXPRO_rvbv

  real, dimension(:), allocatable :: dummy

  !--------------------------------------------------------------
  ! READ generated (stripped) version of input.profiles:
  !
  open(unit=1,file='EXPRO_debug.dat',status='replace')

  open(unit=io,&
       file=trim(path)//'input.profiles.gen',&
       status='old',&
       iostat=ierr)

  read(io,*) EXPRO_ncol
  write(1,*) EXPRO_ncol

  read(io,*) EXPRO_nblock
  write(1,*) EXPRO_nblock

  read(io,*) EXPRO_n_shot
  write(1,*) EXPRO_n_shot

  read(io,*) EXPRO_nion_file
  write(1,*) EXPRO_nion_file

  read(io,*) EXPRO_n_exp
  write(1,*) EXPRO_n_exp

  read(io,*) EXPRO_b_ref
  write(1,*) EXPRO_b_ref

  read(io,*) EXPRO_ip_exp
  write(1,*) EXPRO_ip_exp

  read(io,*) EXPRO_rvbv
  write(1,*) EXPRO_rvbv

  read(io,*) EXPRO_arho
  write(1,*) EXPRO_arho

  allocate(dummy(EXPRO_n_exp))

  ! 1-5
  read(io,*) EXPRO_rho(:)
  write(1,*) EXPRO_rho(:)
  read(io,*) EXPRO_rmin(:)
  write(1,*) EXPRO_rmin(:)
  read(io,*) EXPRO_rmaj(:)
  write(1,*) EXPRO_rmaj(:)
  read(io,*) EXPRO_q(:)     ! |q|
  write(1,*) EXPRO_q(:)
  read(io,*) EXPRO_kappa(:)
  write(1,*) EXPRO_kappa(:)

  ! 6-10
  read(io,*) EXPRO_delta(:)
  write(1,*) EXPRO_delta(:)
  read(io,*) EXPRO_te(:)
  write(1,*) EXPRO_te(:)
  read(io,*) EXPRO_ne(:)
  write(1,*) EXPRO_ne(:)
  read(io,*) EXPRO_z_eff(:)
  write(1,*) EXPRO_z_eff(:)
  read(io,*) EXPRO_w0(:)     ! Note that EXPRO_w0 has GYRO/NEO sign
  write(1,*) EXPRO_w0(:)

  ! 11-15
  read(io,*) EXPRO_flow_mom(:)
  write(1,*) EXPRO_flow_mom(:)
  read(io,*) EXPRO_pow_e(:)
  write(1,*) EXPRO_pow_e(:)
  read(io,*) EXPRO_pow_i(:)
  write(1,*) EXPRO_pow_i(:)
  read(io,*) EXPRO_pow_ei(:)
  write(1,*) EXPRO_pow_ei(:)
  read(io,*) EXPRO_zeta(:)
  write(1,*) EXPRO_zeta(:)

  ! 16-20
  read(io,*) EXPRO_flow_beam(:)
  write(1,*) EXPRO_flow_beam(:)
  read(io,*) EXPRO_flow_wall(:)
  write(1,*) EXPRO_flow_wall(:)
  read(io,*) EXPRO_zmag(:)
  write(1,*) EXPRO_zmag(:)
  read(io,*) EXPRO_ptot(:)
  write(1,*) EXPRO_ptot(:)
  read(io,*) EXPRO_poloidalfluxover2pi(:)
  write(1,*) EXPRO_poloidalfluxover2pi(:)

  ! 21-25
  read(io,*) EXPRO_ni(1,:)
  write(1,*) EXPRO_ni(1,:)
  read(io,*) EXPRO_ni(2,:)
  write(1,*) EXPRO_ni(2,:)
  read(io,*) EXPRO_ni(3,:)
  write(1,*) EXPRO_ni(3,:)
  read(io,*) EXPRO_ni(4,:)
  write(1,*) EXPRO_ni(4,:)
  read(io,*) EXPRO_ni(5,:)
  write(1,*) EXPRO_ni(5,:)

  ! 26-30
  read(io,*) EXPRO_ti(1,:)
  write(1,*) EXPRO_ti(1,:)
  read(io,*) EXPRO_ti(2,:)
  write(1,*) EXPRO_ti(2,:)
  read(io,*) EXPRO_ti(3,:)
  write(1,*) EXPRO_ti(3,:)
  read(io,*) EXPRO_ti(4,:)
  write(1,*) EXPRO_ti(4,:)
  read(io,*) EXPRO_ti(5,:)
  write(1,*) EXPRO_ti(5,:)

  ! 31-35
  read(io,*) EXPRO_vtor(1,:)
  write(1,*) EXPRO_vtor(1,:)
  read(io,*) EXPRO_vtor(2,:)
  write(1,*) EXPRO_vtor(2,:)
  read(io,*) EXPRO_vtor(3,:)
  write(1,*) EXPRO_vtor(3,:)
  read(io,*) EXPRO_vtor(4,:)
  write(1,*) EXPRO_vtor(4,:)
  read(io,*) EXPRO_vtor(5,:)
  write(1,*) EXPRO_vtor(5,:)

  ! 36-40
  read(io,*) EXPRO_vpol(1,:)
  write(1,*) EXPRO_vpol(1,:)
  read(io,*) EXPRO_vpol(2,:)
  write(1,*) EXPRO_vpol(2,:)
  read(io,*) EXPRO_vpol(3,:)
  write(1,*) EXPRO_vpol(3,:)
  read(io,*) EXPRO_vpol(4,:)
  write(1,*) EXPRO_vpol(4,:)
  read(io,*) EXPRO_vpol(5,:)
  write(1,*) EXPRO_vpol(5,:)

  close(io)

  write (1,*) 'Finished read of input.profiles.gen.'

  close(1)

  !--------------------------------------------------------------

  ! Change signs according to orientation control parameters

  EXPRO_b_ref = EXPRO_ctrl_signb*abs(EXPRO_b_ref)
  EXPRO_q(:)  = EXPRO_ctrl_signq*abs(EXPRO_q(:))

  !--------------------------------------------------------------
  ! READ general shape coefficients if they exist: 
  !
  if (EXPRO_nfourier > 0) then

     open(unit=io,&
          file=trim(path)//'input.profiles.geo',&
          status='old')

     read(io,*) i
     do i=1,EXPRO_n_exp
        read(io,*) EXPRO_geo(:,:,i)
     enddo

     close(io)

  endif
  !--------------------------------------------------------------

! Depricated ALPHA profiles functionality.

!  if (EXPRO_ctrl_alpha_profiles_flag == 1) then
!    call Alpha_mainsub
!  !--------------------------------------------------------------------------------
!  ! Mapping equillibrium and alpha profiles onto EXPRO variables from Alpha.
!  ! Density and temperature.
!    call cub_spline(rho_hat,T_i_rho,n_rho_grid,EXPRO_rho,EXPRO_ti(1,:),EXPRO_n_exp)
!    call cub_spline(rho_hat,T_e_rho,n_rho_grid,EXPRO_rho,EXPRO_te,EXPRO_n_exp)
!    call cub_spline(rho_hat,n_i_rho,n_rho_grid,EXPRO_rho,EXPRO_ni(1,:),EXPRO_n_exp)
!    call cub_spline(rho_hat,n_e_rho,n_rho_grid,EXPRO_rho,EXPRO_ne,EXPRO_n_exp)
!
!  ! Shaping and safety factor.
!    call cub_spline(rho_hat,q_rho,n_rho_grid,EXPRO_rho,EXPRO_q,EXPRO_n_exp)
!    call cub_spline(rho_hat,rmaj_rho,n_rho_grid,EXPRO_rho,EXPRO_rmaj,EXPRO_n_exp)
!    call cub_spline(rho_hat,rmin_rho,n_rho_grid,EXPRO_rho,EXPRO_rmin,EXPRO_n_exp)
!    call cub_spline(rho_hat,kappa_rho,n_rho_grid,EXPRO_rho,EXPRO_kappa,EXPRO_n_exp)
!    call cub_spline(rho_hat,delta_rho,n_rho_grid,EXPRO_rho,EXPRO_delta,EXPRO_n_exp)
!
!  ! Crossover energy, slowing-down density, effective temperature.
!    call cub_spline(rho_hat,E_c_hat_rho,n_rho_grid,EXPRO_rho,EXPRO_energy_c(2,:),EXPRO_n_exp)
!    call cub_spline(rho_hat,n_alpha_rho,n_rho_grid,EXPRO_rho,EXPRO_ni(2,:),EXPRO_n_exp)
!    if (use_t_equiv == 1) then
!      call cub_spline(rho_hat,T_alpha_equiv_rho,n_rho_grid,EXPRO_rho,EXPRO_ti(2,:),EXPRO_n_exp)
!    else
!      call cub_spline(rho_hat,E_alpha_grid,n_rho_grid,EXPRO_rho,EXPRO_ti(2,:),EXPRO_n_exp)
!    endif
!
!  else

    if (EXPRO_nion_file .gt. 1) then
      do i = 1, nion_max
        EXPRO_energy_c(i,:) = 33.0344836 * EXPRO_te(:) / EXPRO_ti(2,:)  ! Only appropriate for alpha particles.
      enddo
    endif

!  endif

  deallocate(dummy)

end subroutine EXPRO_read_driver
 
 
