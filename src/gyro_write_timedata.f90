!------------------------------------------------------
! gyro_write_timedata.F90
!
! PURPOSE:
!  This is the master file controlling output of
!  data on a per-timestep basis.  This file also 
!  contains the MPI IO routines 
!
!  - write_distributed_real
!  - write_distributed_complex
!  - write_local_real
!-----------------------------------------------------

subroutine gyro_write_timedata

  use gyro_globals
  use GQLCGM_globals !GQLCGM 2.08.21
  use mpi

  !---------------------------------------------------
  implicit none
  !
  real, dimension(:,:,:), allocatable :: a3
  !
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: n_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: e_plot
  complex, dimension(n_theta_plot,n_x,n_kinetic) :: v_plot
  !---------------------------------------------------

  !---------------------------------------------------
  ! Timestep data:
  !
  if (i_proc == 0) then
     call gyro_write_step(trim(path)//'out.gyro.t',1)
  endif
  !---------------------------------------------------

  !--------------------------------------------------
  ! Output of field-like quantities:
  !
  if (plot_n_flag+plot_e_flag+plot_v_flag > 0) then
     n_plot(:,:,:) = moments_plot(:,:,:,1)
     e_plot(:,:,:) = moments_plot(:,:,:,2)
     v_plot(:,:,:) = moments_plot(:,:,:,3)
  endif
  !
  if (plot_u_flag == 1) then

     ! POTENTIALS

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_u',&
          10,&
          size(phi_plot(:,:,1:n_field)),&
          phi_plot(:,:,1:n_field))

  endif !u_flag==1

  if (plot_epar_flag == 1) then

     ! PARALLEL ELECTRIC FIELD

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_epar',&
          10,&
          size(phi_plot(:,:,n_field+1)),&
          phi_plot(:,:,n_field+1))

  endif !epar_flag==1

  if (plot_n_flag == 1) then

     ! DENSITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_n',&
          10,&
          size(n_plot),&
          n_plot)

  endif !n_flag ==1 

  if (plot_e_flag == 1) then

     ! ENERGY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_e',&
          10,&
          size(e_plot),&
          e_plot)

  endif !e_flag==1

  if (plot_v_flag == 1) then

     ! PARALLEL VELOCITY

     call write_distributed_complex(&
          trim(path)//'out.gyro.moment_v',&
          10,&
          size(v_plot),&
          v_plot)

  endif !v_flag==1

  !--------------------------------------------------

  call gyro_kxky_spectrum

  call write_distributed_real(&
       trim(path)//'out.gyro.kxkyspec',&
       10,&
       size(kxkyspec),&
       kxkyspec)

  if (i_proc == 0) then
     call write_local_real(&
          trim(path)//'out.gyro.k_perp_squared',&
          10,&
          size(k_perp_squared),&
          k_perp_squared)
  endif

  call gyro_field_fluxave

  !-------------------------------------------------------------------
  ! Calculation of fundamental nonlinear fluxes
  !
  call gyro_nonlinear_flux
  call gyro_gbflux
  !-------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Output specific to linear/nonlinear operation:
  !
  if (nonlinear_flag == 0) then

     !=============
     ! BEGIN LINEAR 
     !=============

     call gyro_write_freq(trim(path)//'out.gyro.freq',10)

     if (plot_u_flag == 1) then        

        ! PHI
        call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_phi',10,1,0)

        if (n_field > 1) then
           ! A_PARALLEL 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_a',10,2,0)

        endif

        if (n_field > 2) then
           ! B_PARALLEL 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_aperp',10,3,0)

        endif

        ! E_PARALLEL
        if (eparallel_plot_flag == 1) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_epar',10,n_field+1,0)

        endif

     endif

     if (plot_n_flag == 1) then

        ! DENSITY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_ion',10,5,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_n_elec',10,5,indx_e)

        endif
     endif

     if (plot_e_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_ion',10,6,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_e_elec',10,6,indx_e)

        endif
     endif

     if (plot_v_flag == 1) then

        ! ENERGY
        if (electron_method /= 3) then

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_ion',10,7,1)

        endif
        if (electron_method > 1) then 

           call gyro_ballooning_mode(trim(path)//'out.gyro.balloon_v_elec',10,7,indx_e)

        endif
     endif

     !-----------------------------------------------------------------
     ! Distribution function data:
     !
     if (n_proc == 1 .and. n_n == 1 .and. dist_print == 1) then
        call gyro_write_h(trim(path)//'out.gyro.hp',trim(path)//'out.gyro.ht',10,11)
     endif
     !-----------------------------------------------------------------

     ! GQLCGM  01.27.21
  if(data_step == 0) then
   gbflux(:,:,:) = 0.0
   gbflux_mom(:,:) = 0.0
   gbflux_i(:,:,:,:) = 0.0

   call send_line('data_step=0, gbflux,gbflux_mom,gbflux_i zeroed at t=0')
  endif

     if (i_proc == 0 .and. lindiff_method > 1) then

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',&
                10,size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',&
                10,size(gbflux_i_trapped),gbflux_i_trapped)
        endif

     endif !i_proc ==0 and lindiff >1 

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
     if(i_do_GQLCGM .ge. 1) then  !GQLCGM 2.08.21
             call send_line('**[ GQLCGM_advance]')
        if(i_proc ==0) then
           print *, data_step, t_current
        endif
        if(t_current .gt. 0.) then
          if(i_proc ==0) then
           print *, 'update EP and  GQLCGM reset field matrix '
           print *, '------------------------------------------------------'
          endif
!get flux_EP_s(:,n_x) and flux_EP(:,n_grid_EP)  !GQLCGM 2.10.21
      do i_sEP =1,i_do_GQLCGM

       if(i_proc ==0) then
        print *, 'den_s(1+i_do_GQLCGM,ir_norm)=',den_s(1+i_do_GQLCGM,ir_norm)
        print *, 'den_EP_s(i_sEP,ir_norm)/den_norm=', den_EP_s(i_sEP,ir_norm)/den_norm
        print *, 'den_sd_EP_s(i_sEP,ir_norm)/den_norm=', den_sd_EP_s(i_sEP,ir_norm)/den_norm

        print *, 'alpha_s(1+i_do_GQLCGM,ir_norm)=',alpha_s(1+i_do_GQLCGM,ir_norm)

        print *, 'dlnndr_s(1+i_do_GQLCGM,ir_norm)=',dlnndr_s(1+i_do_GQLCGM,ir_norm)
        print *, 'dlnndr_EP_s(i_sEP,ir_norm)=',dlnndr_EP_s(i_sEP,ir_norm)
       endif

       flux_EP_s(i_sEP,:) = gbflux_i(1+i_sEP,1,1,:) + gbflux_i(1+i_sEP,2,1,:)
       if(i_proc ==0) print *, 'max flux_EP_s(i_sEP,:)=',maxval(flux_EP_s(i_sEP,:))
       flux_EP_s(i_sEP,:) = flux_EP_s(i_sEP,:)*flux_gbnorm_EP
       if(i_proc ==0) print *, 'max metric flux_EP_s(i_sEP,:)=',maxval(flux_EP_s(i_sEP,:))
       !small fine grid to course large grid
         r_EP(:) = r_EP(:)/a_meters  !2.12.21
       call cub_spline_reverse(r_s,flux_EP_s(i_sEP,:),n_x,r_EP,flux_EP(i_sEP,:),n_grid_EP)
         r_EP(:) = r_EP(:)*a_meters
       if(i_proc ==0) print *, 'max flux_EP(i_sEP,:)=',maxval(flux_EP(i_sEP,:))


       if(i_proc==0 ) then
           if(modulo(data_step,10)==0) then
        !!!!   if(modulo(data_step,50)==0) then
               print *, 'metric flux_EP_s'
               print *, '------------------------------------------------------'
               do i=1,n_x
                print *, r_s(i), flux_EP_s(i_sEP,i)
               enddo

               print *, 'metric flux_EP'
               print *, '------------------------------------------------------'
               do i_EP=1,n_grid_EP
                print *, r_EP(i_EP)/a_meters, flux_EP(i_sEP,i_EP)
               enddo
          endif
       endif

      enddo !i_sEP


!EP tranport from last data_step 
 
    call GQLCGM_EPtransport

!    den_EP(:,:) =den_sd_EP(:,:)  !override GQLCGM_EPtransport  2.23.21  tested i_update_gyro_GQLCGM OK 
    
    !add slow down by relaxation GQLCGM 2.22.21
     do i_sEP =1,i_do_GQLCGM
!       den_EP(i_sEP,:) = 0.1*den_EP(i_sEP,:) + 0.9*den_prev_EP(i_sEP,:) 
     den_EP(i_sEP,:) = 0.05*den_EP(i_sEP,:) + 0.95*den_prev_EP(i_sEP,:)
!     den_EP(i_sEP,:) = 0.04*den_EP(i_sEP,:) + 0.96*den_prev_EP(i_sEP,:)
!      den_EP(i_sEP,:) = 0.03*den_EP(i_sEP,:) + 0.97*den_prev_EP(i_sEP,:)
!      den_EP(i_sEP,:) = 0.02*den_EP(i_sEP,:) + 0.98*den_prev_EP(i_sEP,:)
!      den_EP(i_sEP,:) = 0.01*den_EP(i_sEP,:) + 0.99*den_prev_EP(i_sEP,:) 
!       den_EP(i_sEP,:) = 0.001*den_EP(i_sEP,:) + 0.999*den_prev_EP(i_sEP,:) 
       den_prev_EP(i_sEP,:) = den_EP(i_sEP,:)
     enddo !i_sEP

    if(i_proc ==0) then
     !!!if(modulo(data_step,50)==0) then
      if(modulo(data_step,10)==0) then
     !!!if(modulo(data_step,1)==0) then
    do i_sEP =1,i_do_GQLCGM
               print *, 'den_EP den_sd_EP '
               print *, '------------------------------------------------------'
               do i_EP=1,n_grid_EP
                print *, i_EP,den_EP(i_sEP,i_EP), den_sd_EP(i_sEP,i_EP)
               enddo
    enddo !i_sEP
       print *, 'error_EP(i_sEP)=',error_EP(:)
     endif
    endif

    i_update_gyro_GQLCGM = 1   !2.23.21


  if (i_update_gyro_GQLCGM .gt. 0) then  !2.22.21

    do i_sEP =1,i_do_GQLCGM
        r_EP(:) = r_EP(:)/a_meters  !2.12.21
       call cub_spline(r_EP,den_EP(i_sEP,:),n_grid_EP,r_s,den_EP_s(i_sEP,:),n_x)   !note metric units den_EP and deb_EP_s
        r_EP(:) = r_EP(:)*a_meters
    enddo !i_sEP


!set EP factors for field matrix reset

         do i_sEP=1,i_do_GQLCGM

          den_s(1+i_sEP, :) = den_EP_s(i_sEP,:)/den_norm

          do i=1,n_x
           alpha_s(1+i_sEP,i) = den_EP_s(i_sEP,i)/tem_s(1+i_sEP,i)/den_norm
          enddo

          do i=2,n_x-1
            dlnndr_EP_s(i_sEP,i) = -(den_EP_s(i_sEP,i+1)-den_EP_s(i_sEP,i-1))/(r_s(i+1)-r_s(i-1))/den_EP_s(i_sEP,i)
          enddo
          dlnndr_EP_s(i_sEP,1) = dlnndr_EP_s(i_sEP,2)
          dlnndr_EP_s(i_sEP,n_x) = dlnndr_EP_s(i_sEP,n_x-1)

         enddo !i_sEP

  endif !i_update_gyro_GQLCGM

!field matrix reset
        call gyro_make_poisson_matrix
        call gyro_make_ampere_matrix
        call gyro_make_maxwell_matrix

        if(i_proc ==0) print *, '------------------------------------------------------'

        endif !t_current .gt. 0.
        call send_line('**[ GQLCGM_advance_finished]')

     endif !i_do_GQLCGM .ge. 1
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

     !=============
     ! END LINEAR 
     !=============

  else

     !================
     ! BEGIN NONLINEAR 
     !================

     call write_distributed_real(&
          trim(path)//'out.gyro.gbflux_n',&
          10,&
          size(gbflux_n),&
          gbflux_n)

     if (nonlinear_transfer_flag == 1) then
        call write_distributed_real(&
             trim(path)//'out.gyro.nl_transfer',&
             10,&
             size(nl_transfer),&
             nl_transfer)
     endif !nonlinear_transfer_flag ==1 

     if (i_proc == 0) then

        call write_local_real(trim(path)//'out.gyro.field_rms',10,size(ave_phi),ave_phi)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux',10,size(gbflux),gbflux)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_i',10,size(gbflux_i),gbflux_i)

        call write_local_real( &
             trim(path)//'out.gyro.gbflux_mom',10,size(gbflux_mom),gbflux_mom)
        call write_local_real( &
             trim(path)//'out.gyro.gbflux_exc',10,size(gbflux_exc),gbflux_exc)

        if (trapdiff_flag == 1) then
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_trapped',10,&
                size(gbflux_trapped),gbflux_trapped)
           call write_local_real( &
                trim(path)//'out.gyro.gbflux_i_trapped',10,&
                size(gbflux_i_trapped),gbflux_i_trapped)
        endif !trapdiff_flag == 1

        call write_local_real( &
             trim(path)//'out.gyro.zerobar',10,&
             size(field_fluxave),transpose(field_fluxave))


        allocate(a3(n_kinetic,4,n_x))
        do i=1,n_x
           a3(:,1,i) = h0_n(:,i)
           a3(:,2,i) = h0_e(:,i)
           a3(:,3,i) = source_n(:,i)
           a3(:,4,i) = source_e(:,i)
        enddo

        call write_local_real( &
             trim(path)//'out.gyro.source',10,size(a3),a3)

        call write_local_real( &
             trim(path)//'out.gyro.moment_zero',10,&
             size(moments_zero_plot),moments_zero_plot)

        deallocate(a3)

     endif! i_proc ==0

     !================
     ! END NONLINEAR 
     !================

  endif
  !-------------------------------------------------------------------

  call gyro_write_error(trim(path)//'out.gyro.error',10)

  !------------------------------------------------------------
  ! Entropy diagnostics
  !
  if (entropy_flag == 1) then
     call gyro_entropy 
     if (i_proc == 0) then 
        call write_local_real(&
             trim(path)//'out.gyro.entropy.out',10,size(entropy),entropy)
     endif
  endif
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Velocity-space diagnostics
  !
  if (velocity_output_flag == 1) then
     call gyro_nonlinear_flux_velocity
     call write_distributed_real(&
          trim(path)//'out.gyro.flux_velocity',&
          10,&
          size(nonlinear_flux_velocity),&
          nonlinear_flux_velocity)
  endif
  !
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! Write precision-monitoring data
  !
  call gyro_write_precision(10,sum(abs(gbflux)))
  !------------------------------------------------------------

  !------------------------------------------------------------
  call gyro_write_timers(trim(path)//'out.gyro.timing',10)
  !------------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) print *,'[gyro_write_timedata done]'

end subroutine gyro_write_timedata

!===========================================================================
!------------------------------------------------------
! write_distributed_real.f90
!
! PURPOSE:
!  Control merged output of distributed real array.
!------------------------------------------------------

subroutine write_distributed_real(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err,&
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  real :: fn_recv(n_fn)
  !------------------------------------------------------

  select case (io_control)

  case(0)

     return 

  case(1)

     ! Initial open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Output

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_PRECISION,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Reposition after restart

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_real

!------------------------------------------------------
! write_distributed_complex.f90
!
! PURPOSE:
!  Control merged output of complex distributed array.
!------------------------------------------------------

subroutine write_distributed_complex(datafile,io,n_fn,fn)

  use mpi
  use gyro_globals, only : &
       n_n,&
       n_n_1,&
       n_proc_1,&
       recv_status,&
       data_step,&
       GYRO_COMM_WORLD,&
       i_proc,&
       i_err, &
       io_control, &
       fmtstr

  !------------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  complex, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  integer :: i_group_send
  integer :: i_send
  integer :: in
  !
  complex :: fn_recv(n_fn)
  !------------------------------------------------------


  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     if (i_proc == 0) then
        open(unit=io,file=datafile,status='replace')
        close(io)
     endif

  case(2)

     ! Append

     if (i_proc == 0) &
          open(unit=io,file=datafile,status='old',position='append')

     do in=1,n_n

        !-----------------------------------------
        ! Subgroup collector:
        !
        i_group_send = (in-1)/n_n_1

        if (i_group_send /= 0) then

           i_send = i_group_send*n_proc_1

           if (i_proc == 0) then

              call MPI_RECV(fn_recv,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   i_send,&
                   in,&
                   GYRO_COMM_WORLD,&
                   recv_status,&
                   i_err)

           else if (i_proc == i_send) then

              call MPI_SEND(fn,&
                   n_fn,&
                   MPI_DOUBLE_COMPLEX,&
                   0,&
                   in,&
                   GYRO_COMM_WORLD,&
                   i_err)

           endif

        else

           fn_recv(:) = fn(:)

        endif
        !
        !-----------------------------------------

        if (i_proc == 0) then

           write(io,fmtstr) fn_recv(:)

        endif

     enddo ! in

     if (i_proc == 0) close(io)

  case(3)

     ! Rewind

     if (i_proc == 0) then

        open(unit=io,file=datafile,status='old')
        do data_loop=0,data_step

           do in=1,n_n
              read(io,fmtstr) fn_recv(:)
           enddo

        enddo ! data_loop

        endfile(io)
        close(io)

     endif

  end select

end subroutine write_distributed_complex

!------------------------------------------------------
! write_local_real.f90
!
! PURPOSE:
!  This routine write a vector of nondistributed reals.
!------------------------------------------------------

subroutine write_local_real(datafile,io,n_fn,fn)

  use gyro_globals, only : &
       data_step, &
       io_control, &
       fmtstr

  !---------------------------------------------------
  implicit none
  !
  character (len=*), intent(in) :: datafile
  integer, intent(in) :: io
  integer, intent(in) :: n_fn
  real, intent(in) :: fn(n_fn)
  !
  integer :: data_loop
  real :: dummy(n_fn)
  !---------------------------------------------------

  select case (io_control)

  case(0)

     return

  case(1)

     ! Open

     open(unit=io,file=datafile,status='replace')
     close(io)

  case(2)

     ! Append

     open(unit=io,file=datafile,status='old',position='append')
     write(io,fmtstr)  fn(:)
     close(io)

  case(3)

     ! Rewind

     open(unit=io,file=datafile,status='old')

     do data_loop=0,data_step
        read(io,fmtstr) dummy(:)
     enddo

     endfile(io)
     close(io)

  end select

end subroutine write_local_real

