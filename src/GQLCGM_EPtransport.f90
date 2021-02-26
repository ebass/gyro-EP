!--------------------------------------------------------
! GQLCGM_EPtransport.f90
!
! PURPOSE:
!  Iterate  GQLCGM_transport  update den_EP
!--------------------------------------------------------

subroutine GQLCGM_EPtransport

  use mpi
  use GQLCGM_globals
  use gyro_globals

  !----------------------------------------
  implicit none
  !----------------------------------------


  real , dimension(2,n_grid_EP) :: D_EP
  real , dimension(2,n_grid_EP) :: flux_net_source
  real , dimension(2,n_grid_EP) :: den_EP_temp
  integer :: it_EP

  do i_sEP=1,i_do_GQLCGM

!total AE + bkg diffucivity
   do i_EP = 2,n_grid_EP-1
    D_EP(i_sEP,i_EP) = -flux_EP(i_sEP,i_EP)/(den_prev_EP(i_sEP,i_EP+1)-den_prev_EP(i_sEP,i_EP-1))*(r_EP(i_EP+1)-r_EP(i_EP-1)) &
                               + D_EP_bkg(i_sEP,i_EP)
   enddo
    D_EP(i_sEP,1) = D_EP(i_sEP,2)
    D_EP(i_sEP,n_grid_EP) = D_EP(i_sEP,n_grid_EP-1)

!temporary overrid 2.15.21
!2.18.21    D_EP(i_sEP,:) = D_EP_bkg(i_sEP,:)

   do it_EP = 1,10

!find net source flux
    flux_net_source(i_sEP,1) = 0.
    do i_EP=2,n_grid_EP
      flux_net_source(i_sEP,i_EP) = Vprime_EP(i_EP-1)/Vprime_EP(i_EP)*flux_net_source(i_sEP,i_EP-1)  + &
              0.5/Vprime_EP(i_EP)*(r_EP(i_EP)-r_EP(i_EP-1)) &
                 *(den_sd_EP(i_sEP,i_EP)/tau_sd_EP(i_sEP,i_EP)*(1.-den_EP(i_sEP,i_EP)/den_sd_EP(i_sEP,i_EP)) &
                 + den_sd_EP(i_sEP,i_EP-1)/tau_sd_EP(i_sEP,i_EP-1)*(1.-den_EP(i_sEP,i_EP-1)/den_sd_EP(i_sEP,i_EP-1)) )
    enddo

!    do i_EP = 1, n_grid_EP
!     flow_EP(i_sEP,i_EP) = Vprime_EP(i_EP)*flux_net_source(i_sEP,i_EP)
!    enddo

!   if(i_proc ==0) then
!     print *, '-----------------------------------------------------------------'
!     print *, 'r_EP, flow_EP'
!    do i_EP = 1, n_grid_EP
!     print *, r_EP(i_EP)/a_meters, flow_EP(i_sEP,i_EP)
!    enddo
!   endif

!update transported EP density
!    den_EP_temp(i_SEP,n_grid_EP) = 0.01*den_sd_EP(i_SEP,n_grid_EP)  !consistent with fast edge orbit loss
    den_EP_temp(i_SEP,n_grid_EP) = den_sd_EP(i_SEP,n_grid_EP)  
    do i_EP=n_grid_EP-1,1,-1
     den_EP_temp(i_sEP,i_EP) = den_EP_temp(i_sEP,i_EP+1) + (r_EP(i_EP+1)-r_EP(i_EP))  &
              *(flux_net_source(i_sEP,i_EP)+flux_net_source(i_sEP,i_EP+1)) &
              /(D_EP(i_sEP,i_EP) + D_EP(i_sEP,i_EP+1))
    enddo

    den_EP(i_sEP,:) = 0.1*den_EP_temp(i_sEP,:) + 0.9*den_EP(i_sEP,:)  !2.16.21

   enddo !it_EP

!   if(i_proc ==0) then
!     print *, '-----------------------------------------------------------------'
!     print *, 'i_EP, den_EP(i_sEP,i_EP),D_EP(i_sEP,i_EP)'
!    do i_EP = 1, n_grid_EP
!     print *, i_EP, den_EP(i_sEP,i_EP),D_EP(i_sEP,i_EP)
!    enddo
!   endif

!compute covergence error
   error_EP(i_sEP)=0.
   do i_EP=1,n_grid_EP
     error_EP(i_sEP) = error_EP(i_sEP) + ((den_EP(i_sEP,i_EP) - den_prev_EP(i_sEP,i_EP))/den_prev_EP(i_sEP,i_EP))**2
   enddo
   error_EP(i_sEP) = error_EP(i_sEP)/float(n_grid_EP)
   error_EP(i_sEP) = sqrt(error_EP(i_sEP))

!   if(i_proc == 0) print *, 'error_EP(i_sEP)=',error_EP(:)

!start for nextstep
!2.22.21 move to after call to GQLCGM_EPtransport
!den_prev_EP(i_sEP,:) = den_EP(i_sEP,:)



  enddo !i_sEP


end subroutine GQLCGM_EPtransport
