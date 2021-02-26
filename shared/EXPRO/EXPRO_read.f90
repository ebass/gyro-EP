!--------------------------------------------------------------
! EXPRO_read.f90
!
! PURPOSE:
!  Read experimental profiles (serial).
!
! NOTES:
!  Variable naming corresponds to input.profiles 
!  documentation at 
!
!    http://fusion.gat.com/theory/input.profiles
!--------------------------------------------------------------

subroutine EXPRO_read

  use EXPRO_interface

  implicit none

  call EXPRO_read_driver
  call EXPRO_compute_derived

  call EXPRO_write_derived(1,'input.profiles.extra')  ! EMB 1-29-2021

end subroutine EXPRO_read

