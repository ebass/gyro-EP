!---------------------------------------------------------------
! cub_spline.f90
!2.08.21 GQLCGM  cub_spline_reverse.f90
!
! PURPOSE:
!  Take known points x(1:n),y(1:n) and perform cubic spline
!  interpolation at the node points xi(1:ni) to give yi(1:ni).
!
!  Use 'natural' cubic spline interpolation, where natural 
!  means the second derivative is zero at the boundaries.
!
!  INPUT  : x(1:n),y(1:n),n,xi(1:ni),ni
!
!  OUTPUT : yi(1:ni)
!
!  The only requirements are 
!  
!  1. Ordered data: x(i+1) > x(i) and xi(i+1) > xi(i).
!  2. Bounded data: xi(1) < x(1) an1d xi(ni) > x(n).
!
!cub_spline_reversee.f90  does the reveres by simple interpolation  !2.08.21 GQLCGM 
!----------------------------------------------------------------

subroutine cub_spline_reverse(x,y,n,xi,yi,ni)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  real :: x0 
  !
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: x,y
  !
  integer, intent(in) :: ni
  real, dimension(ni) :: xi,yi
  !
  !
  integer :: info
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Check to see that interpolated point is inside data interval
  !
  !2.08.21   !reverse
   if (xi(1) > x(1) .or. xi(ni) < x(n)) then
     print *,'Error in reverse interpolation'
     print *,'xi(1) > x(1)',xi(1),x(1)
     print *,'xi(ni) < x(n)',xi(ni),x(n)
  endif
  !-------------------------------------------------------------
  !
  yi(:)=0.
  do ii=1,ni  !out yi(ii)
   do i=1,n   !in y(i)
    if(xi(ii) .le. x(i+1) .and. xi(ii) .ge. x(i)) then 
            yi(ii) = y(i)+(y(i+1) - y(i))/(x(i+1)-x(i))*(xi(ii)-x(i))
    endif
   enddo
  enddo
  !-------------------------------------------------------------

end subroutine cub_spline_reverse

