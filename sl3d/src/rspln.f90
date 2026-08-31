subroutine rspln(i1,i2,x,y,q,f)
  !
  ! Computes cubic spline interpolation coefficients for y(x) between the
  ! grid points i1 and i2, saving them in q.  The interpolation is
  ! continuous with continuous first and second derivatives.  It agrees
  ! exactly with y at the grid points, and with the three point first
  ! derivatives at both end points (i1 and i2).  x must be monotonic, but
  ! if two successive values of x are equal a discontinuity is assumed and
  ! separate interpolation is done on each strictly monotonic segment.  The
  ! arrays must be dimensioned at least x(i2), y(i2), q(3,i2) and f(3,i2).
  ! f is working storage for rspln.
  !
  !                                                    -RPB
  !
  ! Translated from the original FORTRAN 66 rspln.f.  The algorithm and the
  ! default real kind are unchanged; the arithmetic IFs, the shared DO
  ! terminations and the EQUIVALENCE of y0 into the array yy have been
  ! rewritten.  In the original, yy held (y0,0,0) and was used to fill in
  ! the coefficients of a linear segment, which is what the assignments to
  ! q below now do directly.
  !
  implicit none

  integer, intent(in) :: i1,i2
  real, intent(in) :: x(*),y(*)
  real :: q(3,*),f(3,*)

  real, parameter :: small = 1.e-5
  real :: y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
  integer :: i,j,j1,j2,k,idisc
  logical :: first_is_single

  j1 = i1+1
  y0 = 0.

  ! bail out if there are less than two points total
  if(i2-i1 < 0) return
  first_is_single = (i2-i1 == 0)

  segments: do

     if(first_is_single) then

        j1 = j1-1
        j2 = i2-2
        first_is_single = .false.

     else

        ! search for discontinuities
        a0 = x(j1-1)
        idisc = 0
        do i = j1,i2
           b0 = a0
           a0 = x(i)
           if(abs((a0-b0)/max(a0,b0)) < small) then
              idisc = i
              exit
           end if
        end do
        j1 = j1-1
        if(idisc /= 0) then
           j2 = idisc-3
        else
           j2 = i2-2
        end if

     end if

     ! see if there are enough points to interpolate (at least three)
     if(j2+1-j1 > 0) then

        ! more than two points.  do spline interpolation
        a0 = 0.
        h = x(j1+1)-x(j1)
        h2 = x(j1+2)-x(j1)
        y0 = h*h2*(h2-h)
        h = h*h
        h2 = h2*h2

        ! calculate derivative at near end
        b0 = (y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
        b1 = b0

        ! explicitly reduce banded matrix to an upper banded matrix
        do i = j1,j2
           h = x(i+1)-x(i)
           y0 = y(i+1)-y(i)
           h2 = h*h
           ha = h-a0
           h2a = h-2.*a0
           h3a = 2.*h-3.*a0
           h2b = h2*b0
           q(1,i) = h2/ha
           q(2,i) = -ha/(h2a*h2)
           q(3,i) = -h*h2a/h3a
           f(1,i) = (y0-h*b0)/(h*ha)
           f(2,i) = (h2b-y0*(2.*h-a0))/(h*h2*h2a)
           f(3,i) = -(h2b-3.*y0*ha)/(h*h3a)
           a0 = q(3,i)
           b0 = f(3,i)
        end do

        ! take care of last two rows
        i = j2+1
        h = x(i+1)-x(i)
        y0 = y(i+1)-y(i)
        h2 = h*h
        ha = h-a0
        h2a = h*ha
        h2b = h2*b0-y0*(2.*h-a0)
        q(1,i) = h2/ha
        f(1,i) = (y0-h*b0)/h2a
        ha = x(j2)-x(i+1)
        y0 = -h*ha*(ha+h)
        ha = ha*ha

        ! calculate derivative at far end
        y0 = (y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
        q(3,i) = (y0*h2a+h2b)/(h*h2*(h-2.*a0))
        q(2,i) = f(1,i)-q(1,i)*q(3,i)

        ! solve upper banded matrix by reverse iteration
        do j = j1,j2
           k = i-1
           q(1,i) = f(3,k)-q(3,k)*q(2,i)
           q(3,k) = f(2,k)-q(2,k)*q(1,i)
           q(2,k) = f(1,k)-q(1,k)*q(3,k)
           i = k
        end do
        q(1,i) = b1

        ! fill in the last point with a linear extrapolation
        j2 = j2+2
        q(1,j2) = y0
        q(2,j2) = 0.
        q(3,j2) = 0.

     else if(j2+1-j1 == 0) then

        ! only two points.  use linear interpolation
        j2 = j2+2
        y0 = (y(j2)-y(j1))/(x(j2)-x(j1))
        q(1,j1) = y0
        q(2,j1) = 0.
        q(3,j1) = 0.
        q(1,j2) = y0
        q(2,j2) = 0.
        q(3,j2) = 0.

     else

        ! fill in the last point with a linear extrapolation
        j2 = j2+2
        q(1,j2) = y0
        q(2,j2) = 0.
        q(3,j2) = 0.

     end if

     ! see if this discontinuity is the last
     if(j2-i2 >= 0) return

     ! no.  go back for more
     j1 = j2+2
     if(j1-i2 > 0) then
        ! there is only one point left after the latest discontinuity
        q(1,i2) = y0
        q(2,i2) = 0.
        q(3,i2) = 0.
        return
     end if

  end do segments

end subroutine rspln
