function rsple(i1,i2,x,y,q,s)
  !
  ! Returns the value of the function y(x) evaluated at the point s, using
  ! the cubic spline coefficients computed by rspln and saved in q.  If s
  ! lies outside the interval (x(i1),x(i2)), rsple extrapolates using the
  ! first or last interpolation polynomial.  The arrays must be dimensioned
  ! at least x(i2), y(i2) and q(3,i2).
  !
  !                                                    -RPB
  !
  ! Translated from the original FORTRAN 66 rsple.f.  The algorithm, the
  ! default real kind and the bracketing index saved between calls are all
  ! unchanged; only the arithmetic IF control flow has been rewritten.
  !
  implicit none

  integer, intent(in) :: i1,i2
  real, intent(in) :: x(*),y(*),q(3,*),s
  real :: rsple

  ! the bracketing index is remembered between calls, as in the original
  integer, save :: i = 1
  integer :: ii
  real :: h

  ii = i2-1

  ! guarantee i within bounds
  i = max(i,i1)
  i = min(i,ii)

  if(x(i2) < x(i1)) then

     ! x is decreasing.  change i as necessary
     outer_dec: do
        if(s > x(i)) then
           i = i-1
           if(i < i1) then
              i = i1
              exit outer_dec
           end if
           if(i == i1) exit outer_dec
           cycle outer_dec
        end if
        do
           if(s >= x(i+1)) exit outer_dec
           i = i+1
           if(i > ii) then
              i = ii
              exit outer_dec
           end if
           if(i == ii) exit outer_dec
        end do
     end do outer_dec

  else

     ! x is increasing.  change i as necessary
     outer_inc: do
        if(s > x(i+1)) then
           i = i+1
           if(i > ii) then
              i = ii
              exit outer_inc
           end if
           if(i == ii) exit outer_inc
           cycle outer_inc
        end if
        do
           if(s >= x(i)) exit outer_inc
           i = i-1
           if(i < i1) then
              i = i1
              exit outer_inc
           end if
           if(i == i1) exit outer_inc
        end do
     end do outer_inc

  end if

  ! calculate rsple using the spline coefficients in y and q
  h = s-x(i)
  rsple = y(i)+h*(q(1,i)+h*(q(2,i)+h*q(3,i)))

end function rsple
