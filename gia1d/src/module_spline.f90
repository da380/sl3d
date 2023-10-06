module module_spline


  ! This module contains the Numerical Recipies routines for constructing 
  ! and evaluating cubic splines. The routines have been overloaded for 
  ! double precision real and complex valued functions of a real argument.
  ! Routines have also been added that permit piecewise continuous functions
  ! to be used -- discontinuites are shown by a repetition of the independent
  ! variable in the arrays.


  interface tridag_ser
     module procedure tridag_ser_r,tridag_ser_c
  end interface

  interface spline
     module procedure spline_r,spline_c
  end interface

  interface spline_dis
     module procedure spline_dis_r,spline_dis_c
  end interface

  interface splint
     module procedure splint_r,splint_c
  end interface

  interface splint_dis
     module procedure splint_dis_r,splint_dis_c
  end interface

  interface splint_div1_dis
     module procedure splint_div1_dis_r
  end interface


contains


  subroutine spline_r(x,y,yp1,ypn,y2)
    use nrtype; use nrutil, only : assert_eq
    implicit none
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), intent(in) :: yp1,ypn
    real(dp), dimension(:), intent(out) :: y2
    integer(i4b) :: n
    real(dp), dimension(size(x)) :: a,b,c,r

    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=x(2:n)-x(1:n-1) 
    r(1:n-1)=6._dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0_dp
    b(n)=1.0_dp
    if(yp1 > 0.99e30_dp) then 
       r(1)=0._dp
       c(1)=0._dp
    else
       r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5_dp
    end if
    if(ypn > 0.99e30_dp) then 
       r(n)=0._dp
       a(n)=0._dp
    else
       r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5_dp
    end if
    call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    return
  end subroutine spline_r

  subroutine spline_c(x,y,yp1,ypn,y2)
    use nrtype; use nrutil, only : assert_eq
    implicit none
    real(dp), dimension(:), intent(in) :: x
    complex(dpc), dimension(:), intent(in) :: y
    real(dp), intent(in) :: yp1,ypn
    complex(dpc), dimension(:), intent(out) :: y2
    integer(i4b) :: n
    complex(dpc), dimension(size(x)) :: a,b,c,r

    n=assert_eq(size(x),size(y),size(y2),'spline')
    c(1:n-1)=(x(2:n)-x(1:n-1))*(1.0_dp,0.0_dp)
    r(1:n-1)=6._dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=(1.0_dp,0.0_dp)
    b(n)=(1.0_dp,0.0_dp)
    if(yp1 > 0.99e30_dp) then 
       r(1)=(0.0_dp,0.0_dp)
       c(1)=(0.0_dp,0.0_dp)
    else
       r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=(0.5_dp,0.0_dp)
    end if
    if(ypn > 0.99e30_dp) then 
       r(n)=(0.0_dp,0.0_dp)
       a(n)=(0.0_dp,0.0_dp)
    else
       r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=(0.5_dp,0.0_dp)
    end if
    call tridag_ser(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    return
  end subroutine spline_c


  subroutine spline_dis_r(x,y,yp1,ypn,y2)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: x,y
    real(dp), intent(in) :: yp1,ypn
    real(dp), dimension(:), intent(out) :: y2
    integer(i4b), dimension(:), allocatable :: ind
    integer(i4b) :: i,i1,i2,kount,n
    real(dp), parameter :: tiny=1.0e-4_dp

    n=size(x)
    kount=0
    do i=1,n-1
       if(abs(x(i+1)-x(i)) < tiny) then
          kount=kount+1
          print *, kount
       end if
    end do
    allocate(ind(kount+1))
    kount=0
    do i=1,n-1
       if(abs(x(i+1)-x(i)) < tiny) then
          kount=kount+1
          ind(kount)=i
       end if
    end do
    kount=kount+1
    ind(kount)=n
    do i=1,kount-1
       if(i == 1) then
          i1=1
       else
          i1=ind(i-1)+1
       end if
       i2=ind(i)
       call spline(x(i1:i2),y(i1:i2),yp1,ypn,y2(i1:i2))
    end do
   return
  end subroutine spline_dis_r


  subroutine spline_dis_c(x,y,yp1,ypn,y2)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: x
    complex(dpc), dimension(:), intent(in) :: y
    real(dp), intent(in) :: yp1,ypn
    complex(dpc), dimension(:), intent(out) :: y2
    integer(i4b), dimension(:), allocatable :: ind
    integer(i4b) :: i,i1,i2,kount,n
    real(dp), parameter :: tiny=1.0e-4_dp

    n=size(x)
    kount=0
    do i=1,n-1
       if(abs(x(i+1)-x(i)) < tiny) then
          kount=kount+1
       end if
    end do
    allocate(ind(kount+1))
    kount=0
    do i=1,n-1
       if(abs(x(i+1)-x(i)) < tiny) then
          kount=kount+1
          ind(kount)=i
       end if
    end do
    kount=kount+1
    ind(kount)=n
    do i=1,kount-1
       if(i == 1) then
          i1=1
       else
          i1=ind(i-1)+1
       end if
       i2=ind(i)
       call spline(x(i1:i2),y(i1:i2),yp1,ypn,y2(i1:i2))
    end do
    return
  end subroutine spline_dis_c



  function klo_find(xa,x)
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: xa
    real(dp), intent(in) :: x
    integer(i4b) :: klo, klo_find,n
    n=size(xa)
    klo_find=max(min(locate(xa,x),n-1),1)
    return
  end function klo_find



  function splint_c(xa,ya,y2a,x,kloin)
    use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa
    complex(dpc), dimension(:), intent(in) :: ya,y2a
    real(dp), intent(in) :: x
    complex(dpc) :: splint_c
    integer(i4b), intent(in), optional :: kloin
    integer(i4b) :: khi,klo,n
    real(dp) :: a,b,h

    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    if(present(kloin)) then
       klo=kloin
    else
       klo=max(min(locate(xa,x),n-1),1)
    end if
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_c=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
         *(h**2)/6.0_dp
    return
  end function splint_c
  
  function splint_r(xa,ya,y2a,x,kloin)
    use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    real(dp) :: splint_r
    integer(i4b), intent(in), optional :: kloin
    integer(i4b) :: khi,klo,n
    real(dp) :: a,b,h

    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    if(present(kloin)) then
       klo=kloin
    else
       klo=max(min(locate(xa,x),n-1),1)
    end if
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_r=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
         *(h**2)/6.0_dp
    return
  end function splint_r


  function splint_dis_r(xa,ya,y2a,x,klo)
   use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    real(dp) :: splint_dis_r
    integer(i4b), intent(in) :: klo
    integer(i4b) :: khi
    real(dp) :: a,b,h

    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) then
       splint_dis_r=ya(klo)
       return
    end if
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_dis_r=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
         *(h**2)/6.0_dp
    return
  end function splint_dis_r


  function splint_dis_c(xa,ya,y2a,x,klo)
   use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa
    complex(dpc), dimension(:), intent(in) :: ya,y2a
    real(dp), intent(in) :: x
    complex(dpc) :: splint_dis_c
    integer(i4b), intent(in) :: klo
    integer(i4b) :: khi
    real(dp) :: a,b,h

    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) then
       splint_dis_c=ya(khi)
       return
    end if
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_dis_c=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) &
         *(h**2)/6.0_dp
    return
  end function splint_dis_c


  function splint_div1_dis_r(xa,ya,y2a,x,klo)
   use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    real(dp) :: splint_div1_dis_r
    integer(i4b), intent(in) :: klo
    integer(i4b) :: khi,klot
    real(dp) :: a,b,h

    klot=klo
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) then
       if(klo == 1) then
          splint_div1_dis_r=0.0_dp
          return
       else
          klot=klo-1
          khi=khi-1
          h=xa(khi)-xa(klot)
       end if
    end if
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klot))/h
    splint_div1_dis_r=(ya(khi)-ya(klot))/h+h*(y2a(khi)*(3.0_dp*b**2-1.0_dp) &
         -y2a(klot)*(3.0_dp*a**2-1.0_dp))/6.0_dp
    return
  end function splint_div1_dis_r


  

  
  function splint_div1(xa,ya,y2a,x,kloin)
    use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    integer(i4b), intent(in), optional :: kloin
    real(dp) :: splint_div1
    integer(i4b) :: khi,klo,n
    real(dp) :: a,b,h

    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    if(present(kloin)) then
       klo=kloin
    else
       klo=max(min(locate(xa,x),n-1),1)
    end if
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_div1=(ya(khi)-ya(klo))/h+h*(y2a(khi)*(3.0_dp*b**2-1.0_dp) &
         -y2a(klo)*(3.0_dp*a**2-1.0_dp))/6.0_dp
    return
  end function splint_div1


  function splint_div2(xa,ya,y2a,x,kloin)
    use nrtype; use nrutil, only: assert_eq,nrerror
    implicit none
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), intent(in) :: x
    integer(i4b), intent(in), optional :: kloin
    real(dp) :: splint_div2
    integer(i4b) :: khi,klo,n
    real(dp) :: a,b,h

    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    if(present(kloin)) then
       klo=kloin
    else
       klo=max(min(locate(xa,x),n-1),1)
    end if
    khi=klo+1
    h=xa(khi)-xa(klo)
    if(h == 0.0_dp) call nrerror('bad xa input in splint')
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    splint_div2=a*y2a(klo)+b*y2a(khi)
    return
  end function splint_div2
  



  subroutine spline_int(xa,ya,y2a,sum,seed)
    use nrtype; use nrutil
    implicit none
    integer(i4b) :: i,j,n
    real(dp), dimension(:), intent(in) :: xa,ya,y2a
    real(dp), dimension(:), intent(out) :: sum
    real(dp), intent(in), optional :: seed
    real(dp) :: del
    
    if(present(seed)) then
       sum(1)=seed
    else
       sum(1)=0.0_dp
    end  if
    n=assert_eq(size(xa),size(ya),size(y2a),'spline_int')
    do i=2,n
       del=xa(i)-xa(i-1)
       sum(i)=sum(i-1)+del*(ya(i)+ya(i-1))/2.0_dp + &
            del**3*(y2a(i)+y2a(i-1))/24.0_dp
    end do
    return
  end subroutine spline_int


  subroutine tridag_ser_r(a,b,c,r,u) 
    use nrtype; use nrutil, only : assert_eq,nrerror
    implicit none
      real(dp), dimension(:), intent(in) :: a,b,c,r
      real(dp), dimension(:), intent(out) :: u
      real(dp), dimension(size(b)) :: gam 
      integer(dp) :: n,j
      real(dp) :: bet

      n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
      bet=b(1)
      if(bet == 0.0_dp) call nrerror('tridag_ser: error at code stage 1')
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if(bet == 0.0_dp) &
              call nrerror('tridag_ser: error at code stage 2')
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
    end subroutine tridag_ser_r

    subroutine tridag_ser_c(a,b,c,r,u) 
      use nrtype; use nrutil, only : assert_eq,nrerror
      implicit none
      complex(dpc), dimension(:), intent(in) :: a,b,c,r
      complex(dpc), dimension(:), intent(out) :: u
      complex(dpc), dimension(size(b)) :: gam 
      integer(dp) :: n,j
      complex(dpc) :: bet

      n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
      bet=b(1)
      if(abs(bet) == 0.0_dp) call nrerror('tridag_ser: error at code stage 1')
      u(1)=r(1)/bet
      do j=2,n
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if(abs(bet) == 0.0_dp) &
              call nrerror('tridag_ser: error at code stage 2')
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
    end subroutine tridag_ser_c


    function locate(xx,x)
      use nrtype
      implicit none
      real(dp), dimension(:), intent(in) :: xx
      real(dp), intent(in) :: x
      integer(i4b) :: locate
      integer(i4b) :: n,jl,jm,ju
      logical :: ascnd

      n=size(xx)
      ascnd= (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do 
         if(ju-jl <= 1) exit
         jm=(ju+jl)/2
         if(ascnd .eqv. (x >= xx(jm))) then
            jl=jm
         else
            ju=jm
         end if
      end do
      if(x == xx(1)) then
         locate=1
      else if(x == xx(n)) then
         locate=n-1
      else
         locate=jl
      end if
      return
    end function locate


  end module module_spline
