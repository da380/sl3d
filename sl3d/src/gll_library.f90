!=======================================================================
!
!  Library to compute the Gauss-Lobatto-Legendre points and weights
!  Based on Gauss-Lobatto routines from M.I.T.
!  Department of Mechanical Engineering
!
!=======================================================================

function endw1(n,alpha,beta)
  use nrtype
  implicit none

  ! function value
  real(dp) :: endw1

  ! function inputs
  integer(i4b), intent(in) ::  n
  real(dp), intent(in) ::  alpha,beta


  ! local variables
  real(dp), parameter :: zero=0.0_dp,one=1.0_dp,two=2.0_dp, &
  three=3.0_dp,four=4.0_dp
  real(dp) ::  apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  real(dp), external :: gammaf
  integer(i4b) ::  i

  f3 = zero
  apb   = alpha+beta
  if (n == 0) then
     endw1 = zero
     return
  endif
  f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
     endw1 = f1
     return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
     endw1 = f2
     return
  endif
  do i=3,n
     di   = real(i-1)
     abn  = alpha+beta+di
     abnn = abn+di
     a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
     a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
     a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
     f3   =  -(a2*f2+a1*f1)/a3
     f1   = f2
     f2   = f3
  enddo
  endw1  = f3
  
end function endw1

!
!=======================================================================
!

function endw2(n,alpha,beta)
  use nrtype  
  implicit none

  ! function value
  real(dp) :: endw2

  ! function inputs
  integer(i4b), intent(in) ::  n
  real(dp), intent(in) ::  alpha,beta

  ! local variables
  real(dp), parameter :: zero=0.0_dp,one=1.0_dp,two=2.0_dp, & 
       three=3.0_dp,four=4.0_dp
  real(dp) ::  apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3
  real(dp), external :: gammaf
  integer(i4b) ::  i

  apb   = alpha+beta
  f3 = zero
  if (n == 0) then
     endw2 = zero
     return
  endif
  f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
     endw2 = f1
     return
  endif
  fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
     endw2 = f2
     return
  endif
  do i=3,n
     di   = real(i-1)
     abn  = alpha+beta+di
     abnn = abn+di
     a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
     a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
     a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
     f3   =  -(a2*f2+a1*f1)/a3
     f1   = f2
     f2   = f3
  enddo
  endw2  = f3

end function endw2

!
!=======================================================================
!

function gammaf (x)
  use nrtype
  implicit none

  ! function value
  real(dp) :: gammaf
  
  
  
  ! function argument
  real(dp), intent(in) :: x

  ! local variables
  real(dp), parameter :: half=0.5_dp,one=1.0_dp,two=2.0_dp

  gammaf = one

  if (x == -half) gammaf = -two*sqrt(pi_d)
  if (x ==  half) gammaf =  sqrt(pi_d)
  if (x ==  one ) gammaf =  one
  if (x ==  two ) gammaf =  one
  if (x ==  1.5_dp) gammaf =  sqrt(pi_d)/2.0_dp
  if (x ==  2.5_dp) gammaf =  1.5_dp*sqrt(pi_d)/2.0_dp
  if (x ==  3.5_dp) gammaf =  2.5_dp*1.5_dp*sqrt(pi_d)/2.0_dp
  if (x ==  3.0_dp ) gammaf =  2._dp
  if (x ==  4.0_dp ) gammaf = 6.0_dp
  if (x ==  5.0_dp ) gammaf = 24.0_dp
  if (x ==  6.0_dp ) gammaf = 120.0_dp

  end function gammaf

!
!=====================================================================
!

  subroutine jacg (xjac,np,alpha,beta)

!=======================================================================!
!                                                                       !
! computes np Gauss points, which are the zeros of the                  !
! Jacobi polynomial with parameters alpha and beta                      !
!                                                                       !
!                  .alpha = beta =  0.0  ->  Legendre points            !
!                  .alpha = beta = -0.5  ->  Chebyshev points           !
!                                                                       !
!=======================================================================!

    use nrtype
    implicit none

    ! input output variables
    integer(i4b), intent(in) ::  np
    real(dp), intent(in) :: alpha,beta
    real(dp), dimension(np),  intent(out) ::  xjac
    
    ! local variables
    integer(i4b) ::  k,j,i,jmin,jm,n
    real(dp) ::  xlast,dth,x,x1,x2,recsum,delx,xmin,swap
    real(dp) ::  p,pd,pm1,pdm1,pm2,pdm2
    
    integer(i4b) , parameter :: K_MAX_ITER = 10
    real(dp), parameter :: zero = 0.0_dp, eps = 1.0e-12_dp
    
    pm1 = zero
    pm2 = zero
    pdm1 = zero
    pdm2 = zero
    
    xlast = 0.0_dp
    n   = np-1
    dth = 4.0_dp*atan(1.0_dp)/(2.0_dp*real(n)+2.0_dp)
    p = 0.0_dp
    pd = 0.0_dp
    jmin = 0
    do j=1,np
       if(j == 1) then
          x = cos((2.0_dp*(real(j)-1.0_dp)+1.0_dp)*dth)
       else
          x1 = cos((2.0_dp*(real(j)-1.0_dp)+1.0_dp)*dth)
          x2 = xlast
          x  = 0.5_dp*(x1+x2)
       endif
       do k=1,K_MAX_ITER
          call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
          recsum = 0.0_dp
          jm = j-1
          do i=1,jm
             recsum = recsum+1.0_dp/(x-xjac(np-i+1))
          enddo
          delx = -p/(pd-recsum*p)
          x    = x+delx
          if(abs(delx) < eps) exit
       enddo
       xjac(np-j+1) = x
       xlast        = x
    enddo
    do i=1,np
       xmin = 2.0_dp
       do j=i,np
          if(xjac(j) < xmin) then
             xmin = xjac(j)
             jmin = j
          endif
       enddo
       if(jmin /= i) then
          swap = xjac(i)
          xjac(i) = xjac(jmin)
          xjac(jmin) = swap
       endif
    enddo
    
  end subroutine jacg

!
!=====================================================================
!

  subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,bet,x)

!=======================================================================!
!                                                                       !
! Computes the Jacobi polynomial of degree n and its derivative at x    !
!                                                                       !
!=======================================================================!
    use nrtype
    implicit none

    ! input output variables
    integer(i4b), intent(in) ::  n
    real(dp), intent(in) :: x,alp,bet
    real(dp), intent(out) ::  poly,pder,polym1,pderm1,polym2,pderm2

    

    ! local variables
    real(dp) ::  apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave
    integer(i4b) ::  k

    apb  = alp+bet
    poly = 1.0_dp
    pder = 0.0_dp
    psave = 0.0_dp
    pdsave = 0.0_dp

    if (n == 0) return
    
    polyl = poly
    pderl = pder
    poly  = 0.5_dp*(alp-bet+(apb+2.0_dp)*x)
    pder  = 0.5_dp*(apb+2.d0)
    if (n == 1) return

    do k=2,n
       dk = real(k)
       a1 = 2.0_dp*dk*(dk+apb)*(2.0_dp*dk+apb-2.0_dp)
       a2 = (2.0_dp*dk+apb-1.0_dp)*(alp**2-bet**2)
       b3 = (2.0_dp*dk+apb-2.0_dp)
       a3 = b3*(b3+1.0_dp)*(b3+2.0_dp)
       a4 = 2.0_dp*(dk+alp-1.0_dp)*(dk+bet-1.0_dp)*(2.0_dp*dk+apb)
       polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
       pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
       psave  = polyl
       pdsave = pderl
       polyl  = poly
       poly   = polyn
       pderl  = pder
       pder   = pdern
    enddo
    
    polym1 = polyl
    pderm1 = pderl
    polym2 = psave
    pderm2 = pdsave
    
  end subroutine jacobf

!
!------------------------------------------------------------------------
!

  FUNCTION PNDLEG (Z,N)
    
    !------------------------------------------------------------------------!
    !                                                                        !
    !     Compute the derivative of the Nth order Legendre polynomial at Z.  !
    !     Based on the recursion formula for the Legendre polynomials.       !
    !                                                                        !
    !------------------------------------------------------------------------!
    use nrtype
    implicit none

    real(dp) :: PNDLEG
    real(dp), intent(in) :: z
    integer(i4b), intent(in) ::  n

    real(dp) ::  P1,P2,P1D,P2D,P3D,FK,P3
    integer(i4b) ::  k
    
    P1   = 1.0_dp
    P2   = Z
    P1D  = 0.0_dp
    P2D  = 1.0_dp
    P3D  = 1.0_dp
    
    do K = 1, N-1
       FK  = real(K)
       P3  = ((2.0_dp*FK+1.0_dp)*Z*P2 - FK*P1)/(FK+1.0_dp)
       P3D = ((2.0_dp*FK+1.0_dp)*P2 + (2.0_dp*FK+1.0_dp)*Z*P2D - FK*P1D) / (FK+1.0_dp)
       P1  = P2
       P2  = P3
       P1D = P2D
       P2D = P3D
    end do

    PNDLEG = P3D

  end function pndleg

!
!------------------------------------------------------------------------
!

  FUNCTION PNLEG (Z,N)

!------------------------------------------------------------------------!
!                                                                        !
!     Compute the value of the Nth order Legendre polynomial at Z.       !
!     Based on the recursion formula for the Legendre polynomials.       ! 
!                                                                        !
!------------------------------------------------------------------------!
    use nrtype
    implicit none

    real(dp) :: PNLEG
    real(dp), intent(in) :: z
    integer(i4b), intent(in) ::  n
    
    real(dp) ::  P1,P2,P3,FK
    integer(i4b) ::  k

    P1   = 1.0_dp
    P2   = Z
    P3   = P2

    do K = 1, N-1
       FK  = real(K)
       P3  = ((2.0_dp*FK+1.0_dp)*Z*P2 - FK*P1)/(FK+1.0_dp)
       P1  = P2
       P2  = P3
    enddo
    
    PNLEG = P3
    
  end function pnleg

!
!------------------------------------------------------------------------
!

  function pnormj (n,alpha,beta)
    
    use nrtype
    implicit none

    real(dp) :: pnormj

    real(dp), intent(in) ::  alpha,beta
    integer(i4b), intent(in) :: n

    real(dp) ::  one,two,dn,const,prod,dindx,frac
    real(dp), external :: gammaf
    integer(i4b) ::  i

    one   = 1.0_dp
    two   = 2.0_dp
    dn    = real(n)
    const = alpha+beta+one

    if (n <= 1) then
       prod   = gammaf(dn+alpha)*gammaf(dn+beta)
       prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
       pnormj = prod * two**const/(two*dn+const)
       return
    endif

    prod  = gammaf(alpha+one)*gammaf(beta+one)
    prod  = prod/(two*(one+const)*gammaf(const+one))
    prod  = prod*(one+alpha)*(two+alpha)
    prod  = prod*(one+beta)*(two+beta)

    do i=3,n
       dindx = real(i)
       frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
       prod  = prod*frac
    enddo

    pnormj = prod * two**const/(two*dn+const)
    
  end function pnormj

!
!------------------------------------------------------------------------
!

  subroutine zwgjd(z,w,np,alpha,beta)

!=======================================================================!
!                                                                       !
!     Z w g j d : Generate np Gauss-Jacobi points and weights           !
!                 associated with Jacobi polynomial of degree n = np-1  !
!                                                                       !
!     Note : Coefficients alpha and beta must be greater than -1.       !
!     ----                                                              !
!=======================================================================!
    
    use nrtype
    implicit none
    
    real(dp), parameter :: zero=0.0_dp,one=1.0_dp,two=2.0_dp
    
    integer(i4b), intent(in) ::  np
    real(dp), dimension(np), intent(out) ::  z,w
    real(dp), intent(in) ::  alpha,beta

    integer(i4b) ::  n,np1,np2,i
    real(dp) ::  p,pd,pm1,pdm1,pm2,pdm2
    real(dp) ::  apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef
    real(dp),  external :: gammaf,pnormj

    pd = zero
    pm1 = zero
    pm2 = zero
    pdm1 = zero
    pdm2 = zero
    
    n    = np-1
    apb  = alpha+beta
    p    = zero
    pdm1 = zero
    
    if (np <= 0) stop 'minimum number of Gauss points is 1'
    
    if ((alpha <= -one) .or. (beta <= -one)) stop 'alpha and beta must be greater than -1'
    
    if (np == 1) then
       z(1) = (beta-alpha)/(apb+two)
       w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
       return
    endif
    
    call jacg(z,np,alpha,beta)

    np1   = n+1
    np2   = n+2
    dnp1  = real(np1)
    dnp2  = real(np2)
    fac1  = dnp1+alpha+beta+one
    fac2  = fac1+dnp1
    fac3  = fac2+one
    fnorm = pnormj(np1,alpha,beta)
    rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
    do i=1,np
       call jacobf(p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
       w(i) = -rcoef/(p*pdm1)
    enddo
    
  end subroutine zwgjd

!
!------------------------------------------------------------------------
!

  subroutine zwgljd(z,w,np,alpha,beta)

!========================================================================!
!                                                                        !
!     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the      !
!     -----------   weights associated with Jacobi polynomials of degree !
!                   n = np-1.                                            !
!                                                                        !
!     Note : alpha and beta coefficients must be greater than -1.        !
!            Legendre polynomials are special case of Jacobi polynomials !
!            just by setting alpha and beta to 0.                        !
!                                                                        !
!========================================================================!

    use nrtype
    implicit none
  

    real(dp), parameter :: zero=0.0_dp,one=1.0_dp,two=2.0_dp
    
    integer(i4b), intent(in) ::  np
    real(dp), intent(in) ::  alpha,beta
    real(dp), dimension(np), intent(out) ::  z,w

    integer(i4b) :: n,nm1,i
    real(dp) :: p,pd,pm1,pdm1,pm2,pdm2
    real(dp) :: alpg,betg
    real(dp), external :: endw1,endw2
    
    p = zero
    pm1 = zero
    pm2 = zero
    pdm1 = zero
    pdm2 = zero

    n   = np-1
    nm1 = n-1
    pd  = zero
    
    if (np <= 1) stop 'minimum number of Gauss-Lobatto points is 2'
    
    ! with spectral elements, use at least 3 points
    if (np <= 2) stop 'minimum number of Gauss-Lobatto points for the SEM is  3 or more'
    
    if ((alpha <= -one) .or. (beta <= -one)) stop 'alpha and beta must be greater than -1'

    if (nm1 > 0) then
       alpg  = alpha+one
       betg  = beta+one
       call zwgjd(z(2),w(2),nm1,alpg,betg)
    endif
    
    z(1)  = - one
    z(np) =  one
    
    do i=2,np-1
       w(i) = w(i)/(one-z(i)**2)
    enddo
    
    call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
    w(1)  = endw1(n,alpha,beta)/(two*pd)
    call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
    w(np) = endw2(n,alpha,beta)/(two*pd)
    
  end subroutine zwgljd

