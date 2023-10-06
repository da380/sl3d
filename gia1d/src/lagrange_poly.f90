
subroutine lagrange_any(xi,NGLL,xigll,h,hprime)
  !-----------------------------------------------------------------------!
  !                                                                       !
  ! subroutine to compute the Lagrange interpolants based upon the        !
  ! GLL points and their first derivatives at any point xi in [-1,1]      !
  !                                                                       !
  !-----------------------------------------------------------------------!

    
    use nrtype
    implicit none

    ! input and output variables
    integer(i4b), intent(in) :: NGLL
    real(dp), intent(in) :: xi
    real(dp), dimension(NGLL), intent(in) :: xigll
    real(dp), dimension (NGLL), intent(out) :: h,hprime
    
    ! local variables
    integer(i4b) ::  dgr,i,j
    real(dp) ::  prod1,prod2



    do dgr=1,NGLL
       
       prod1 = 1.0_dp
       prod2 = 1.0_dp
       do i=1,NGLL
          if(i /= dgr) then
             prod1 = prod1*(xi-xigll(i))
             prod2 = prod2*(xigll(dgr)-xigll(i))
          endif
       enddo
       h(dgr)=prod1/prod2

       hprime(dgr)=0.0d0
       do i=1,NGLL
          if(i /= dgr) then
             prod1=1.0d0
             do j=1,NGLL
                if(j /= dgr .and. j /= i) prod1 = prod1*(xi-xigll(j))
             enddo
             hprime(dgr) = hprime(dgr)+prod1
          endif
       enddo
       hprime(dgr) = hprime(dgr)/prod2
       
    enddo

  end subroutine lagrange_any



  
  function lagrange_deriv_GLL(I,j,ZGLL,NZ)

    !------------------------------------------------------------------------!
    !                                                                        !
    !     Compute the value of the derivative of the I-th                    !
    !     Lagrange interpolant through the                                   !
    !     NZ Gauss-Lobatto Legendre points ZGLL at point ZGLL(j)             !
    !                                                                        !
    !------------------------------------------------------------------------!

    use nrtype
    implicit none

    ! function value
    real(dp) :: lagrange_deriv_GLL
    ! function inputs
    integer(i4b), intent(in)  ::  i,j,nz    
    real(dp) ::  zgll(0:nz-1)

    ! local variables
    integer(i4b) ::  degpoly

 !   double precision, external :: pnleg,pndleg
    real(dp), external :: pnleg,pndleg

    degpoly = nz - 1
    if (i == 0 .and. j == 0) then
       lagrange_deriv_GLL = - real(degpoly)*(real(degpoly)+1.0_dp) / 4.0_dp
    else if (i == degpoly .and. j == degpoly) then
       lagrange_deriv_GLL = real(degpoly)*(real(degpoly)+1.0_dp) / 4.0_dp
    else if (i == j) then
       lagrange_deriv_GLL = 0.0_dp
    else
       lagrange_deriv_GLL = pnleg(zgll(j),degpoly) / &
            (pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))) &
            + (1.0_dp-zgll(j)*zgll(j))*pndleg(zgll(j),degpoly) & 
            / (real(degpoly)*(real(degpoly)+1.0_dp)* & 
            pnleg(zgll(i),degpoly)*(zgll(j)-zgll(i))* & 
            (zgll(j)-zgll(i)))
    endif
    
  end function lagrange_deriv_GLL

