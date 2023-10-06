module module_gia

  use nrtype
  use module_util
  use module_mat
  use module_model
  use module_function
  use module_fourier
  use module_spline
  
  implicit none

  integer(i4b), protected :: nlat_ice5g = 180
  integer(i4b), protected :: nlong_ice5g = 360

  integer(i4b), protected :: nt_ice5g = 39
  integer(i4b), protected :: nt_ice6g = 58

  complex(dpc), dimension(:), allocatable, save :: hall

  real(dp), dimension(-18:18,2), save :: rob_coefs,crob_coefs
  real(dp), dimension(-18:18), save :: rob_lats


  type displ_point

     real(dp) :: time
     integer(i4b) :: dir
     !(1) = radial, (2) = theta, (3) = phi, (4) = gravity pert
     real(dp) :: displ

  end type displ_point


  type loc_data

     real(dp) :: lat,long
     integer(i4b) :: ntd
     type(displ_point), dimension(:), allocatable :: displs

  end type loc_data


contains

  !=============================================!
  !        Routines for calculating RHSs        !
  !=============================================!


  subroutine memdot_bvec_all(memm1,UVpm1,Wm1,dUm1,dVm1,dWm1,memdot,bUVp,bW,load)
    
    use nrtype
    implicit none


    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(in) :: dUm1,dVm1,dWm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: bUVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: bW
    complex(dpc), dimension((lmax+1)**2), intent(in), optional :: load

    integer(i4b) :: inode,ispec,ilayer,pa,ua,va,wa,l,m,ispec1,lmtot,lmcur,i,jnode,ub,vb,wb
    real(dp) :: zeta2,zeta,zeta2m2,rl,rr,rmu
    complex(dpc) :: uu,vv,ww,dU,dV,dW
    complex(dpc), dimension(5,(lmax+1)**2) :: dten
    complex(dpc), dimension(5,ngl,nphi) :: Aspat
    complex(dpc), dimension(5,(lmax+1)**2) :: Alm

    ! A(1) = A(-,-)
    ! A(2) = A(-,0)
    ! A(3) = A(-,+) = 0.5*A(0,0)
    ! A(4) = A(0,+)
    ! A(5) = A(+,+)

    memdot = 0.0_dp
    bUVp = 0.0_dp
    bW = 0.0_dp

    ! Calculate spatial d - m
    do ilayer = 1,nsect

       ! No terms in elastic or fluid regions
       if (sect_ind(3,ilayer) /=2) cycle

       Alm = 0.0_dp

       !$OMP PARALLEL DO PRIVATE(pa,ua,va,wa,ispec1,lmtot,lmcur,i,zeta2,zeta,zeta2m2), &
       !$OMP& PRIVATE(rl,rr,rmu,uu,vv,ww,dU,dV,dW,ub,vb,wb), &
       !$OMP& FIRSTPRIVATE(dten,Aspat,Alm), &
       !$OMP& SHARED(lmax,memm1,UVpm1,Wm1,dUm1,dVm1,dWm1,memdot,bUVp,bW,load)
       do ispec = sect_ind(1,ilayer), sect_ind(2,ilayer)
          Alm = 0.0_dp
          do inode = 1,ngll
             dten = 0.0_dp
             Aspat = 0.0_dp
             rr = r_node(inode,ispec)
             rmu = mui_node(1,inode,ispec)
             if ((ispec == 1) .and. (inode == 1)) rr = r_node(2,1)
             do l = 0,lmax
                lmtot = l**2 + 1
                call get_depth(l,rl,ispec1)
                if (ispec1 > ispec) then
                   exit
                end if

                ua = gvn_ssg(inode,ispec,ispec1,2)
                va = gvn_ssg(inode,ispec,ispec1,3)
                wa = gvn_tor(inode,ispec,ispec1)

                zeta2 = l*(l+1)
                zeta = sqrt(zeta2)
                zeta2m2 = zeta2 - 2.0_dp

                do m = -l,l  
                   lmcur = lmtot + l + m
                   uu = UVpm1(lmcur,ua)
                   vv = UVpm1(lmcur,va)
                   ww = Wm1(lmcur,wa)
                   dU = dUm1(lmcur,inode,ispec)
                   dV = dVm1(lmcur,inode,ispec)
                   dW = dWm1(lmcur,inode,ispec)

                   if (l >= 2) then

                      dten(1,lmcur) = 0.5_dp*zeta*sqrt(zeta2m2)*(vv - ii*ww)/rr
                      dten(2,lmcur) = zeta*(rr*dV - vv + uu + ii*(ww - rr*dW)) &
                                      /(2.0_dp*sqrt(2.0_dp)*rr)
                      dten(3,lmcur) = (rr*dU - uu + 0.5_dp*zeta2*vv)/(3.0_dp*rr)
                      dten(4,lmcur) = zeta*(rr*dV - vv + uu - ii*(ww - rr*dW)) &
                                      /(2.0_dp*sqrt(2.0_dp)*rr)
                      dten(5,lmcur) = 0.5_dp*zeta*sqrt(zeta2m2)*(vv + ii*ww)/rr
                      
                   else if (l == 1) then
                      dten(2,lmcur) = zeta*(rr*dV - vv + uu + ii*(ww - rr*dW)) &
                                      /(2.0_dp*sqrt(2.0_dp)*rr)
                      dten(3,lmcur) = (rr*dU - uu + 0.5_dp*zeta2*vv)/(3.0_dp*rr)
                      dten(4,lmcur) = zeta*(rr*dV - vv + uu - ii*(ww - rr*dW)) &
                                      /(2.0_dp*sqrt(2.0_dp)*rr)
                   else 
                      dten(3,lmcur) = (rr*dU - uu + 0.5_dp*zeta2*vv)/(3.0_dp*rr)
                   end if
                end do
             end do

             call fun_from_coefs_ten_tl(dten,Aspat)
             do i = 1,5
                Aspat(i,:,:) = si_spat_node(:,:,1,inode,ispec) &
                               *(Aspat(i,:,:) - memm1(i,:,:,1,inode,ispec))
                ! Calculate memdot
                memdot(i,:,:,1,inode,ispec) = Aspat(i,:,:)
             end do
  
             ! Calculate Alm
             call coefs_from_fun_ten_tl(Aspat,Alm(:,:))


             ! Calculate bvecs
             do l = 1,lmax
                lmtot = l**2 + 1
                call get_depth(l,rl,ispec1)
                if (ispec1 > ispec) then
                   exit
                end if

                ua = gvn_ssg(inode,ispec,ispec1,2)
                va = ua + 1
                wa = gvn_tor(inode,ispec,ispec1)

                zeta2 = l*(l+1)
                zeta = sqrt(zeta2)
                zeta2m2 = zeta2 - 2.0_dp

                rr = r_node(inode,ispec)
                rmu = mui_node(1,inode,ispec)

                do m = -l,l
                   lmcur = lmtot + l + m

                   bUVp(lmcur,ua) = bUVp(lmcur,ua) + jac_element(ispec)*wgll(inode)*rr*rmu &
                                                     *(sqrt(2.0_dp)*zeta &
                                                       *(Alm(2,lmcur)+Alm(4,lmcur)) &
                                                      - 4.0_dp*Alm(3,lmcur))

                   bUVp(lmcur,va) = bUVp(lmcur,va) + jac_element(ispec)*wgll(inode)*rr*rmu &
                                                     *(zeta*sqrt(zeta2m2) &
                                                       *(Alm(1,lmcur)+Alm(5,lmcur)) &
                                                      + 2.0_dp*zeta2*Alm(3,lmcur) &
                                                      - sqrt(2.0_dp)*zeta &
                                                        *(Alm(2,lmcur) + Alm(4,lmcur)))

                   bW(lmcur,wa) = bW(lmcur,wa) + ii*jac_element(ispec)*wgll(inode)*rr*rmu &
                                                 *(zeta*sqrt(zeta2m2) &
                                                   *(Alm(1,lmcur) - Alm(5,lmcur)) &
                                                  - sqrt(2.0_dp)*zeta &
                                                    *(Alm(2,lmcur) - Alm(4,lmcur)))

                   ! Derivative terms
                   do jnode = 1,ngll

                      ub = gvn_ssg(jnode,ispec,ispec1,2)
                      vb = ub + 1
                      wb = gvn_tor(jnode,ispec,ispec1)

                      bUVp(lmcur,ub) = bUVp(lmcur,ub) + wgll(inode)*hprime(inode,jnode)*rr*rmu &
                                                        *4.0_dp*rr*Alm(3,lmcur)
                      bUVp(lmcur,vb) = bUVp(lmcur,vb) + wgll(inode)*hprime(inode,jnode)*rr*rmu &
                                                        *sqrt(2.0_dp)*zeta*rr &
                                                        *(Alm(2,lmcur) + Alm(4,lmcur))
                      bW(lmcur,wb) = bW(lmcur,wb) + ii*wgll(inode)*hprime(inode,jnode)*rr*rmu &
                                                    *sqrt(2.0_dp)*zeta*rr &
                                                    *(Alm(2,lmcur) - Alm(4,lmcur))
                      
                   end do

                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end do

    ! Add load terms
    lmtot = 1
    if (present(load)) then
       do l = 1,lmax
          pa = pa_surf(l)
          ua = pa + 1          
          do m = -l,l
             lmcur = lmtot + l + m
             bUVp(lmcur,pa) =  bUVp(lmcur,pa)  - load(lmcur)*(r_node(ngll,nspec))**2
             bUVp(lmcur,ua) = bUVp(lmcur,ua) - load(lmcur)*grav_node(ngll,nspec)*(r_node(ngll,nspec))**2
          end do
          lmtot = lmtot + 2*l + 1
       end do
    end if


  end subroutine memdot_bvec_all


  subroutine bvec_sph_mem_l(ispec1,l,nrows,mem,UVptm1,dUm1,dVm1,b)

    use nrtype
    use module_util
    use module_model

    implicit none

    integer(i4b), intent(in) :: ispec1,l,nrows
    complex(dpc), dimension(5,-lmax:lmax,ntau,ngll,nspec), intent(in) :: mem
    complex(dpc), dimension(-lmax:lmax,nrows), intent(in) :: UVptm1
    complex(dpc), dimension(-lmax:lmax,ngll,nspec), intent(in) :: dUm1,dVm1
    complex(dpc), dimension(-lmax:lmax,nrows), intent(out) :: b
    integer(i4b) :: inode,ispec,iglob,jnode, &
                    jglob,ia,ib,ja,jb,i,ilayer,ispec11,m
    real(dp) :: rr,tmp,zeta2,zeta2m2,zeta,c      
    real(dp), dimension(ntau) :: rmui,rsi

    !open(newunit=io,file='togtest.out')
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zeta = sqrt(zeta2)
    c = sqrt((2*l+1)/(fourpi_d))
    
    b = 0.0_dp

    ia = 1
    do ilayer = 1,nsect
       if (ispec1 > sect_ind(2,ilayer)) cycle

       ! No terms in elastic or fluid regions
       if (sect_ind(3,ilayer) /= 2) cycle

       ispec11 = max(ispec1,sect_ind(1,ilayer))

       ! begin loop over the spectral elements
       do ispec = ispec11,sect_ind(2,ilayer)

          do inode = 1,ngll
             
             iglob = ibool(inode,ispec)
             rr = r_node(inode,ispec)
             rmui(:) = mui_node(:,inode,ispec)
             rsi(:) = si_node(:,inode,ispec)
  
             ! non-derivative terms
             if (inode == 1) then
                ia = gvn_ssg(inode,ispec,ispec1,2)
             else
                ia = ia + 3
             end if
             ib = ia + 1

             ! memory terms
             ! U' terms

             do i = 1,ntau

                do m = -l,l

                   !if (m > l) exit
                   !if (abs(m) > l) cycle
             
                   b(m,ia) = b(m,ia) + zeta2*jac_element(ispec)*wgll(inode) &
                                       *(rmui(i)*rsi(i) &
                                       *(rr*dVm1(m,inode,ispec) &
                                       - UVptm1(m,ib) + UVptm1(m,ia) &
                                       - mem(4,m,i,inode,ispec)))


                   ! V' terms          
                   b(m,ib) = b(m,ib) + (2.0_dp*zeta2/3.0_dp) &
                                       *jac_element(ispec)*wgll(inode) &
                                       *(rmui(i)*rsi(i) &
                                       *(rr*dUm1(m,inode,ispec) &
                                       - UVptm1(m,ia) &
                                       + 0.5_dp*zeta2*UVptm1(m,ib) &
                                       - mem(3,m,i,inode,ispec)))

                   b(m,ib) = b(m,ib) + zeta2*zeta2m2 &
                                       *jac_element(ispec)*wgll(inode) &
                                       *(rmui(i)*rsi(i) &
                                       *(UVptm1(m,ib) &
                                       - mem(1,m,i,inode,ispec)))

                end do

             end do

             ! derivative terms
             
             do jnode = 1,ngll

                if (jnode == 1) then
                   ja = gvn_ssg(jnode,ispec,ispec1,2)
                else
                   ja = ja + 3
                end if
                jb = ja + 1

                tmp = hprime(jnode,inode)/jac_element(ispec)
                jglob = ibool(jnode,ispec)
                rr = r_node(jnode,ispec)
                rmui(:) = mui_node(:,jnode,ispec)
                rsi(:) = si_node(:,jnode,ispec)

                ! memory terms
                tmp = rr*tmp
                if (jnode == inode) tmp = tmp - 1.0_dp

                do i = 1,ntau

                   do m = -l,l
                      !if (m > l) exit
                      !if (abs(m) > l) cycle

                      ! U' terms
                      b(m,ia) = b(m,ia) + (4.0_dp/3.0_dp)*tmp*wgll(jnode) &
                                          *jac_element(ispec) &
                                          *(rmui(i)*rsi(i) &
                                          *(rr*dUm1(m,jnode,ispec) &
                                          - UVptm1(m,ja) &
                                          + 0.5_dp*zeta2*UVptm1(m,jb) &
                                          - mem(3,m,i,jnode,ispec)))


                      ! V' terms
                      b(m,ib) = b(m,ib) + zeta2*jac_element(ispec) &
                                          *wgll(jnode)*tmp &
                                          *(rmui(i)*rsi(i) &
                                          *(rr*dVm1(m,jnode,ispec) &
                                          - UVptm1(m,jb) &
                                          + UVptm1(m,ja) &
                                          - mem(4,m,i,jnode,ispec)))

                   end do

                end do
             
             end do

          end do

       end do

    end do

  end subroutine bvec_sph_mem_l



  !=============================================!
  !      Routines for forward calculations      !
  !=============================================!

 
  subroutine time_derivs_full(memm1,UVpm1,Wm1,memdot,UVpdot,Wdot,load_coefs)
  ! routine to calculate velocities from displacements

    implicit none

    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVpdot
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: Wdot
    complex(dpc), dimension((lmax+1)**2), intent(in), optional :: load_coefs


    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dUm1,dVm1,dWm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: bUVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: bW
    integer(i4b) :: kl,ku,ldab,nrows,info,m,kd,ispec1, &
                    inode,ispec,ua,va,i,ilayer,ispec11,l
    real(dp) :: zeta,zeta2,rr,rl
    real(dp), dimension(ntau) :: rsi
    complex(dpc), dimension(:,:), allocatable :: UVp_trans,W_trans

    UVpdot = 0.0_dp
    Wdot = 0.0_dp
    memdot = 0.0_dp

    !-------------------------------------!
    !         Calculate stresses          !
    !-------------------------------------!

    call calculate_dUVW_full(UVpm1,Wm1,dUm1,dVm1,dWm1)


    !-------------------------------------!
    !      Calculate memdot and RHSs      !
    !-------------------------------------!

    call memdot_bvec_all(memm1,UVpm1,Wm1,dUm1,dVm1,dWm1,memdot,bUVp,bW,load_coefs)


    allocate(UVp_trans(nglob_ssg,-lmax:lmax),W_trans(nglob_tor,-lmax:lmax))

    UVp_trans = 0.0_dp
    W_trans = 0.0_dp


    do l = 1,lmax

       zeta2 = l*(l+1)
       zeta = sqrt(zeta2)

       call get_depth(l,rl,ispec1)


       !---------------------------------------!
       !      Calculate spheroidal modes       !
       !---------------------------------------!
     
       kl = 3*(ngll-1) + 2
       ku = 3*(ngll-1) + 2
       kd = ku
    
       nrows = gvn_ssg(ngll,nspec,ispec1,3)

       UVp_trans(1:nrows,-l:l) = transpose(bUVp(l**2 + 1:(l+1)**2,1:nrows))

       if (l == 1) then
          ldab = 2*kl + ku + 1
          call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                      ipiv_ssg_1(1:nrows),UVp_trans(1:nrows,-l:l),nrows,info)
       
       else
          ldab = kd + 1             
          call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                      UVp_trans(1:nrows,-l:l),nrows,info)
       end if

       !print *, l,UVp_trans(nrows-1,0)

       UVpdot(l**2 + 1:(l+1)**2,1:nrows) = transpose(UVp_trans(1:nrows,-l:l))   

       !---------------------------------------!
       !       Calculate toroidal modes        !
       !---------------------------------------!
     
       kl = ngll-1
       ku = ngll-1
       kd = ku
    
       nrows = gvn_tor(ngll,nspec,ispec1)

       W_trans(1:nrows,-l:l) = transpose(bW(l**2 + 1:(l+1)**2,1:nrows))

       if (l == 1) then
          ldab = 2*kl + ku + 1
          call zgbtrs('N',nrows,kl,ku,2*l+1,aa_tor_1(:,1:nrows),ldab, &
                      ipiv_tor_1(1:nrows),W_trans(1:nrows,-l:l),nrows,info)
       
       else
          ldab = kd + 1             
          call zpbtrs('U',nrows,kd,2*l+1,aa_tor(:,1:nrows,l),ldab, &
                      W_trans(1:nrows,-l:l),nrows,info)
       end if


       Wdot(l**2 + 1:(l+1)**2,1:nrows) = transpose(W_trans(1:nrows,-l:l))   

       UVp_trans = 0.0_dp
       W_trans = 0.0_dp

    end do
    

  end subroutine time_derivs_full



  subroutine update_full_rk(dt,memm1,UVpm1,Wm1,mem,UVp,W,load_coefs)
  ! routine to calculate the displacements from the previous displacements
  ! using second-order Runge-Kutta time-stepping  

    use nrtype

    implicit none

    real(dp), intent(in) :: dt
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: mem
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: W
    complex(dpc), dimension((lmax+1)**2), intent(in), optional :: load_coefs

    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memtmp
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wtmp


    memdot = 0.0_dp
    UVpdot = 0.0_dp
    Wdot = 0.0_dp

    UVptmp = 0.0_dp
    Wtmp = 0.0_dp
    memtmp = 0.0_dp
    
    if (present(load_coefs)) then   
       call time_derivs_full(memm1,UVpm1,Wm1,memdot,UVpdot,Wdot,load_coefs)


       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       Wtmp = Wm1 + Wdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp


       call time_derivs_full(memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot,load_coefs)


       UVp = UVpm1 + UVpdot*dt
       W = Wm1 + Wdot*dt
       mem = memm1 + memdot*dt

    else

       call time_derivs_full(memm1,UVpm1,Wm1,memdot,UVpdot,Wdot)


       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       Wtmp = Wm1 + Wdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp


       call time_derivs_full(memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot)


       UVp = UVpm1 + UVpdot*dt
       W = Wm1 + Wdot*dt
       mem = memm1 + memdot*dt

    end if

  
  end subroutine update_full_rk


  !=========================================!
  !            Useful routines              !
  !=========================================!

  subroutine get_lat_long(ilat,ilong,lat,long)

    implicit none

    integer(i4b), intent(in) :: ilat,ilong
    real(dp), intent(out) :: lat,long

    lat = ilat - 90.5_dp
    long = ilong - 181.0_dp

  end subroutine get_lat_long


  subroutine get_ilat_ilong(lat,long,ilat1,ilong1)

    implicit none
    
    real(dp), intent(in) :: lat,long
    integer(i4b), intent(out) :: ilat1,ilong1

    integer(i4b) :: ilat,ilong
    real(dp) :: lati,longi

    ilat1 = floor(lat + 90.5_dp)
    ilong1 = floor(long + 181.0_dp) 

    if (ilong1 == 361) ilong1 = 1
    
  end subroutine get_ilat_ilong

    
  function interp(func_arr,lat,long)

    implicit none

    real(dp) :: interp
    real(dp), dimension(nlong_ice5g,nlat_ice5g), intent(in) :: func_arr
    real(dp), intent(in) :: lat,long

    integer(i4b) :: ilat1,ilong1,ilat2,ilong2
    real(dp) :: xx,yy,f11,f12,f21,f22

    xx = long - floor(long)
    yy = lat - floor(lat)

    !if (xx > 0.5) then
    !   xx = xx - 0.5
    !else
    !   xx = xx + 0.5
    !end if

    !if (yy > 0.5) then
    !   yy = yy - 0.5
    !else
    !   yy = yy + 0.5
    !end if

    !ilong1 = long - xx

    call get_ilat_ilong(lat,long,ilat1,ilong1)

    ilat2 = ilat1 + 1
    ilong2 = ilong1 + 1
    if (ilat2 == 181) ilat2 = 180
    if (ilong2 == 361) ilong2 = 1
    if (ilat1 == 0) ilat1 = 1
      
    f11 = func_arr(ilong1,ilat1)
    f12 = func_arr(ilong1,ilat2)
    f21 = func_arr(ilong2,ilat1)
    f22 = func_arr(ilong2,ilat2)

    interp = f11 + xx*(f21 - f11) + yy*(f12 - f11) + xx*yy*(f11 + f22 - f12 - f21)

  end function interp

  
  function time_point(it,dt)
    ! finds time point of ice_5g which is (just) below time point of calculation

    implicit none

    integer(i4b) :: time_point
    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt

    real(dp) :: t_act

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       time_point = floor(t_act/(1000.0_dp*yr2sec/t_norm)) + 1
    else
       time_point = nt_ice5g - ceiling(((21000.0_dp*yr2sec/t_norm) - t_act)/(500.0_dp*yr2sec/t_norm))
    end if

    !print *, it,time_point,t_act*t_norm*sec2yr

  end function time_point



  subroutine time_point_interp(it,dt,ice,ice_it,icedot_it)

    implicit none

    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice
    complex(dpc), dimension(ngl,nphi), intent(out) :: ice_it
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: icedot_it

    integer(i4b) :: t_point
    real(dp) :: t_act,t_rem

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       t_point = floor(t_act/(1000.0_dp*yr2sec/t_norm)) + 1
       t_rem = t_act - (t_point - 1)*(1000_dp*yr2sec/t_norm)
       if (t_rem < 0.0_dp) t_rem = 0.0_dp
       if (t_rem > 999.5_dp*yr2sec/t_norm) then
          t_point = t_point + 1
          t_rem = 0.0_dp
       end if
          
       ice_it = ice(:,:,t_point) + (t_norm*t_rem/(1000.0_dp*yr2sec)) &
                                   *(ice(:,:,t_point+1) - ice(:,:,t_point))
       
       if (present(icedot_it)) then
          icedot_it = (t_norm/(1000.0_dp*yr2sec))*(ice(:,:,t_point+1) - ice(:,:,t_point))
       end if

       
    else
       t_point = nt_ice5g - ceiling(((21000.0_dp*yr2sec/t_norm) - t_act)/(500.0_dp*yr2sec/t_norm))
       if (t_point == nt_ice5g) then
          ice_it = ice(:,:,nt_ice5g)
          if (present(icedot_it)) then
             icedot_it = (t_norm/(500.0_dp*yr2sec))*(ice(:,:,t_point) - ice(:,:,t_point-1))
          end if
       else
          t_rem = t_act - 4000.0_dp*yr2sec/t_norm - (t_point - 5)*(500.0_dp*yr2sec/t_norm)
          if (t_rem < 0.0_dp) t_rem = 0.0_dp
          if (t_rem > 499.5_dp*yr2sec/t_norm) then
             t_point = t_point + 1
             t_rem = 0.0_dp
          end if
          ice_it = ice(:,:,t_point) + (t_norm*t_rem/(500.0_dp*yr2sec)) &
                                      *(ice(:,:,t_point+1) - ice(:,:,t_point))
          if (present(icedot_it)) then
             icedot_it = (t_norm/(500.0_dp*yr2sec))*(ice(:,:,t_point+1) - ice(:,:,t_point))
          end if
       end if
    end if

    !print *, t_act*t_norm*sec2yr,t_point,t_rem*t_norm*sec2yr

  end subroutine time_point_interp




  subroutine load_interp(it,dt,ice,ice_it)

    implicit none

    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,2), intent(in) :: ice
    complex(dpc), dimension(ngl,nphi), intent(out) :: ice_it

    real(dp) :: t_act,t_diff

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       t_diff = t_act - (1000.0_dp*yr2sec/t_norm)*floor(t_act*t_norm/(1000.0_dp*yr2sec))
       ice_it = ice(:,:,1) + (t_norm*t_diff/(1000.0_dp*yr2sec))*(ice(:,:,2) - ice(:,:,1))
    else
       t_diff = t_act - (500.0_dp*yr2sec)*floor(t_act*t_norm/(500.0_dp*yr2sec))/t_norm
       ice_it = ice(:,:,1) + (t_norm*t_diff/(500.0_dp*yr2sec))*(ice(:,:,2) - ice(:,:,1))
    end if

    print *, it,t_diff*t_norm*sec2yr

  end subroutine load_interp


  subroutine dload_interp(it,dt,ice,icedot_it)

    implicit none

    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(:,:), intent(in) :: ice
    complex(dpc), dimension(:), intent(out) :: icedot_it

    real(dp) :: t_act

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       icedot_it = (t_norm/(1000.0_dp*yr2sec))*(ice(:,2) - ice(:,1))
    else
       icedot_it = (t_norm/(500.0_dp*yr2sec))*(ice(:,2) - ice(:,1))
    end if

  end subroutine dload_interp


  subroutine dloadt_interp(it,dt,ice,icedot_it)

    implicit none

    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,2), intent(in) :: ice
    complex(dpc), dimension(ngl,nphi), intent(out) :: icedot_it

    real(dp) :: t_act

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       icedot_it = (t_norm/(1000.0_dp*yr2sec))*(ice(:,:,2) - ice(:,:,1))
    else
       icedot_it = (t_norm/(500.0_dp*yr2sec))*(ice(:,:,2) - ice(:,:,1))
    end if

  end subroutine dloadt_interp


  subroutine time_point_rem(it,dt)
    ! finds time point of ice_5g which is (just) below time point of calculation

    implicit none

    integer(i4b) :: time_point
    integer(i4b), intent(in) :: it
    real(dp), intent(in) :: dt

    real(dp) :: t_act

    t_act = (it-1)*dt

    if (t_act < 4000*yr2sec/t_norm) then
       time_point = floor(t_act/(1000.0_dp*yr2sec/t_norm)) + 1
    else
       time_point = nt_ice5g - ceiling(((21000.0_dp*yr2sec/t_norm) - t_act)/(500.0_dp*yr2sec/t_norm))
    end if

  end subroutine time_point_rem


  subroutine calculate_dUV(ispec1,l,UVp,dU,dV)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: ispec1,l
    complex(dpc), dimension(-l:l,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension(-l:l,ngll,nspec), intent(out) :: dU,dV

    integer(i4b) :: ispec,inode,jnode,ilayer,ispec11,m,ua,va

    real(dp) :: hp,jac

    dU = 0.0_dp
    dV = 0.0_dp

    do ilayer = 1,nsect
       if (ispec1 > sect_ind(2,ilayer)) cycle

       ! No terms in fluid regions
       if (sect_ind(3,ilayer) == 3) cycle

       ispec11 = max(ispec1,sect_ind(1,ilayer))

       ! begin loop over the spectral elements
       do ispec = ispec11,sect_ind(2,ilayer)
          do inode = 1,ngll
             do jnode = 1,ngll
                hp = hprime(inode,jnode)
                jac = jac_element(ispec)

                ! non-derivative terms
                if (jnode == 1) then
                   ua = gvn_ssg(jnode,ispec,ispec1,2)
                else
                   ua = ua + 3
                end if
                va = ua + 1

                do m = -l,l

                   dU(m,inode,ispec) = dU(m,inode,ispec) &
                                       + UVp(m,ua)*hp/jac

                   dV(m,inode,ispec) = dV(m,inode,ispec) &
                                       + UVp(m,va)*hp/jac

                end do


             end do
          end do
       end do
    end do


  end subroutine calculate_dUV


  subroutine calculate_dUV_all(UVp,dU,dV)

    use nrtype

    implicit none

    complex(dpc), dimension(-lmax:lmax,nglob_ssg,lmax), intent(in) :: UVp
    complex(dpc), dimension(-lmax:lmax,ngll,nspec,lmax), intent(out) :: dU,dV

    integer(i4b) :: ispec,inode,jnode,ilayer,ispec11,m,ua,va,l,ispec1,lmcur

    real(dp) :: hp,jac,rl

    dU = 0.0_dp
    dV = 0.0_dp
    
    lmcur = 2
    do l = 1,lmax
       call get_depth(l,rl,ispec1)       

       do ilayer = 1,nsect
          if (ispec1 > sect_ind(2,ilayer)) cycle
          
          ! No terms in fluid regions
          if (sect_ind(3,ilayer) == 3) cycle

          ispec11 = max(ispec1,sect_ind(1,ilayer))

          ! begin loop over the spectral elements
          do ispec = ispec11,sect_ind(2,ilayer)
             do inode = 1,ngll
                do jnode = 1,ngll
                   hp = hprime(inode,jnode)
                   jac = jac_element(ispec)

                   ! non-derivative terms
                   if (jnode == 1) then
                      ua = gvn_ssg(jnode,ispec,ispec1,2)
                   else
                      ua = ua + 3
                   end if
                   va = ua + 1

                   do m = -lmax,lmax
                      if(m > l) exit
                      if(abs(m) > l) cycle

                      dU(m,inode,ispec,l) = dU(m,inode,ispec,l) &
                                            + UVp(m,ua,l)*hp/jac

                      dV(m,inode,ispec,l) = dV(m,inode,ispec,l) &
                                            + UVp(m,va,l)*hp/jac

                   end do

                end do


             end do
          end do
       end do
       lmcur = lmcur + 2*l + 1
    end do


  end subroutine calculate_dUV_all


  subroutine calculate_dUVW_full(UVp,W,dU,dV,dW)

    use nrtype

    implicit none

    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: W
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(out) :: dU,dV,dW

    integer(i4b) :: ispec,inode,jnode,ilayer,ispec11,l,m,ua,va,wa,lmcur,ispec1
    real(dp) :: hp,jac,rl

    dU = 0.0_dp
    dV = 0.0_dp
    dW = 0.0_dp

    do l = 0,lmax
       call get_depth(l,rl,ispec1)

       do ilayer = 1,nsect
          if (ispec1 > sect_ind(2,ilayer)) cycle

          ! No terms in fluid regions
          if (sect_ind(3,ilayer) == 3) cycle

          ispec11 = max(ispec1,sect_ind(1,ilayer))

          ! begin loop over the spectral elements
          do ispec = ispec11,sect_ind(2,ilayer)
             do inode = 1,ngll
                do jnode = 1,ngll
                   hp = hprime(inode,jnode)
                   jac = jac_element(ispec)

                   ! non-derivative terms
                   if (jnode == 1) then
                      ua = gvn_ssg(jnode,ispec,ispec1,2)
                      wa = gvn_tor(jnode,ispec,ispec1)
                   else
                      ua = ua + 3
                      wa = wa + 1
                   end if
                   va = ua + 1

                   do m = -l,l
                      lmcur = l**2 + 1 + l + m

                      dU(lmcur,inode,ispec) = dU(lmcur,inode,ispec) &
                                              + UVp(lmcur,ua)*hp/jac

                      dV(lmcur,inode,ispec) = dV(lmcur,inode,ispec) &
                                              + UVp(lmcur,va)*hp/jac

                      dW(lmcur,inode,ispec) = dW(lmcur,inode,ispec) &
                                              + W(lmcur,wa)*hp/jac

                   end do

                end do


             end do
          end do
       end do
    end do


  end subroutine calculate_dUVW_full



  subroutine calculate_dU0(Up0,dU)

    use nrtype

    implicit none

    complex(dpc), dimension(:), intent(in) :: Up0
    complex(dpc), dimension(ngll,nspec), intent(out) :: dU

    integer(i4b) :: ispec,inode,jnode,ilayer,ispec11,m,ua,va,l,lmcur

    real(dp) :: hp,jac,rl

    dU = 0.0_dp
    
    do ilayer = 1,nsect
          
       ispec11 = max(1,sect_ind(1,ilayer))

       ! begin loop over the spectral elements
       do ispec = ispec11,sect_ind(2,ilayer)
          do inode = 1,ngll
             do jnode = 1,ngll
                hp = hprime(inode,jnode)
                jac = jac_element(ispec)

                ! non-derivative terms
                if (jnode == 1) then
                   ua = gvn_ssg_0(jnode,ispec,2)
                else
                   ua = ua + 2
                end if

                dU(inode,ispec) = dU(inode,ispec) + Up0(ua)*hp/jac

             end do

          end do

       end do
    end do

  end subroutine calculate_dU0


  subroutine surface_fields(UVp,ps,Us)

    implicit none

    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2), intent(out) :: ps,Us

    integer(i4b) :: l,ispec1,ua,pa,lmcur,m
    real(dp) :: rl

    ps = 0.0_dp
    Us = 0.0_dp

    lmcur = 2
    do l = 1,lmax

       call get_depth(l,rl,ispec1)

       pa = pa_surf(l)
       ua = pa + 1

       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle
          
          ps(lmcur+l+m) = UVp(lmcur+l+m,pa)
          Us(lmcur+l+m) = UVp(lmcur+l+m,ua)

       end do          

       lmcur = lmcur + 2*l + 1
    end do

  end subroutine surface_fields

  subroutine surface_fields_full(UVp,ps,Us)

    implicit none

    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2), intent(out) :: ps,Us

    integer(i4b) :: l,ispec1,ua,pa,lmcur,m
    real(dp) :: rl

    ps = 0.0_dp
    Us = 0.0_dp

    lmcur = 1
    do l = 0,lmax
       call get_depth(l,rl,ispec1)
       pa = gvn_ssg(ngll,nspec,ispec1,1)
       ua = pa + 1

       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle
          
          ps(lmcur+l+m) = UVp(lmcur+l+m,pa)
          Us(lmcur+l+m) = UVp(lmcur+l+m,ua)

       end do          

       lmcur = lmcur + 2*l + 1
    end do

  end subroutine surface_fields_full


  subroutine surface_field(UVp,Us)

    implicit none

    complex(dpc), dimension(-lmax:lmax,nglob_ssg,lmax), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2), intent(out) :: Us

    integer(i4b) :: l,ispec1,ua,lmcur,m
    real(dp) :: rl

    Us = 0.0_dp

    lmcur = 2
    do l = 1,lmax
       call get_depth(l,rl,ispec1)
       ua = gvn_ssg(ngll,nspec,ispec1,2)

       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle
    
          Us(lmcur+l+m) = UVp(m,ua,l)

       end do          

       lmcur = lmcur + 2*l + 1
    end do


  end subroutine surface_field


  subroutine rob_projection

    implicit none

    integer(i4b) :: i

    rob_coefs(0,1) = 1.0000_dp
    rob_coefs(0,2) = 0.0000_dp
    rob_coefs(1,1) = 0.9986_dp
    rob_coefs(1,2) = 0.0620_dp
    rob_coefs(2,1) = 0.9954_dp
    rob_coefs(2,2) = 0.1240_dp
    rob_coefs(3,1) = 0.9900_dp
    rob_coefs(3,2) = 0.1860_dp
    rob_coefs(4,1) = 0.9822_dp
    rob_coefs(4,2) = 0.2480_dp
    rob_coefs(5,1) = 0.9730_dp
    rob_coefs(5,2) = 0.3100_dp
    rob_coefs(6,1) = 0.9600_dp
    rob_coefs(6,2) = 0.3720_dp
    rob_coefs(7,1) = 0.9427_dp
    rob_coefs(7,2) = 0.4340_dp
    rob_coefs(8,1) = 0.9216_dp
    rob_coefs(8,2) = 0.4958_dp
    rob_coefs(9,1) = 0.8962_dp
    rob_coefs(9,2) = 0.5571_dp
    rob_coefs(10,1) = 0.8679_dp
    rob_coefs(10,2) = 0.6176_dp
    rob_coefs(11,1) = 0.8350_dp
    rob_coefs(11,2) = 0.6769_dp
    rob_coefs(12,1) = 0.7986_dp
    rob_coefs(12,2) = 0.7346_dp
    rob_coefs(13,1) = 0.7597_dp
    rob_coefs(13,2) = 0.7903_dp
    rob_coefs(14,1) = 0.7186_dp
    rob_coefs(14,2) = 0.8435_dp
    rob_coefs(15,1) = 0.6732_dp
    rob_coefs(15,2) = 0.8936_dp
    rob_coefs(16,1) = 0.6213_dp
    rob_coefs(16,2) = 0.9394_dp
    rob_coefs(17,1) = 0.5722_dp
    rob_coefs(17,2) = 0.9761_dp
    rob_coefs(18,1) = 0.5322_dp
    rob_coefs(18,2) = 1.0000_dp

    rob_lats = 0.0_dp

    do i = -18,-1
       rob_coefs(i,1) = rob_coefs(-i,1)
       rob_coefs(i,2) = -rob_coefs(-i,2)
       rob_lats(i) = 5.0_dp*i
       rob_lats(-i) = -5.0_dp*i
    end do

    call spline(rob_lats,rob_coefs(:,1),big,big,crob_coefs(:,1))
    call spline(rob_lats,rob_coefs(:,2),big,big,crob_coefs(:,2))



  end subroutine rob_projection


  subroutine get_xy(lat,long,xx,yy)

    implicit none

    real(dp), intent(in) :: lat,long
    real(dp), intent(out) :: xx,yy

    real(dp) :: plen,pdfe

    plen = splint(rob_lats,rob_coefs(:,1),crob_coefs(:,1),lat)
    pdfe = splint(rob_lats,rob_coefs(:,2),crob_coefs(:,2),lat)

    yy = pdfe*0.5072
    xx = long*plen/360.0_dp


  end subroutine get_xy


end module module_gia
