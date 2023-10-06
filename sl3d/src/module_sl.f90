module module_sl

  use nrtype
  use module_util
  use module_mat
  use module_model
  use module_function
  use module_fourier
  use module_gia
  
  implicit none

  real(dp), protected :: min_change = 0.01_dp
  complex(dpc), dimension(:), allocatable, save :: hall_load

  type meas

     integer(i4b) :: type
     ! 1 = SL, 2 = vertical displacement, 3 = gravity pert

     complex(dpc), dimension(:), allocatable :: load_lm

  end type meas


  type point_meas

     ! ALL SHOULD BE NORMALISED

     ! 1 = SL, 2 = vertical displacement, 3 = gravity pert
     integer(i4b) :: type

     ! number of times of measurements at location
     integer(i4b) :: nt

     ! location of measurements
     real(dp) :: lat,long

     ! time before present - should be in order from oldest
     real(dp), dimension(:), allocatable :: times

     complex(dpc), dimension(:), allocatable :: values

  end type point_meas


  type rsh_meas

     ! 1 = SL, 2 = vertical displacement, 3 = gravity pert
     integer(i4b) :: type

     ! number of times of measurements
     integer(i4b) :: nt

     ! time before present - should be in order from oldest
     real(dp), dimension(:), allocatable :: times

     ! size ((lmax+1)**2,nt)
     real(dp), dimension(:,:), allocatable :: values     

  end type rsh_meas


  type comb_meas

     logical(lgt) :: pres_point,pres_rsh
     integer(i4b) :: n_point,n_rsh
     type(point_meas), dimension(:), allocatable :: point
     type(rsh_meas), dimension(:), allocatable :: rsh

  end type comb_meas


contains


  subroutine set_ice_model(ice)

    implicit none

    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(out) :: ice

    integer(i4b) :: it,ilat,ilong,ilonged,io,igl,iphi
    real(dp) :: tmp,long,lat
    character(len=256) :: file1,icefile
    real(dp), dimension(nlong_ice5g,nlat_ice5g,nt_ice5g) :: ice_t

    do it = 1,nt_ice5g
       ! every 1000 years from 21000 to 17000, then every 500
     
       if (it < 6) then
          call string_cat_int('/raid/oc251/data/ice5g/ascii/ice5g_v1.2_',22-it,file1)
       else if (it < 20) then
          call string_cat_int('/raid/oc251/data/ice5g/ascii/ice5g_v1.2_',(nt_ice5g-it)/2,file1)
       else
          call string_cat_int('/raid/oc251/data/ice5g/ascii/ice5g_v1.2_0',(nt_ice5g-it)/2,file1)
       end if

       if ((mod(it,2) == 0) .and. (it > 5)) then
          icefile = trim(file1) // trim('.5k_1deg.sftgit.ascii')
       else
          icefile = trim(file1) // trim('.0k_1deg.sftgit.ascii')
       end if

       open(newunit=io,file=icefile,action='read',form='formatted')

       do ilat = 1,nlat_ice5g
          do ilong = 1,nlong_ice5g
             ilonged = ilong + (nlong_ice5g/2)
             if (ilonged > nlong_ice5g) ilonged = ilonged - nlong_ice5g
             !ilonged = ilong
             read(io,*) tmp
             ice_t(ilonged,ilat,it) = tmp

          end do
       end do
       close(io)

    end do

    ! normalise
    ice_t = ice_t/r_norm  

    ice = 0.0_dp

    do it = 1,nt_ice5g

       ! interpolate function
       do iphi = 1,nphi
          long = aphi(iphi)*rad2deg - 180.0_dp
          do igl = 1,ngl
             lat = (pio2_d - tgl(igl))*rad2deg
             ice(igl,iphi,it) = interp(ice_t(:,:,it),lat,long)
          end do
       end do

    end do
    

  end subroutine set_ice_model



  subroutine set_ice_model_6g(ice)

    implicit none

    complex(dpc), dimension(ngl,nphi,nt_ice6g), intent(out) :: ice

    integer(i4b) :: it,ilat,ilong,io,igl,iphi
    real(dp) :: tmp,long,lat
    character(len=256) :: file1,icefile
    real(dp), dimension(nlong_ice5g,nlat_ice5g,nt_ice6g) :: ice_t

    do it = 1,nt_ice6g
       ! every 1000 years from 26000 to 21000, then every 500
     
       if (it < 6) then
          call string_cat_int('/raid/oc251/data/ice5g/ascii/ICE6G.',27-it,file1)
       else
          call string_cat_int('/raid/oc251/data/ice5g/ascii/ICE6G.',(nt_ice5g-it)/2,file1)
       end if

       if (mod(it,2) == 0) then
          icefile = trim(file1) // trim('.0')
       else
          icefile = trim(file1) // trim('.5')
       end if

       open(newunit=io,file=icefile,action='read',form='formatted')

       do ilat = 1,nlat_ice5g
          do ilong = 1,nlong_ice5g
             read(io,*) tmp
             ice_t(ilong,ilat,it) = tmp

          end do
       end do
       close(io)

    end do

    ! normalise
    ice_t = ice_t/r_norm  

    ice = 0.0_dp

    do it = 1,nt_ice5g

       ! interpolate function
       do iphi = 1,nphi
          long = aphi(iphi)*rad2deg - 180.0_dp
          do igl = 1,ngl
             lat = (pio2_d - tgl(igl))*rad2deg
             ice(igl,iphi,it) = interp(ice_t(:,:,it),lat,long)
          end do
       end do

    end do
    

  end subroutine set_ice_model_6g



  subroutine memdot_bvec_sl(memm1,UVpm1,Wm1,dUm1,dVm1,dWm1,memdot,bUVp,bW,b_lm,ocean_lm)
    
    use nrtype
    implicit none


    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(in) :: dUm1,dVm1,dWm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: bUVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: bW
    complex(dpc), dimension((lmax+1)**2), intent(in), optional :: b_lm,ocean_lm

    integer(i4b) :: inode,ispec,ilayer,pa,ua,va,wa,l,m,ispec1,lmtot,lmcur,i,jnode,ub,vb,wb,lm,lmm
    real(dp) :: zeta2,zeta,zeta2m2,rl,rr,rmu,gsurf,rsurf
    complex(dpc) :: uu,vv,ww,dU,dV,dW,bint,area
    complex(dpc) :: tmp
    complex(dpc), dimension(5,(lmax+1)**2) :: dten
    complex(dpc), dimension(5,ngl,nphi) :: Aspat
    complex(dpc), dimension(5,(lmax+1)**2) :: Alm

    ! A(1) = A(-,-)
    ! A(2) = A(-,0)
    ! A(3) = A(-,+) = 0.5*A(0,0)
    ! A(4) = A(0,+)
    ! A(5) = A(+,+)

    rsurf = r_node(ngll,nspec)
    gsurf = grav_node(ngll,nspec)

    memdot = 0.0_dp
    bUVp = 0.0_dp
    bW = 0.0_dp

    ! Calculate spatial d - m
    do ilayer = 1,nsect

       ! No terms in elastic or fluid regions
       if (sect_ind(3,ilayer) /=2) cycle

       Alm = 0.0_dp

       !$OMP PARALLEL DO PRIVATE(pa,ua,va,wa,ispec1,lmtot,lmcur,i,zeta2,zeta,zeta2m2), &
       !$OMP& PRIVATE(rl,rr,rmu,uu,vv,ww,dU,dV,dW,ub,vb,wb,inode,ispec,l,m,jnode,lm,lmm,tmp), &
       !$OMP& FIRSTPRIVATE(dten,Aspat,Alm), &
       !$OMP& SHARED(lmax,memm1,UVpm1,Wm1,dUm1,dVm1,dWm1,memdot,bUVp,bW,ilayer), &
       !$OMP& SHARED(rsurf,gsurf)
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

                   !if ((real(ww) == 0.0_dp) .and. (imag(ww) == 0.0_dp)) ww = 0.0_dp
                   !if ((real(dW) == 0.0_dp) .and. (imag(dW) == 0.0_dp)) dW = 0.0_dp

                   !ww = 0.0_dp
                   !dW = 0.0_dp

                   !if ((ispec == nspec/2) .and. (inode == 1) .and. (l==2)) then
                   !   print *, m,uu,vv,dU,dV,ww,dW
                   !end if

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

             !if (.not. struct3d) then
             !   do lm = 1,(lmax+1)**2
             !      tmp = 0.5_dp*(dten(1,lm) + dten(5,lm))

             !      dten(1,lm) = tmp
             !      dten(5,lm) = tmp
                   
             !      tmp = 0.5_dp*(dten(2,lm) + dten(4,lm))
             !      dten(2,lm) = tmp
             !      dten(4,lm) = tmp
             !   end do

             !end if

             call fun_from_coefs_ten_tl(dten,Aspat)


             do i = 1,5
                Aspat(i,:,:) = si_spat_node(:,:,1,inode,ispec) &
                               *(Aspat(i,:,:) - memm1(i,:,:,1,inode,ispec))


                ! Calculate memdot
                memdot(i,:,:,1,inode,ispec) = Aspat(i,:,:)

             end do
 
             ! Calculate Alm
             call coefs_from_fun_ten_tl(Aspat,Alm(:,:))

             do l = 1,lmax
                do m = -l,l
                   lm = l**2 + l + m + 1
                   lmm = l**2 + l - m + 1

                   tmp = 0.5_dp*(Alm(3,lm) + ((-1)**m)*conjg(Alm(3,lmm)))
                   Alm(3,lm) = tmp
                   Alm(3,lmm) = conjg(tmp)*((-1)**m)

                   tmp = 0.5_dp*(Alm(1,lm) + ((-1)**m)*conjg(Alm(5,lmm)))
                   Alm(1,lm) = tmp
                   Alm(5,lmm) = conjg(tmp)*((-1)**m)

                   tmp = 0.5_dp*(Alm(2,lm) + ((-1)**m)*conjg(Alm(4,lmm)))
                   Alm(2,lm) = tmp
                   Alm(4,lmm) = conjg(tmp)*((-1)**m)
                
                end do
             end do

             !if ((ispec == nspec/2) .and. (inode == 1)) then
             !   do m = -2,2
             !      print *, m,rr*(Alm(1,7+m) + Alm(5,7+m))/(2.0_dp*sqrt(6.0_dp)), &
             !               ii*rr*(Alm(1,7+m) - Alm(5,7+m))/(2.0_dp*sqrt(6.0_dp)), &
             !               3.0_dp*rr*Alm(3,7+m), &
             !               rr*sqrt(2.0_dp)*(Alm(2,7+m) + Alm(4,7+m))/(sqrt(6.0_dp)), &
             !               rr*ii*sqrt(2.0_dp)*(Alm(2,7+m) - Alm(4,7+m))/(sqrt(6.0_dp))
             !   end do
             !end if


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

    ! Add load terms if present
    if (present(b_lm)) then
       bint = b_lm(1)*sqrt(fourpi_d)*(rsurf**2)
       area = ocean_lm(1)*sqrt(fourpi_d)*(rsurf**2)
       lmtot = 2
       do l = 1,lmax
          pa = pa_surf(l)
          ua = pa + 1          
          do m = -l,l
             lmcur = lmtot + l + m
             bUVp(lmcur,pa) = bUVp(lmcur,pa) - rice*b_lm(lmcur)
             bUVp(lmcur,pa) = bUVp(lmcur,pa) + rice*bint*ocean_lm(lmcur)/area
             bUVp(lmcur,ua) = bUVp(lmcur,ua) - rice*b_lm(lmcur)*gsurf
             bUVp(lmcur,ua) = bUVp(lmcur,ua) + rice*bint*ocean_lm(lmcur)*gsurf/area

          end do
          lmtot = lmtot + 2*l + 1
       end do
    end if

  end subroutine memdot_bvec_sl



  subroutine time_derivs_full_sl(ice_load,ice_loaddot,sea_level,memm1,UVpm1,Wm1, &
                                sldot,memdot,UVpdot,Wdot)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,ice_loaddot,sea_level
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVpdot
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: Wdot
    complex(dpc), dimension(ngl,nphi), intent(out) :: sldot

    real(dp) :: area
    complex(dpc), dimension(ngl,nphi) :: ocean_func,B_func,U_func,p_func,CUp
    complex(dpc), dimension((lmax+1)**2) :: ocean_func_lm,B_func_lm,CUp_lm,U_lm,p_lm
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: bUVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: bW

    integer(i4b) :: l,ispec1,lmcur,lmact,kd,kl,ku,ldab,nrows,ua,va,pa,info,m,ispec11,ilayer, &
                    ispec,inode,i,iter,io,igl,iphi
    real(dp) :: rl,rr,gsurf,ratio,rsurf
    real(dp), dimension(ntau) :: rsi
    complex(dpc) :: int,bint
    complex(dpc), dimension(:,:), allocatable :: A1u,A1u_trans
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpcur,UVpnew,UVp0

    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV,dW

    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: UVp_trans
    complex(dpc), dimension(nglob_tor,-lmax:lmax) :: W_trans

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    ! Calculate C
    ocean_func = ocean_function(ice_load,sea_level)

    ! Calculate spherical harmonic coefficients of C
    call coefs_from_fun(ocean_func,ocean_func_lm)

    ! Calculate A
    area = sqrt(fourpi_d)*ocean_func_lm(1)*(rsurf**2)

    ! Calculate B = (1-C)*Idot
    B_func = B_function(ice_loaddot,ocean_func)

    ! Calculate spherical harmonic coefficients of B
    call coefs_from_fun(B_func,B_func_lm)

    bint = B_func_lm(1)*sqrt(fourpi_d)*(rsurf**2)

    ! Calculate stresses
    call calculate_dUVW_full(UVpm1,Wm1,dU,dV,dW)

    !-------------------------------------!
    !      Calculate memdot and RHSs      !
    !-------------------------------------!

    call memdot_bvec_sl(memm1,UVpm1,Wm1,dU,dV,dW,memdot,bUVp,bW, &
                         B_func_lm,ocean_func_lm)

    print *, "Calculating u0"

    do l = 1,lmax

       call get_depth(l,rl,ispec1)
       
       !---------------------------------------!
       !      Calculate spheroidal modes       !
       !---------------------------------------!

       kl = 3*(ngll-1) + 2
       ku = 3*(ngll-1) + 2
       kd = ku
    
       nrows = gvn_ssg(ngll,nspec,ispec1,3)

       UVp_trans = 0.0_dp
       UVp_trans(1:nrows,-l:l) = transpose(bUVp(l**2+1:(l+1)**2,1:nrows))

       pa = pa_surf(l)
       ua = pa+1
       if (l == 2) then
          call get_depth(2,rl,ispec1)
          ua = gvn_ssg(ngll,nspec/2,ispec1,2)
          print *, UVp_trans(ua,0)

          open(newunit=io,file='bvectest10in.out')
          call get_depth(2,rl,ispec1)
          do ilayer = 3,nsect
             do ispec = max(sect_ind(1,ilayer),ispec1), sect_ind(2,ilayer)
                do inode = 1,ngll
                   ua = gvn_ssg(inode,ispec,ispec1,2)
                   write(io,*) r_node(inode,ispec),real(UVp_trans(ua,0)),imag(UVp_trans(ua,0))
                end do
             end do
          end do
          close(io)
       end if

       UVp_trans = UVp_trans*(10.0_dp**15)

       if (l == 1) then
          ldab = 2*kl + ku + 1
          call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                      ipiv_ssg_1(1:nrows),UVp_trans(1:nrows,-l:l),nrows,info)
          

       else
          ldab = kd + 1             
          call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                      UVp_trans(1:nrows,-l:l),nrows,info)
       end if

       UVp_trans = UVp_trans/(10.0_dp**15)

       if (l == 2) then
          call get_depth(2,rl,ispec1)
          ua = gvn_ssg(ngll,nspec/2,ispec1,2)
          print *, UVp_trans(ua,0)

          open(newunit=io,file='bvectest10out.out')
          call get_depth(2,rl,ispec1)
          do ilayer = 3,nsect
             do ispec = max(sect_ind(1,ilayer),ispec1), sect_ind(2,ilayer)
                do inode = 1,ngll
                   ua = gvn_ssg(inode,ispec,ispec1,2)
                   write(io,*) r_node(inode,ispec),real(UVp_trans(ua,0)),imag(UVp_trans(ua,0))
                end do
             end do
          end do
          close(io)
       end if

       UVpcur(l**2+1:(l+1)**2,1:nrows) = transpose(UVp_trans(1:nrows,-l:l))

       !---------------------------------------!
       !       Calculate toroidal modes        !
       !---------------------------------------!

       if (struct3d) then
     
          kl = ngll-1
          ku = ngll-1
          kd = ku
    
          nrows = gvn_tor(ngll,nspec,ispec1)
          
          W_trans(1:nrows,-l:l) = transpose(bW(l**2+1:(l+1)**2,1:nrows))

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
       else
          
          Wdot = 0.0_dp

       end if

    end do


    !================================================!
    !       ITERATE TO FIND ACTUAL DERIVATIVES       !
    !================================================!
    
    UVp0 = UVpcur

    allocate(A1u(-lmax:lmax,nglob_ssg))
    allocate(A1u_trans(nglob_ssg,-lmax:lmax))

    A1u = 0.0_dp
    A1u_trans = 0.0_dp

    ratio = 1.0_dp
    iter = 0
    do while (ratio > min_change)
       iter = iter + 1
       
       !-------------------------------------------!
       !         Calculate A perturbation          !
       !-------------------------------------------!
       
       ! spatial U and phi
       U_lm = 0.0_dp
       p_lm = 0.0_dp   

       lmcur = 2
       do l = 1,lmax

          pa = pa_surf(l)
          ua = pa + 1
          
          do m = -lmax,lmax
             if (m > l) exit
             if (abs(m) > l) cycle
             lmact = lmcur + l + m

             p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
             U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

          end do
          lmcur = lmcur + 2*l + 1
       end do

       p_func = 0.0_dp
       U_func = 0.0_dp  
      
       call fun_from_coefs(p_lm,p_func)
       call fun_from_coefs(U_lm,U_func)  

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func

       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun(CUp,CUp_lm)

       ! Calculate int(C*Updot)
       int = sqrt(fourpi_d)*CUp_lm(1)*(rsurf**2)


       !--------------------------------!
       !          Find new UVp          !
       !--------------------------------!

       UVpnew = UVp0

       ! l = 1 to lmax
       lmcur = 2

       kl = 3*(ngll-1) + 2
       ku = kl
       kd = kl

       do l = 1,lmax

          pa = pa_surf(l)
          ua = pa + 1
          nrows = ua + 1

          A1u_trans = 0.0_dp

          do m = -lmax,lmax
             if (m > l) exit
             if (abs(m) > l) cycle
             lmact = lmcur + l + m
             A1u_trans(pa,m) = A1u_trans(pa,m) + ocean_func_lm(lmact)*int*roce/(area*gsurf)
             A1u_trans(pa,m) = A1u_trans(pa,m) - CUp_lm(lmact)*roce/gsurf
  
             A1u_trans(ua,m) = A1u_trans(ua,m) + ocean_func_lm(lmact)*int*roce/(area)
             A1u_trans(ua,m) = A1u_trans(ua,m) - CUp_lm(lmact)*roce

          end do

          ! Perform A0-1 A1 u

          if (l == 1) then
             ldab = 2*kl + ku + 1
             call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                         ipiv_ssg_1(1:nrows),A1u_trans(1:nrows,-l:l),nrows,info)
          else
             ldab = kd + 1             
             call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                         A1u_trans(1:nrows,-l:l),nrows,info)
          end if


          UVpnew(lmcur:lmcur+2*l,:) = UVpnew(lmcur:lmcur+2*l,:) &
                                        - transpose(A1u_trans(:,-l:l))

          lmcur = lmcur + 2*l + 1

       end do

       ua = pa_surf(lmax/2) + 1

       ratio = abs(UVpcur((lmax/2)**2+(lmax/2)+1,ua) - UVpnew((lmax/2)**2+(lmax/2)+1,ua)) &
               /abs(UVpnew((lmax/2)**2+(lmax/2)+1,ua) - UVp0((lmax/2)**2+(lmax/2)+1,ua))

       print *, "Iteration = ",iter, ", Ratio = ",ratio

       UVpcur = UVpnew

    end do

    UVpdot = UVpcur
    

    !============================================!
    !    CALCULATE RATE OF CHANGE OF SEA LEVEL   !
    !============================================!

    gsurf = grav_node(ngll,nspec)

    ! spatial U and phi
    U_lm = 0.0_dp
    p_lm = 0.0_dp

    lmcur = 2     

    do l = 1,lmax

       pa = pa_surf(l)
       ua = pa + 1
          
       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle
          lmact = lmcur + l + m

          p_lm(lmact) = UVpdot(lmact,pa)
          U_lm(lmact) = UVpdot(lmact,ua)

       end do
       lmcur = lmcur + 2*l + 1
    end do

    p_func = 0.0_dp
    U_func = 0.0_dp  
      
    call fun_from_coefs(p_lm,p_func)
    call fun_from_coefs(U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func

    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(CUp,CUp_lm(1))
    
    sldot = real(sqrt(fourpi_d)*(rsurf**2)*((CUp_lm(1)/gsurf) - (rice*B_func_lm(1)/roce))/area &
                 - U_func - p_func/gsurf)

    !print *, CUp_lm(1),B_func_lm(1),area,U_func(20,20),p_func(20,20),sldot(20,20)

  end subroutine time_derivs_full_sl


  subroutine update_full_rk_sl(dt,ice_load,ice_loaddot,slm1,memm1,UVpm1,Wm1,sl,mem,UVp,W)
  ! routine to calculate the displacements from the previous displacements
  ! using Euler time-stepping  

    use nrtype

    implicit none

    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,ice_loaddot
    complex(dpc), dimension(ngl,nphi), intent(in) :: slm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(ngl,nphi), intent(out) :: sl
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: mem
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: W

    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot,Wtmp
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot,memtmp
    complex(dpc), dimension(ngl,nphi) :: sldot,sltmp,ice_loadtmp

    integer(i4b) :: ua

    UVpdot = 0.0_dp
    Wdot = 0.0_dp
    memdot = 0.0_dp

    ! Calculate time derivatives
    call time_derivs_full_sl(ice_load,ice_loaddot,slm1,memm1,UVpm1,Wm1, &
                             sldot,memdot,UVpdot,Wdot)

    UVptmp = UVpm1 + UVpdot*dt/2.0_dp

    Wtmp = Wm1 + Wdot*dt/2.0_dp

    memtmp = memm1 + memdot*dt/2.0_dp

    sltmp = slm1 + sldot*dt/2.0_dp

    ice_loadtmp = ice_load + ice_loaddot*dt/2.0_dp

    ! Calculate time derivatives
    call time_derivs_full_sl(ice_loadtmp,ice_loaddot,sltmp,memtmp,UVptmp,Wtmp, &
                             sldot,memdot,UVpdot,Wdot)

    UVp = UVpm1 + UVpdot*dt

    W = Wm1 + Wdot*dt

    mem = memm1 + memdot*dt

    sl = slm1 + sldot*dt
  
  end subroutine update_full_rk_sl


  subroutine forward_full_sl(nt,dt,nout,ice_data,sl_0,sl_out,UVp_out,W_out,mem_out,nout_opt)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_0
    complex(dpc), dimension(ngl,nphi,nout), intent(out) :: sl_out
    complex(dpc), dimension(:,:,:), intent(out), optional :: UVp_out
    complex(dpc), dimension(:,:,:), intent(out), optional :: W_out
    complex(dpc), dimension(:,:,:,:,:,:,:), intent(out), optional :: mem_out
    integer(i4b), intent(in), optional :: nout_opt

    integer(i4b) :: it,it_out,iout,it_out_uwm,iout_uwm
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor,2) :: W
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,2) :: mem
    complex(dpc), dimension(ngl,nphi,2) :: sl_f
    complex(dpc), dimension(ngl,nphi) :: ice_loaddot,ice_load

    UVp = 0.0_dp
    W = 0.0_dp
    mem = 0.0_dp
    sl_f = 0.0_dp
    ice_loaddot = 0.0_dp
    ice_load = 0.0_dp

    sl_out = 0.0_dp

    if (present(UVp_out)) UVp_out = 0.0_dp
    if (present(W_out)) W_out = 0.0_dp
    if (present(mem_out)) mem_out = 0.0_dp

    ! Calculate it_out
    it_out = (nt-1)/(nout-1)

    if (present(nout_opt)) then
       it_out_uwm = (nt-1)/(nout_opt-1)
    else
       it_out_uwm = it_out
    end if

    sl_f(:,:,1) = sl_0
    sl_out(:,:,1) = sl_0

    iout = 1
    iout_uwm = 1
    do it = 2,nt
       print *, "Forward time step: ",it," of ",nt 
       call time_point_interp(it-1,dt,ice_data,ice_load,ice_loaddot)
       call update_full_rk_sl(dt,ice_load,ice_loaddot,sl_f(:,:,1),mem(:,:,:,:,:,:,1), &
                              UVp(:,:,1),W(:,:,1),sl_f(:,:,2),mem(:,:,:,:,:,:,2), &
                              UVp(:,:,2),W(:,:,2))
       
       if (mod(it-1,it_out) == 0) then
          iout = iout + 1
          sl_out(:,:,iout) = sl_f(:,:,2)

       end if

       if (mod(it-1,it_out_uwm) == 0) then
          iout_uwm = iout_uwm + 1

          if (present(UVp_out)) UVp_out(:,:,iout_uwm) = UVp(:,:,2)
          if (present(W_out)) W_out(:,:,iout_uwm) = W(:,:,2)
          if (present(mem_out)) mem_out(:,:,:,:,:,:,iout_uwm) = mem(:,:,:,:,:,:,2)

       end if

       mem(:,:,:,:,:,:,1) = mem(:,:,:,:,:,:,2)
       UVp(:,:,1) = UVp(:,:,2)
       W(:,:,1) = W(:,:,2)
       sl_f(:,:,1) = sl_f(:,:,2)

    end do


  end subroutine forward_full_sl


  !=======================================================!
  !                   ADJOINT ROUTINES                    !
  !=======================================================!


  subroutine point_load(lat,long,damp,load_lm)

    integer(i4b), intent(in) :: damp
    real(dp), intent(in) :: lat,long
    complex(dpc), dimension((lmax+1)**2), intent(out) :: load_lm

    integer(i4b) :: l,m,lmcur
    real(dp) :: theta,phi,c1,c2

    real(dp), dimension(1,2*lmax+1) :: dnm

    theta = (90.0_dp - lat)*deg2rad
    phi = (long + 180.0_dp)*deg2rad

    lmcur = 1

    load_lm = 0.0_dp

    do l = 0,lmax

       dnm = 0.0_dp
       call rotmx2(0,l,theta,dnm(1,1:2*l+1),1,2*l+1)

       c1 = sqrt((2*l+1)/(fourpi_d))

       if (damp == 0) then
          c2 = 1.0_dp
       else
          c2 = exp(-twopi_d*(l+1)/(damp + 0.5_dp))
       end if

       do m = -l,l
          load_lm(lmcur+l+m) = c1*dnm(1,l+m+1)*exp(-ii*m*phi)*c2
       end do

       lmcur = lmcur + 2*l + 1

    end do

  end subroutine point_load


  subroutine stat_adj(ice_load,sl_for,load,UVp_adj,Q_adj,K_adj)

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_for
    type(meas), intent(in) :: load
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp_adj
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: Q_adj
    complex(dpc), intent(out), optional :: K_adj

    complex(dpc), dimension(ngl,nphi) :: ocean_func,p_func,U_func,CUp,hall_spat
    complex(dpc), dimension((lmax+1)**2) :: ocean_func_lm,CUp_lm,p_lm,U_lm
    integer(i4b) :: l,m,igl,iphi,io,lmcur,pa,ua,kl,ku,kd,ldab,nrows,ispec1,info,iter
    real(dp) :: c1,theta,phi,gsurf,area,rl,ilat,ilong,rsurf,ratio,c2
    complex(dpc) :: int
    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: bvec
    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: A1u_trans
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpcur,UVpnew,UVp0

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    ! Calculate C
    ocean_func = ocean_function(ice_load,sl_for)

    ! Calculate spherical harmonic coefficients of C
    call coefs_from_fun(ocean_func,ocean_func_lm)

    ! Calculate A
    area = sqrt(fourpi_d)*ocean_func_lm(1)*(rsurf**2)

    UVp_adj = 0.0_dp
    UVpcur = 0.0_dp
    UVpnew = 0.0_dp
    UVp0 = 0.0_dp

    lmcur = 2

    kd = 3*(ngll-1) + 2
    kl = kd
    ku = kd

    do l = 1,lmax
       bvec = 0.0_dp

       pa = pa_surf(l)
       ua = pa + 1
       nrows = ua + 1


       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle

          if (load%type == 1) then

             bvec(pa,m) = load%load_lm(lmcur+l+m)/gsurf &
                          - sqrt(fourpi_d)*load%load_lm(1)*ocean_func_lm(lmcur+l+m)/(area*gsurf)
             bvec(ua,m) = load%load_lm(lmcur+l+m) &
                          - sqrt(fourpi_d)*load%load_lm(1)*ocean_func_lm(lmcur+l+m)/area

          else if (load%type == 2) then

             bvec(ua,m) = -load%load_lm(lmcur+l+m)

          else if (load%type == 3) then

             bvec(pa,m) = -load%load_lm(lmcur+l+m)

          end if

       end do

       if (l == 1) then
          ldab = 2*kl + ku + 1
          call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                      ipiv_ssg_1(1:nrows),bvec(1:nrows,-l:l),nrows,info)
       else
          ldab = kd + 1             
          call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                      bvec(1:nrows,-l:l),nrows,info)
       end if
          
       UVpcur(l**2+1:(l+1)**2,1:nrows) = transpose(bvec(1:nrows,-l:l))
      
       lmcur = lmcur + 2*l + 1

    end do


    !================================================!
    !       ITERATE TO FIND ACTUAL DISPLACEMENT      !
    !================================================!
    
    UVp0 = UVpcur

    ratio = 1.0_dp
    iter = 0
    do while (ratio > min_change)
       iter = iter + 1
       
       !-------------------------------------------!
       !         Calculate A perturbation          !
       !-------------------------------------------!
       
       ! spatial U and phi
       U_lm = 0.0_dp
       p_lm = 0.0_dp  

       lmcur = 2

       do l = 1,lmax
          pa = pa_surf(l)
          ua = pa + 1
          
          do m = -lmax,lmax
             if (m > l) exit
             if (abs(m) > l) cycle

             p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
             U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

          end do
          lmcur = lmcur + 2*l + 1
       end do

       p_func = 0.0_dp
       U_func = 0.0_dp  
      
       call fun_from_coefs(p_lm,p_func)
       call fun_from_coefs(U_lm,U_func)  

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func

       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun(CUp,CUp_lm)

       ! Calculate int(C*Updot)
       int = sqrt(fourpi_d)*CUp_lm(1)*(rsurf**2)


       !--------------------------------!
       !          Find new UVp          !
       !--------------------------------!

       UVpnew = UVp0

       kl = 3*(ngll-1) + 2
       ku = kl
       kd = kl

       lmcur = 2
       do l = 1,lmax

          pa = pa_surf(l)
          ua = pa + 1
          nrows = ua + 1

          A1u_trans = 0.0_dp

          do m = -lmax,lmax
             if (m > l) exit
             if (abs(m) > l) cycle  
             A1u_trans(pa,m) = A1u_trans(pa,m) + ocean_func_lm(lmcur+l+m)*int*roce/(area*gsurf)
             A1u_trans(pa,m) = A1u_trans(pa,m) - CUp_lm(lmcur+l+m)*roce/gsurf

             A1u_trans(ua,m) = A1u_trans(ua,m) + ocean_func_lm(lmcur+l+m)*int*roce/(area)
             A1u_trans(ua,m) = A1u_trans(ua,m) - CUp_lm(lmcur+l+m)*roce

          end do

          ! Perform A0-1 A1 u
          if (l == 1) then
             ldab = 2*kl + ku + 1
             call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                         ipiv_ssg_1(1:nrows),A1u_trans(1:nrows,-l:l),nrows,info)
          else
             ldab = kd + 1             
             call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                         A1u_trans(1:nrows,-l:l),nrows,info)
          end if

          UVpnew(l**2+1:(l+1)**2,1:nrows) = UVpnew(l**2+1:(l+1)**2,1:nrows) &
                                            - transpose(A1u_trans(1:nrows,-l:l))

          lmcur = lmcur + 2*l + 1

       end do

       ua = pa_surf(lmax) + 1

       ratio = abs(UVpcur(lmax**2+lmax+1,ua) - UVpnew(lmax**2+lmax+1,ua)) &
               /abs(UVpnew(lmax**2+lmax+1,ua) - UVp0(lmax**2+lmax+1,ua))

       UVpcur = UVpnew

       print *, "Iteration = ",iter,", Ratio = ",ratio 

    end do

    UVp_adj = UVpcur


    if (.not. present(Q_adj)) return

    print *, "Here"

    ! spatial U and phi
    U_lm = 0.0_dp
    p_lm = 0.0_dp  
    
    lmcur = 2

    do l = 1,lmax
       pa = pa_surf(l)
       ua = pa + 1
          
       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle

          p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
          U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

       end do
       lmcur = lmcur + 2*l + 1
    end do

    p_func = 0.0_dp
    U_func = 0.0_dp  
      
    call fun_from_coefs(p_lm,p_func)
    call fun_from_coefs(U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func
    
    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(CUp,CUp_lm(1))

    if (load%type == 1) then

       call fun_from_coefs(load%load_lm,hall_spat)
    
       ! Calculate K and Q
       K_adj = sqrt(fourpi_d)*CUp_lm(1)*(rsurf**2)/(gsurf*area) &
               + sqrt(fourpi_d)*load%load_lm(1)*(rsurf**2)/(gsurf*area*roce)

       Q_adj = CUp/gsurf &
               - ocean_func*sqrt(fourpi_d)**(rsurf**2)*(load%load_lm(1)/roce + CUp_lm(1))/(gsurf*area) &
               + hall_spat/(roce*gsurf)

    else if (load%type == 2 .or. load%type == 3) then

       ! Calculate K and Q
       K_adj = sqrt(fourpi_d)*(rsurf**2)*CUp_lm(1)/(gsurf*area)

       Q_adj = CUp/gsurf - ocean_func*sqrt(fourpi_d)*(rsurf**2)*(CUp_lm(1))/(gsurf*area)

    end if

  end subroutine stat_adj


  subroutine time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,Wm1,memdot,UVpdot,Wdot,Qdot,Kdot)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,sl_for
    complex(dpc), dimension(5,ngl,nphi,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVpdot
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: Wdot
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: Qdot
    complex(dpc), intent(out), optional :: Kdot

    real(dp) :: area
    complex(dpc), dimension(ngl,nphi) :: U_func,p_func,CUp,ocean_func
    complex(dpc), dimension((lmax+1)**2) :: ocean_func_lm,CUp_lm,U_lm,p_lm
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: bUVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: bW

    integer(i4b) :: l,ispec1,lmcur,kl,ku,kd,ldab,nrows,ua,va,pa,info,m,ispec11,ilayer, &
                    ispec,inode,i,iter,io,igl,iphi
    real(dp) :: zeta2,zeta,c,rl,rr,gsurf,long,lat,ratio,rsurf
    real(dp), dimension(ntau) :: rsi
    complex(dpc) :: int
    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: A1u_trans
    complex(dpc), dimension(-lmax:lmax,nglob_ssg) :: A1u
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpcur,UVpnew,UVp0

    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV,dW

    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: UVp_trans
    complex(dpc), dimension(nglob_tor,-lmax:lmax) :: W_trans

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    ! Calculate C
    ocean_func = ocean_function(ice_load,sl_for)

    ! Calculate spherical harmonic coefficients of C
    call coefs_from_fun(ocean_func,ocean_func_lm)

    ! Calculate A
    area = sqrt(fourpi_d)*ocean_func_lm(1)*(rsurf**2)

    ! Calculate stresses
    call calculate_dUVW_full(UVpm1,Wm1,dU,dV,dW)

    !-------------------------------------!
    !      Calculate memdot and RHSs      !
    !-------------------------------------!

    bUVp = 0.0_dp
    bW = 0.0_dp

    memdot = 0.0_dp

    call memdot_bvec_sl(memm1,UVpm1,Wm1,dU,dV,dW,memdot,bUVp,bW)

    UVpcur = 0.0_dp       

    print *, "Calculating u0"

    do l = 1,lmax

       call get_depth(l,rl,ispec1)

       !---------------------------------------!
       !      Calculate spheroidal modes       !
       !---------------------------------------!

       kl = 3*(ngll-1) + 2
       ku = kl
       kd = kl
   
       nrows = gvn_ssg(ngll,nspec,ispec1,3)

       ! Solve
       UVp_trans = 0.0_dp
       UVp_trans(1:nrows,-l:l) = transpose(bUVp(l**2+1:(l+1)**2,1:nrows))

       !if (l == 1) print *, nrows,bUVp(l**2+1:(l+1)**2,nrows/2)

       !if (l==1) then
       !   print *, nrows, bUVp(2:4,nrows/2)
       !end if

       if (l == 1) then
          ldab = 2*kl + ku + 1
          call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                      ipiv_ssg_1(1:nrows),UVp_trans(1:nrows,-l:l),nrows,info)

       else
          ldab = kd + 1             
          call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                      UVp_trans(1:nrows,-l:l),nrows,info)
          
       end if

       UVpcur(l**2+1:(l+1)**2,1:nrows) = transpose(UVp_trans(1:nrows,-l:l))

       !---------------------------------------!
       !       Calculate toroidal modes        !
       !---------------------------------------!

       if (struct3d) then
     
          kl = ngll-1
          ku = ngll-1
          kd = ku
    
          nrows = gvn_tor(ngll,nspec,ispec1)
          
          W_trans = 0.0_dp
          W_trans(1:nrows,-l:l) = transpose(bW(l**2+1:(l+1)**2,1:nrows))

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

       else

          Wdot = 0.0_dp

       end if

    end do

    !================================================!
    !       ITERATE TO FIND ACTUAL DERIVATIVES       !
    !================================================!
    
    UVp0 = UVpcur

    A1u = 0.0_dp
    A1u_trans = 0.0_dp

    ratio = 1.0_dp
    iter = 0
    do while (ratio > min_change)
       iter = iter + 1
       
       !-------------------------------------------!
       !         Calculate A perturbation          !
       !-------------------------------------------!
       
       ! spatial U and phi
       U_lm = 0.0_dp
       p_lm = 0.0_dp  

       lmcur = 2

       do l = 1,lmax
          pa = pa_surf(l)
          ua = pa + 1
          
          do m = -l,l
             !if (m > l) exit
             !if (abs(m) > l) cycle

             p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
             U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

          end do
          lmcur = lmcur + 2*l + 1
       end do

       p_func = 0.0_dp
       U_func = 0.0_dp  
      
       call fun_from_coefs(p_lm,p_func)
       call fun_from_coefs(U_lm,U_func) 

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func


       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun(CUp,CUp_lm)

       !print *, real(U_lm(100:102)),real(p_lm(100:102)),real(ocean_func_lm(100:102)),&
       !         real(CUp_lm(100:102))

       ! Calculate int(C*Updot)
       int = sqrt(fourpi_d)*CUp_lm(1)*(rsurf**2)


       !--------------------------------!
       !          Find new UVp          !
       !--------------------------------!

       UVpnew = UVp0

       ! l = 1 to lmax
       lmcur = 2

       kl = 3*(ngll-1) + 2
       ku = kl
       kd = kl

       do l = 1,lmax

          pa = pa_surf(l)
          ua = pa + 1
          nrows = ua + 1

          A1u_trans = 0.0_dp

          do m = -l,l
             !if (m > l) exit
             !if (abs(m) > l) cycle  
             A1u_trans(pa,m) = A1u_trans(pa,m) + ocean_func_lm(lmcur+l+m)*int*roce/(area*gsurf)
             A1u_trans(pa,m) = A1u_trans(pa,m) - CUp_lm(lmcur+l+m)*roce/gsurf

             A1u_trans(ua,m) = A1u_trans(ua,m) + ocean_func_lm(lmcur+l+m)*int*roce/(area)
             A1u_trans(ua,m) = A1u_trans(ua,m) - CUp_lm(lmcur+l+m)*roce

          end do

          ! Perform A0-1 A1 u
          if (l == 1) then
             ldab = 2*kl + ku + 1
             call zgbtrs('N',nrows,kl,ku,2*l+1,aa_ssg_1(:,1:nrows),ldab, &
                         ipiv_ssg_1(1:nrows),A1u_trans(1:nrows,-l:l),nrows,info)
          else
             ldab = kd + 1             
             call zpbtrs('U',nrows,kd,2*l+1,aa_ssg(:,1:nrows,l),ldab, &
                         A1u_trans(1:nrows,-l:l),nrows,info)
          end if

          UVpnew(lmcur:lmcur+2*l,:) = UVpnew(lmcur:lmcur+2*l,:) - transpose(A1u_trans(:,-l:l))

          lmcur = lmcur + 2*l + 1

       end do

       ua = pa_surf(lmax/2) + 1

       ratio = abs(UVpcur((lmax/2)**2+(lmax/2)+1,ua) - UVpnew((lmax/2)**2+(lmax/2)+1,ua)) &
               /abs(UVpnew((lmax/2)**2+(lmax/2)+1,ua) - UVp0((lmax/2)**2+(lmax/2)+1,ua))

       !ua = pa_surf(2) + 1

       !ratio = abs(UVpcur(0,ua,2) - UVpnew(0,ua,2)) &
       !        /abs(UVpnew(0,ua,2) - UVp0(0,ua,2))

       print *, "Iteration = ",iter,", Ratio = ",ratio 


       UVpcur = UVpnew


    end do

    UVpdot = UVpcur

    if(.not. present(Qdot)) return

    Qdot = 0.0_dp
    Kdot = 0.0_dp

    ! spatial U and phi
    U_lm = 0.0_dp
    p_lm = 0.0_dp   

    lmcur = 2

    do l = 1,lmax
       call get_depth(l,rl,ispec1)
       pa = gvn_ssg(ngll,nspec,ispec1,1)
       ua = pa + 1
          
       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle

          p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
          U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

       end do
       lmcur = lmcur + 2*l + 1
    end do

    p_func = 0.0_dp
    U_func = 0.0_dp  
      
    call fun_from_coefs(p_lm,p_func)
    call fun_from_coefs(U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func

    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(CUp,CUp_lm(1))

    ! AGAIN ASSUMING NO TIME DEPENDENCE OF ADJOINT LOAD

    Kdot = sqrt(fourpi_d)*(rsurf**2)*CUp_lm(1)/(gsurf*area)

    Qdot = CUp/gsurf - ocean_func*sqrt(fourpi_d)*(rsurf**2)*CUp_lm(1)/(gsurf*area)
                                      

  end subroutine time_derivs_full_adj


  subroutine update_full_rk_adj(dt,ice_load,sl_for,memm1,UVpm1,Wm1,mem,UVp,W,Qm1,Km1,Qt,Kt)
  ! routine to calculate the displacements from the previous displacements
  ! using Euler time-stepping  

    use nrtype

    implicit none

    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_for
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: Wm1
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(out) :: mem
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(out) :: W
    complex(dpc), dimension(ngl,nphi), intent(in), optional :: Qm1
    complex(dpc), intent(in), optional :: Km1
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: Qt
    complex(dpc), intent(out),optional :: Kt

    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot,Wtmp
    complex(dpc), dimension(ngl,nphi) :: Qdot
    complex(dpc) :: Kdot
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot,memtmp
    complex(dpc), dimension(5,(lmax+1)**2) :: memcoefs

    integer(i4b) :: io,igl,iphi


    if (present(Qm1)) then

       call time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,Wm1,memdot,UVpdot,Wdot)

       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       Wtmp = Wm1 + Wdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp

       call time_derivs_full_adj(ice_load,sl_for,memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot, &
                                 Qdot,Kdot) 

       UVp = UVpm1 + UVpdot*dt
       W = Wm1 + Wdot*dt
       mem = memm1 + memdot*dt
       Qt = Qm1 + Qdot*dt
       Kt = Km1 + Kdot*dt

    else

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,Wm1,memdot,UVpdot,Wdot)
       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       Wtmp = Wm1 + Wdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp

       call time_derivs_full_adj(ice_load,sl_for,memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot)
       UVp = UVpm1 + UVpdot*dt
       W = Wm1 + Wdot*dt
       mem = memm1 + memdot*dt

    end if

 
  end subroutine update_full_rk_adj



  subroutine adj_sl_ice(ice_load,sl_for,adj_load,Us_adj,ps_adj,Q_adj,K_adj,sl_adj,ice_adj)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load, sl_for
    type(meas), intent(in) :: adj_load
    complex(dpc), dimension((lmax+1)**2), intent(in) :: Us_adj,ps_adj
    complex(dpc), dimension(ngl,nphi), intent(in) :: Q_adj
    complex(dpc), intent(in) :: K_adj
    complex(dpc), dimension(ngl,nphi), intent(out) :: sl_adj,ice_adj

    real(dp) :: gsurf
    complex(dpc), dimension(ngl,nphi) :: ocean_func,CK,CUp,U_spat,p_spat,hspat

    integer(i4b) :: io,iphi,igl
    real(dp) :: long,lat

    sl_adj = 0.0_dp
    ice_adj = 0.0_dp

    gsurf = grav_node(ngll,nspec)

    ! Calculate C
    ocean_func = ocean_function(ice_load,sl_for)

    ! Calculate CK
    CK = ocean_func*K_adj

    ! Calcalate spatial U and p
    call fun_from_coefs(Us_adj,U_spat)
    call fun_from_coefs(ps_adj,p_spat)

    ! Calculate CUp
    CUp = ocean_func*(gsurf*U_spat + p_spat)

    ! Calculate spatial load
    call fun_from_coefs(adj_load%load_lm,hspat)

    ! ASSUMING ADJOINT LOAD IS ZERO OTHER THAN AT T=0

    ! Calculate adjoint sea level and ice
    sl_adj = Q_adj - CUp/gsurf + CK
    ice_adj = Q_adj - CUp/gsurf + CK - hspat/(roce*gsurf)

  end subroutine adj_sl_ice


  !=========================================!
  !             Sea level stuff             !
  !=========================================!

  function ocean_function(ice_load,sea_level)

    implicit none

    complex(dpc), dimension(ngl,nphi) :: ocean_function
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,sea_level

    integer(i4b) :: iphi,igl,io
    real(dp) :: test

    ocean_function = 0.0_dp

    do iphi = 1,nphi
       do igl = 1,ngl
          test = roce*sea_level(igl,iphi) - rice*ice_load(igl,iphi)
          if (test > 0.0_dp) then
             ocean_function(igl,iphi) = 1.0_dp
          else
             ocean_function(igl,iphi) = 0.0_dp
          end if
       end do
    end do


  end function ocean_function


  function B_function(ice_loaddot,ocean_func)

    implicit none

    complex(dpc), dimension(ngl,nphi) :: B_function
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_loaddot,ocean_func

    integer(i4b) :: iphi,igl
    real(dp) :: test

    B_function = 0.0_dp

    do iphi = 1,nphi
       do igl = 1,ngl
          B_function(igl,iphi) = (1.0_dp - ocean_func(igl,iphi))*ice_loaddot(igl,iphi)
       end do
    end do


  end function B_function


  subroutine extract_sl(file1,inv_data)

    implicit none

    character(len=256), intent(in) :: file1
    type(point_meas), dimension(:), intent(inout) :: inv_data

    integer(i4b) :: nloc,iloc,nt,it,t_suff,io,igl,iphi
    character(len=256) :: filet,fileo
    real(dp) :: long,lat,longr,latr,sl1,sl2
    

    nloc = size(inv_data)

    do iloc = 1,nloc
       nt = inv_data(iloc)%nt
       long = inv_data(iloc)%long
       lat = inv_data(iloc)%lat

       !print *, long,lat

       do it = 1,nt

          t_suff = floor(inv_data(iloc)%times(it))/100
          call string_cat_int(file1,t_suff/10,filet)
          if (mod(t_suff,10) == 0) then
             fileo = trim(filet) // '.0'
          else
             fileo = trim(filet) // '.5'
          end if

          open(newunit=io,file=fileo)
          do iphi = 1,nphi
             do igl = 1,ngl
                read(io,*) longr,latr,sl1,sl2

                if ((longr == long) .and. (latr == lat)) then
                   inv_data(iloc)%values(it) = sl1
                end if

             end do
             read(io,*)
          end do
          close(io)
       end do
    end do 


  end subroutine extract_sl


  subroutine extract_rsl(file1,inv_data)

    implicit none

    character(len=256), intent(in) :: file1
    type(point_meas), dimension(:), intent(inout) :: inv_data

    integer(i4b) :: nloc,iloc,nt,it,t_suff,io,igl,iphi,iop
    character(len=256) :: filet,fileo
    real(dp) :: long,lat,longr,latr,sl1,sl2,slp   

    nloc = size(inv_data)

    do iloc = 1,nloc
       nt = inv_data(iloc)%nt
       long = inv_data(iloc)%long
       lat = inv_data(iloc)%lat

       fileo = trim(file1) // trim('0.0')
       open(newunit=iop,file=fileo) 
       do iphi = 1,nphi
          do igl = 1,ngl
             read(iop,*) longr,latr,sl1,sl2

             if ((longr == long) .and. (latr == lat)) then
                slp = sl1
             end if

          end do
          read(iop,*)
       end do
       close(iop)

       ! Present day sea level
       inv_data(iloc)%values(nt) = slp

       do it = 1,nt-1

          t_suff = floor(inv_data(iloc)%times(it))/100
          call string_cat_int(file1,t_suff/10,filet)
          if (mod(t_suff,10) == 0) then
             fileo = trim(filet) // '.0'
          else
             fileo = trim(filet) // '.5'
          end if

          open(newunit=io,file=fileo)
          do iphi = 1,nphi
             do igl = 1,ngl
                read(io,*) longr,latr,sl1,sl2

                if ((longr == long) .and. (latr == lat)) then
                   inv_data(iloc)%values(it) = sl1 - slp
                end if

             end do
             read(io,*)
          end do
          close(io)
       end do
    end do 


  end subroutine extract_rsl


  subroutine extract_rsl_full(file1,nt,t2,inv_data)

    implicit none

    character(len=256), intent(in) :: file1
    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: t2
    type(point_meas), dimension(:), intent(inout) :: inv_data

    integer(i4b) :: nloc,iloc,it,t_suff,igl,iphi,iop,itn,igl1,igl2,iphi1,iphi2,it1,it2,ntl
    character(len=256) :: filet,fileo
    real(dp) :: long,lat,sl1,sl2,slp,xx,yy,f111,f112,f121,f122,f211,f212,f221,f222,dt,tt,time
    real(dp), dimension(ngl,nphi,nt) :: sl_in

    dt = t2/(nt-1)

    fileo = trim(file1)
    open(newunit=iop,file=fileo,action='read') 
    do it = 1,nt
       do iphi = 1,nphi
          do igl = 1,ngl
             read(iop,*) itn,long,lat,xx,yy,sl1,sl2
             sl_in(igl,iphi,it) = sl1
          end do
          read(iop,*)
       end do
       read(iop,*)
    end do
    close(iop)

    nloc = size(inv_data)

    do iloc = 1,nloc
       ntl = inv_data(iloc)%nt
       long = inv_data(iloc)%long
       lat = inv_data(iloc)%lat

       call get_igl_iphi(lat,long,igl1,iphi1,xx,yy)       
       igl2 = igl1 + 1
       iphi2 = iphi1 + 1

       ! Present day sea level
       f111 = sl_in(igl1,iphi1,nt)
       f112 = sl_in(igl2,iphi1,nt)
       f121 = sl_in(igl1,iphi2,nt)
       f122 = sl_in(igl2,iphi2,nt)

       slp = f111 + (f121 - f111)*xx + (f112 - f111)*yy &
             + (f122 - f112 - f121 + f111)*xx*yy

       inv_data(iloc)%values(ntl) = slp

       do it = 1,ntl-1

          time = inv_data(iloc)%times(it)
          
          call get_it(nt,dt,time,it1,tt)
          it2 = it1 + 1

          f111 = sl_in(igl1,iphi1,it1)
          f112 = sl_in(igl2,iphi1,it1)
          f121 = sl_in(igl1,iphi2,it1)
          f122 = sl_in(igl2,iphi2,it1)
          f211 = sl_in(igl1,iphi1,it2)
          f212 = sl_in(igl2,iphi1,it2)
          f221 = sl_in(igl1,iphi2,it2)
          f222 = sl_in(igl2,iphi2,it2)

          sl1 = f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                + (f122 - f112 - f121 + f111)*xx*yy &
                + tt*(f211+ (f221 - f211)*xx + (f212 - f211)*yy &
                + (f222 - f212 - f221 + f211)*xx*yy &
                - (f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                + (f122 - f112 - f121 + f111)*xx*yy))

          inv_data(iloc)%values(it) = sl1 - slp          

       end do
    end do 


  end subroutine extract_rsl_full



  subroutine extract_grace(file1,inv_data)

    implicit none

    character(len=256), intent(in) :: file1
    type(rsh_meas), intent(inout) :: inv_data

    integer(i4b) :: nt,it,t_suff,io,l,m,lt,mt,lmcur
    character(len=256) :: filet,fileo
    real(dp) :: pot
    

    nt = inv_data%nt
    
    do it = 1,nt

       t_suff = floor(inv_data%times(it))/100
       call string_cat_int(file1,t_suff/10,filet)
       if (mod(t_suff,10) == 0) then
          fileo = trim(filet) // '.0'
       else
          fileo = trim(filet) // '.5'
       end if

       open(newunit=io,file=fileo)
       do l = 2,lmax
          do m = -l,l
             lmcur = l**2 + l + m + 1
             read(io,*) lt,mt,pot
             inv_data%values(lmcur,it) = pot
          end do
       end do
       close(io)
    end do

  end subroutine extract_grace


  subroutine write_data(t2,tp,nout,file1,type,form,sl,UVp)

    implicit none
    
    integer(i4b), intent(in) :: nout,type,form
    real(dp), intent(in) :: t2,tp
    character(len=256), intent(in) :: file1
    complex(dpc), dimension(ngl,nphi,nout), intent(in), optional :: sl
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nout), intent(in), optional :: UVp

    ! type
    ! 1 = sea level
    ! 2 = vertical displacement
    ! 3 = gravitational potential

    ! form
    ! 1 = spatial
    ! 2 = complex spherical harmonic
    ! 3 = real spherical harmonic

    integer(i4b) :: it,num1,num2,io,iphi,igl,l,m,lmc
    real(dp) :: dt_out,t2y,tpy,lat,long,xx,yy,out1,out2
    character(len=256) :: file2,file3

    t2y = t2*t_norm*sec2yr
    tpy = tp*t_norm*sec2yr

    do it = 1,nout
       dt_out = tp - t2*(it - 1)/(nout - 1)
       num1 = floor(dt_out/1000)
       call string_cat_int(file1,num1,file2)
       file2 = trim(file3) // trim('.')
       num2 = floor((dt_out - 1000*num1)/100)
       call string_cat_int(file3,num2,file2)

       open(newunit=io,file=file2)
       if (form == 1) then
          do iphi = 1,nphi
             long = aphi(iphi)*rad2deg - 180.0_dp
             do igl = 1,ngl
                lat = (pio2_d - tgl(igl))*rad2deg
                call get_xy(lat,long,xx,yy)
                if (type == 1) then
                   out1 = real(sl(igl,iphi,it))*r_norm
                   out2 = (real(sl(igl,iphi,it)) - real(sl(igl,iphi,1)))*r_norm
                end if
                write(io,*) long,lat,xx,yy,out1,out2
             end do
             write(io,*)
          end do

       else if (form == 2) then
          do l = 0,lmax
             do m = -l,l
                lmc = l**2 + l + m + 1
                if (type == 2) then
                   out1 = real(UVp(lmc,pa_surf(l)+1,it))*r_norm
                   out2 = imag(UVp(lmc,pa_surf(l)+1,it))*r_norm
                else if (type == 3) then
                   out1 = real(UVp(lmc,pa_surf(l),it))*pot_norm
                   out2 = imag(UVp(lmc,pa_surf(l),it))*pot_norm
                end if
                write(io,*) l,m,out1,out2
             end do
          end do

       else if (form == 3) then

          do l = 0,lmax
             do m = -l,l
                if (m < 0) then
                   lmc = l**2 + l - m + 1
                   if (type == 2) then
                      out1 = sqrt(2.0_dp)*real(UVp(lmc,pa_surf(l)+1,it))*r_norm
                   else if (type == 3) then
                      out1 = sqrt(2.0_dp)*real(UVp(lmc,pa_surf(l),it))*pot_norm
                   end if
                else if (m == 0) then
                   if (type == 2) then
                      out1 = real(UVp(lmc,pa_surf(l)+1,it))*r_norm
                   else if (type == 3) then
                      out1 = real(UVp(lmc,pa_surf(l),it))*pot_norm
                   end if
                else
                   if (type == 2) then
                      out1 = -sqrt(2.0_dp)*imag(UVp(lmc,pa_surf(l)+1,it))*r_norm
                   else if (type == 3) then
                      out1 = -sqrt(2.0_dp)*imag(UVp(lmc,pa_surf(l),it))*pot_norm
                   end if
                end if
                write(io,*) l,m,out1,out2
             end do
          end do

       else
          print *, "No such form"
          return
       end if

       close(io)
    end do


  end subroutine write_data


  subroutine get_igl_iphi(lat,long,igl1,iphi1,xx,yy)

    implicit none

    real(dp), intent(in) :: lat,long
    integer(i4b), intent(out) :: igl1,iphi1
    real(dp), intent(out) :: xx,yy

    integer(i4b) :: igl,iphi
    real(dp) :: lati,longi

    do iphi = 1,nphi
       longi = aphi(iphi)*rad2deg - 180.0_dp
       if (longi > long) then
          iphi1 = iphi - 1
          xx = (long - (aphi(iphi1)*rad2deg - 180.0_dp))/(longi - aphi(iphi1)*rad2deg)
          exit
       end if
    end do

    do igl = 1,ngl
       lati = (pio2_d - tgl(igl))*rad2deg
       if (lati > lat) then
          igl1 = igl - 1
          yy = (lat - (pio2_d - tgl(igl1))*rad2deg)/(tgl(igl1)*rad2deg - lati)
          exit
       end if
    end do
          

  end subroutine get_igl_iphi


  subroutine get_it(nt,dt,time,it1,tt)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt,time
    integer(i4b), intent(out) :: it1
    real(dp), intent(out) :: tt

    integer(i4b) :: it
    real(dp) :: itime

    it1 = 0
    tt = 0.0_dp

    ! REMEMBER TIME IS BEFORE PRESENT
    do it = 1,nt
       itime = 21000.0_dp*yr2sec/t_norm - (it-1)*dt
       if (itime < time) then
          it1 = it - 1
          itime = 21000.0_dp*yr2sec/t_norm - (it1-1)*dt
          tt = (itime - time)/dt
          exit
       end if
    end do

    if (it1 == 0) then
       it1 = nt - 1
       tt = 1.0_dp
    end if

  end subroutine get_it


  !subroutine UVp_l_0(nrows,bvec,UVpcur)

  !  use nrtype

  !  implicit none

  !  integer(i4b), intent(in) :: nrows
  !  complex(dpc), dimension(1:nrows), intent(in) :: bvec
  !  complex(dpc), dimension(1:nrows), intent(out) :: UVpcur

  !  integer(i4b) :: irow
  !  complex(dpc) :: fi,ui

  !  UVpcur = 0.0_dp

  !  do irow = 1,nrows
        
  !     fi = dot_product(conjg(bvec(:)),Umat(:,irow))
  !     ui = fi/s_vals(irow)

  !     UVpcur(:) = UVpcur(:) + ui*(Vmat(irow,:))

  !  end do

  !end subroutine UVp_l_0




end module module_sl
