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
          call string_cat_int('/home/oc251/raid/data/ice5g/ascii/ice5g_v1.2_',22-it,file1)
       else if (it < 20) then
          call string_cat_int('/home/oc251/raid/data/ice5g/ascii/ice5g_v1.2_',(nt_ice5g-it)/2,file1)
       else
          call string_cat_int('/home/oc251/raid/data/ice5g/ascii/ice5g_v1.2_0',(nt_ice5g-it)/2,file1)
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



  subroutine memdot_bvec_sl(memm1,UVpm1,dUm1,dVm1,memdot,bUVp,b_lm,ocean_lm)
    
    use nrtype
    implicit none

    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(in) :: dUm1,dVm1
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: bUVp
    complex(dpc), dimension((lmax+1)**2), intent(in), optional :: b_lm,ocean_lm

    integer(i4b) :: inode,ispec,ilayer,pa,ua,va,l,m,ispec1,lmtot,lmcur,i,jnode,ub,vb
    real(dp) :: zeta2,zeta,zeta2m2,rl,rr,rmu,gsurf,rsurf,rsi,tmp
    complex(dpc) :: bint,area,comb1,comb2,comb3

    ! A(1) = A(-,-)
    ! A(2) = A(-,0)
    ! A(3) = A(-,+) = 0.5*A(0,0)
    ! A(4) = A(0,+)
    ! A(5) = A(+,+)

    rsurf = r_node(ngll,nspec)
    gsurf = grav_node(ngll,nspec)

    memdot = 0.0_dp
    bUVp = 0.0_dp

    do ilayer = 1,nsect

       ! No terms in elastic or fluid regions
       if (sect_ind(3,ilayer) /=2) cycle


       !$OMP PARALLEL DO PRIVATE(pa,ua,va,ispec1,lmtot,lmcur,i,zeta2,zeta,zeta2m2), &
       !$OMP& PRIVATE(rl,rr,rmu,rsi,ub,vb,tmp,comb1,comb2,comb3,inode,jnode,l,m), &
       !$OMP& SHARED(lmax,memm1,UVpm1,dUm1,dVm1,memdot,bUVp,b_lm,ocean_lm)
       do ispec = sect_ind(1,ilayer), sect_ind(2,ilayer)

          do inode = 1,ngll
             rr = r_node(inode,ispec)
             rmu = mui_node(1,inode,ispec)
             rsi = si_node(1,inode,ispec)

             if (rr == 0.0_dp) rr = r_node(2,1)

             do l = 0,lmax
                lmtot = l**2 + 1
                call get_depth(l,rl,ispec1)
                if (ispec1 > ispec) then
                   exit
                end if

                ua = gvn_ssg(inode,ispec,ispec1,2)
                va = ua + 1

                zeta2 = l*(l+1)
                zeta = sqrt(zeta2)
                zeta2m2 = zeta2 - 2.0_dp

                do m = -l,l  
                   lmcur = lmtot + l + m

                   comb1 = rsi*(UVpm1(lmcur,va) - memm1(1,lmcur,1,inode,ispec))
                   comb2 = rsi*(rr*dUm1(lmcur,inode,ispec) - UVpm1(lmcur,ua) &
                           + 0.5_dp*zeta2*UVpm1(lmcur,va) - memm1(3,lmcur,1,inode,ispec))
                   comb3 = rsi*(rr*dVm1(lmcur,inode,ispec) - UVpm1(lmcur,va) &
                           + UVpm1(lmcur,ua) - memm1(4,lmcur,1,inode,ispec))

                   memdot(1,lmcur,1,inode,ispec) = comb1
                   memdot(3,lmcur,1,inode,ispec) = comb2
                   memdot(4,lmcur,1,inode,ispec) = comb3


                   !if ((ispec == nspec/2) .and. (inode == 1) .and. (l==2)) then
                      !print *, m,uu,vv,dU,dV,mm,rm,sm
                   !   print *, m,memdot(1,lmcur,1,inode,ispec),memdot(3,lmcur,1,inode,ispec)
                   !end if

                   bUVp(lmcur,ua) = bUVp(lmcur,ua) + zeta2*jac_element(ispec)*wgll(inode) &
                                                     *rmu*comb3

                   bUVp(lmcur,va) = bUVp(lmcur,va) + zeta2*jac_element(ispec)*wgll(inode)*rmu &
                                                     *((2.0_dp/3.0_dp)*comb2 + zeta2m2*comb1)



                   ! Derivative terms
                   do jnode = 1,ngll

                      tmp = hprime(inode,jnode)*rr
                      if (jnode == inode) tmp = tmp - jac_element(ispec)

                      if (jnode /= 1) then
                         ub = ub + 3
                      else
                         ub = gvn_ssg(jnode,ispec,ispec1,2)
                      end if
                      !ub = gvn_ssg(jnode,ispec,ispec1,2)
                      vb = ub + 1

                      bUVp(lmcur,ub) = bUVp(lmcur,ub) + (4.0_dp/3.0_dp)*wgll(inode)*tmp*rmu*comb2
                      bUVp(lmcur,vb) = bUVp(lmcur,vb) + zeta2*wgll(inode)*tmp*rmu*comb3
                      
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



  subroutine time_derivs_full_sl(ice_load,ice_loaddot,sea_level,memm1,UVpm1, &
                                sldot,memdot,UVpdot)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,ice_loaddot,sea_level
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVpdot
    complex(dpc), dimension(ngl,nphi), intent(out) :: sldot

    real(dp) :: area
    complex(dpc), dimension(ngl,nphi) :: ocean_func,B_func,U_func,p_func,CUp
    complex(dpc), dimension((lmax+1)**2) :: ocean_func_lm,B_func_lm,CUp_lm,U_lm,p_lm
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: bUVp

    integer(i4b) :: l,ispec1,lmcur,lmact,kd,kl,ku,ldab,nrows,ua,va,pa,info,m,ispec11,ilayer, &
                    ispec,inode,i,iter,io,igl,iphi
    real(dp) :: rl,rr,gsurf,ratio,rsurf
    real(dp), dimension(ntau) :: rsi
    complex(dpc) :: int,bint
    complex(dpc), dimension(:,:), allocatable :: A1u,A1u_trans
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpcur,UVpnew,UVp0

    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV

    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: UVp_trans

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    kl = 3*(ngll-1) + 2
    ku = 3*(ngll-1) + 2
    kd = ku

    ! Calculate C
    ocean_func = ocean_function(ice_load,sea_level)

    ! Calculate spherical harmonic coefficients of C
    call coefs_from_fun_parr(lmax,ocean_func,ocean_func_lm)

    ! Calculate A
    area = sqrt(fourpi_d)*ocean_func_lm(1)*(rsurf**2)

    ! Calculate B = (1-C)*Idot
    B_func = B_function(ice_loaddot,ocean_func)

    ! Calculate spherical harmonic coefficients of B
    call coefs_from_fun_parr(lmax,B_func,B_func_lm)

    bint = B_func_lm(1)*sqrt(fourpi_d)*(rsurf**2)

    ! Calculate stresses
    call calculate_dUV_full(lmax,UVpm1,dU,dV)

    !-------------------------------------!
    !      Calculate memdot and RHSs      !
    !-------------------------------------!

    call memdot_bvec_sl(memm1,UVpm1,dU,dV,memdot,bUVp,B_func_lm,ocean_func_lm)

    print *, "Calculating u0"

    do l = 1,lmax

       call get_depth(l,rl,ispec1)
       
       !---------------------------------------!
       !      Calculate spheroidal modes       !
       !---------------------------------------!
    
       nrows = gvn_ssg(ngll,nspec,ispec1,3)

       UVp_trans(1:nrows,-l:l) = transpose(bUVp(l**2+1:(l+1)**2,1:nrows))

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

             p_lm(lmcur+l+m) = UVpcur(lmcur+l+m,pa)
             U_lm(lmcur+l+m) = UVpcur(lmcur+l+m,ua)

          end do
          lmcur = lmcur + 2*l + 1
       end do

       p_func = 0.0_dp
       U_func = 0.0_dp  
      
       call fun_from_coefs_parr(lmax,p_lm,p_func)
       call fun_from_coefs_parr(lmax,U_lm,U_func)  

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func

       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun_parr(lmax,CUp,CUp_lm)

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
      
    call fun_from_coefs_parr(lmax,p_lm,p_func)
    call fun_from_coefs_parr(lmax,U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func

    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(lmax,CUp,CUp_lm(1))
    
    sldot = sqrt(fourpi_d)*(rsurf**2)*((CUp_lm(1)/gsurf) - (rice*B_func_lm(1)/roce))/area &
            - U_func - p_func/gsurf

  end subroutine time_derivs_full_sl


  subroutine update_full_rk_sl(dt,ice_load,ice_loaddot,slm1,memm1,UVpm1,sl,mem,UVp)
  ! routine to calculate the displacements from the previous displacements
  ! using Euler time-stepping  

    use nrtype

    implicit none

    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,ice_loaddot
    complex(dpc), dimension(ngl,nphi), intent(in) :: slm1
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension(ngl,nphi), intent(out) :: sl
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(out) :: mem
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp

    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec) :: memdot,memtmp
    complex(dpc), dimension(ngl,nphi) :: sldot,sltmp,ice_loadtmp

    integer(i4b) :: ua

    UVpdot = 0.0_dp
    memdot = 0.0_dp

    ! Calculate time derivatives
    call time_derivs_full_sl(ice_load,ice_loaddot,slm1,memm1,UVpm1,sldot,memdot,UVpdot)

    UVptmp = UVpm1 + UVpdot*dt/2.0_dp

    memtmp = memm1 + memdot*dt/2.0_dp

    sltmp = slm1 + sldot*dt/2.0_dp

    ice_loadtmp = ice_load + ice_loaddot*dt/2.0_dp

    ! Calculate time derivatives
    call time_derivs_full_sl(ice_loadtmp,ice_loaddot,sltmp,memtmp,UVptmp,sldot,memdot,UVpdot)

    UVp = UVpm1 + UVpdot*dt

    mem = memm1 + memdot*dt

    sl = slm1 + sldot*dt
  
  end subroutine update_full_rk_sl


  subroutine forward_full_sl(nt,dt,nout,ice_data,sl_0,sl_out,UVp_out,mem_out)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_0
    complex(dpc), dimension(ngl,nphi,nout), intent(out) :: sl_out
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nout), intent(out), optional :: UVp_out
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec,nout), intent(out), optional :: mem_out

    integer(i4b) :: it,it_out,iout
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec,2) :: mem
    complex(dpc), dimension(ngl,nphi,2) :: sl_f
    complex(dpc), dimension(ngl,nphi) :: ice_loaddot,ice_load

    UVp = 0.0_dp
    mem = 0.0_dp
    sl_f = 0.0_dp
    ice_loaddot = 0.0_dp
    ice_load = 0.0_dp

    sl_out = 0.0_dp

    if (present(UVp_out)) UVp_out = 0.0_dp
    if (present(mem_out)) mem_out = 0.0_dp

    ! Calculate it_out
    it_out = (nt-1)/(nout-1)

    sl_f(:,:,1) = sl_0
    sl_out(:,:,1) = sl_0

    iout = 1
    do it = 2,nt
       print *, "Forward time step: ",it," of ",nt 
       call time_point_interp(it-1,dt,ice_data,ice_load,ice_loaddot)
       call update_full_rk_sl(dt,ice_load,ice_loaddot,sl_f(:,:,1),mem(:,:,:,:,:,1), &
                              UVp(:,:,1),sl_f(:,:,2),mem(:,:,:,:,:,2),UVp(:,:,2))
       
       if (mod(it-1,it_out) == 0) then
          iout = iout + 1
          sl_out(:,:,iout) = sl_f(:,:,2)

          if (present(UVp_out)) UVp_out(:,:,iout) = UVp(:,:,2)
          if (present(mem_out)) mem_out(:,:,:,:,:,iout) = mem(:,:,:,:,:,2)

       end if

       mem(:,:,:,:,:,1) = mem(:,:,:,:,:,2)
       UVp(:,:,1) = UVp(:,:,2)
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


  subroutine real_sph_load(l,m,load_lm)

    integer(i4b), intent(in) :: l,m
    complex(dpc), dimension((lmax+1)**2), intent(out) :: load_lm

    integer(i4b) :: lmcur

    load_lm = 0.0_dp

    if (m < 0) then

       lmcur = l**2 + l + m + 1
       load_lm(lmcur) = ((-1)**m)/sqrt(2.0_dp)

       lmcur = l**2 + l - m + 1
       load_lm(lmcur) = 1.0_dp/sqrt(2.0_dp)


    else if (m == 0) then

       lmcur = l**2 + l + m + 1
       load_lm(lmcur) = 1.0_dp

    else

       lmcur = l**2 + l + m + 1
       load_lm(lmcur) = -ii/sqrt(2.0_dp)

       lmcur = l**2 + l - m + 1
       load_lm(lmcur) = ii*((-1)**m)/sqrt(2.0_dp)

    end if


  end subroutine real_sph_load


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
    call coefs_from_fun(lmax,ocean_func,ocean_func_lm)

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
      
       call fun_from_coefs(lmax,p_lm,p_func)
       call fun_from_coefs(lmax,U_lm,U_func)  

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func

       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun(lmax,CUp,CUp_lm)

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
      
    call fun_from_coefs(lmax,p_lm,p_func)
    call fun_from_coefs(lmax,U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func
    
    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(lmax,CUp,CUp_lm(1))

    if (load%type == 1) then

       call fun_from_coefs(lmax,load%load_lm,hall_spat)
    
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


  subroutine time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,memdot,UVpdot,Qdot,Kdot)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,sl_for
    complex(dpc), dimension(5,(lmax+1)**2,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(out) :: memdot
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVpdot
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: Qdot
    complex(dpc), intent(out), optional :: Kdot

    real(dp) :: area
    complex(dpc), dimension(ngl,nphi) :: U_func,p_func,CUp,ocean_func
    complex(dpc), dimension((lmax+1)**2) :: ocean_func_lm,CUp_lm,U_lm,p_lm
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: bUVp

    integer(i4b) :: l,ispec1,lmcur,kl,ku,kd,ldab,nrows,ua,va,pa,info,m,ispec11,ilayer, &
                    ispec,inode,i,iter,io,igl,iphi
    real(dp) :: zeta2,zeta,c,rl,rr,gsurf,long,lat,ratio,rsurf
    real(dp), dimension(ntau) :: rsi
    complex(dpc) :: int
    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: A1u_trans
    complex(dpc), dimension(-lmax:lmax,nglob_ssg) :: A1u
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpcur,UVpnew,UVp0

    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV

    complex(dpc), dimension(nglob_ssg,-lmax:lmax) :: UVp_trans

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    ! Calculate C
    ocean_func = ocean_function(ice_load,sl_for)

    ! Calculate spherical harmonic coefficients of C
    call coefs_from_fun(lmax,ocean_func,ocean_func_lm)

    ! Calculate A
    area = sqrt(fourpi_d)*ocean_func_lm(1)*(rsurf**2)

    ! Calculate stresses
    call calculate_dUV_full(lmax,UVpm1,dU,dV)

    !-------------------------------------!
    !      Calculate memdot and RHSs      !
    !-------------------------------------!

    bUVp = 0.0_dp

    memdot = 0.0_dp

    call memdot_bvec_sl(memm1,UVpm1,dU,dV,memdot,bUVp)

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
      
       call fun_from_coefs(lmax,p_lm,p_func)
       call fun_from_coefs(lmax,U_lm,U_func) 

       ! Construct U,P,C product
       CUp = 0.0_dp
       CUp = (gsurf*U_func + p_func)*ocean_func


       ! Calculate spherical harmonic coefficients of CUp
       CUp_lm = 0.0_dp
       call coefs_from_fun(lmax,CUp,CUp_lm)

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
      
    call fun_from_coefs(lmax,p_lm,p_func)
    call fun_from_coefs(lmax,U_lm,U_func)  

    ! Construct U,P,C product
    CUp = 0.0_dp
    CUp = (gsurf*U_func + p_func)*ocean_func

    ! Calculate spherical harmonic coefficients of CUp
    CUp_lm = 0.0_dp
    call zero_coef_from_fun(lmax,CUp,CUp_lm(1))

    ! AGAIN ASSUMING NO TIME DEPENDENCE OF ADJOINT LOAD

    Kdot = sqrt(fourpi_d)*(rsurf**2)*CUp_lm(1)/(gsurf*area)

    Qdot = CUp/gsurf - ocean_func*sqrt(fourpi_d)*(rsurf**2)*CUp_lm(1)/(gsurf*area)
                                      

  end subroutine time_derivs_full_adj


  subroutine update_full_rk_adj(dt,ice_load,sl_for,memm1,UVpm1,mem,UVp,Qm1,Km1,Qt,Kt)
  ! routine to calculate the displacements from the previous displacements
  ! using Euler time-stepping  

    use nrtype

    implicit none

    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_for
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(in) :: memm1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVpm1
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(out) :: mem
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(out) :: UVp
    complex(dpc), dimension(ngl,nphi), intent(in), optional :: Qm1
    complex(dpc), intent(in), optional :: Km1
    complex(dpc), dimension(ngl,nphi), intent(out), optional :: Qt
    complex(dpc), intent(out),optional :: Kt

    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension(ngl,nphi) :: Qdot
    complex(dpc) :: Kdot
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec) :: memdot,memtmp
    complex(dpc), dimension(5,(lmax+1)**2) :: memcoefs

    integer(i4b) :: io,igl,iphi


    if (present(Qm1)) then

       call time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,memdot,UVpdot)

       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp

       call time_derivs_full_adj(ice_load,sl_for,memtmp,UVptmp,memdot,UVpdot,Qdot,Kdot) 

       UVp = UVpm1 + UVpdot*dt
       mem = memm1 + memdot*dt
       Qt = Qm1 + Qdot*dt
       Kt = Km1 + Kdot*dt

    else

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_for,memm1,UVpm1,memdot,UVpdot)
       UVptmp = UVpm1 + UVpdot*dt/2.0_dp
       memtmp = memm1 + memdot*dt/2.0_dp

       call time_derivs_full_adj(ice_load,sl_for,memtmp,UVptmp,memdot,UVpdot)
       UVp = UVpm1 + UVpdot*dt
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
    call fun_from_coefs(lmax,Us_adj,U_spat)
    call fun_from_coefs(lmax,ps_adj,p_spat)

    ! Calculate CUp
    CUp = ocean_func*(gsurf*U_spat + p_spat)

    ! Calculate spatial load
    call fun_from_coefs(lmax,adj_load%load_lm,hspat)

    ! ASSUMING ADJOINT LOAD IS ZERO OTHER THAN AT T=0

    ! Calculate adjoint sea level and ice
    sl_adj = Q_adj - CUp/gsurf + CK
    ice_adj = Q_adj - CUp/gsurf + CK - hspat/(roce*gsurf)

    open(newunit=io,file='sl_init_kern.out')
    do iphi = 1,nphi
       do igl = 1,ngl
          write(io,*) iphi,igl,Q_adj(igl,iphi),sl_adj(igl,iphi),ice_adj(igl,iphi),U_spat(igl,iphi),p_spat(igl,iphi),hspat(igl,iphi)
       end do
       write(io,*)
    end do
    close(io)
    print *, K_adj

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
