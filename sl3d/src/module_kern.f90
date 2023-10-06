module module_kern

  use nrtype
  use module_util
  use module_mat
  use module_model
  use module_function
  use module_fourier
  use module_gia
  use module_sl

  implicit none


contains

  !===================================================!
  !   Adjoint and viscosity kernel full calculation   !
  !===================================================!
  
  subroutine adj_kern_visc_full(nt,dt,nout,ice_data,sl_fin,adj_load,UVp,W,mem,kern)

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_fin
    type(meas), intent(in) :: adj_load
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nout), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor,nout), intent(in) :: W
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,nout), intent(in) :: mem
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(out) :: kern

    integer(i4b) :: it,it1,it_out,iout,rem
    real(dp) :: fac
    complex(dpc), dimension(ngl,nphi) :: ice_load,sltmp,ice_tmp
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor,2) :: W_adj
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,2) :: mem_adj
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot,Wtmp
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot,memtmp

    UVp_adj = 0.0_dp
    W_adj = 0.0_dp
    mem_adj = 0.0_dp
    kern = 0.0_dp

    it_out = (nt-1)/(nout-1)

    call time_point_interp(nt,dt,ice_data,ice_load)

    call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1))

    call kern_visc(1,nt,dt,UVp(:,:,nout),W(:,:,nout),UVp_adj(:,:,1),W_adj(:,:,1), &
                   mem(:,:,:,:,:,:,nout),mem_adj(:,:,:,:,:,:,1),kern)

    iout = nout
    do it = 2,nt

       print *, "adjoint time step = ",it," of ",nt

       ! Interpolate to find rate of change of load
       ! We want load at previous time step and time reversed
       it1 = time_point(nt-it+2,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+2,dt,ice_data(:,:,it1:it1+1),ice_load)
       end if

       UVpdot = 0.0_dp
       Wdot = 0.0_dp
       memdot = 0.0_dp

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_fin(:,:,nt-it+2),mem_adj(:,:,:,:,:,:,1), &
                                 UVp_adj(:,:,1),W_adj(:,:,1),memdot,UVpdot,Wdot)

       UVptmp = UVp_adj(:,:,1) + UVpdot*dt/2.0_dp
       Wtmp = W_adj(:,:,1) + Wdot*dt/2.0_dp
       memtmp = mem_adj(:,:,:,:,:,:,1) + memdot*dt/2.0_dp
       sltmp = 0.5_dp*(sl_fin(:,:,nt-it+2) + sl_fin(:,:,nt-it+1))

       it1 = time_point(nt-it+1,dt)
       if (it1 == nt_ice5g) then
          ice_tmp = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+1,dt,ice_data(:,:,it1:it1+1),ice_tmp)
       end if
       ice_tmp = 0.5_dp*(ice_tmp + ice_load)

       call time_derivs_full_adj(ice_tmp,sltmp,memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot)
       UVp_adj(:,:,2) = UVp_adj(:,:,1) + UVpdot*dt
       W_adj(:,:,2) = W_adj(:,:,1) + Wdot*dt
       mem_adj(:,:,:,:,:,:,2) = mem_adj(:,:,:,:,:,:,1) + memdot*dt

       rem = mod(it-1,it_out)
       if (rem == 0) then
          iout = iout - 1
          UVptmp = UVp(:,:,iout)
          Wtmp = W(:,:,iout)
          memtmp = mem(:,:,:,:,:,:,iout)
       else
          fac = (it_out - rem)/it_out
          UVptmp = UVp(:,:,iout-1) + fac*(UVp(:,:,iout) - UVp(:,:,iout-1))
          Wtmp = W(:,:,iout-1) + fac*(W(:,:,iout) - W(:,:,iout-1))
          memtmp = mem(:,:,:,:,:,:,iout-1) + fac*(mem(:,:,:,:,:,:,iout) - mem(:,:,:,:,:,:,iout-1))
       end if

       call kern_visc(it,nt,dt,UVptmp,Wtmp,UVp_adj(:,:,2),W_adj(:,:,2), &
                      memtmp,mem_adj(:,:,:,:,:,:,2),kern)

       UVp_adj(:,:,1) = UVp_adj(:,:,2)
       W_adj(:,:,1) = W_adj(:,:,2)
       mem_adj(:,:,:,:,:,:,1) = mem_adj(:,:,:,:,:,:,2)
       
    end do



  end subroutine adj_kern_visc_full


  subroutine adj_kern_visc_misfit(nt,dt,nout,ice_data,sl_fin,UVp,W,mem,kern,point_data,rsh_data)

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_fin
    type(point_meas), dimension(:), intent(in), optional :: point_data
    type(rsh_meas), dimension(:), intent(in), optional :: rsh_data
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nt), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor,nt), intent(in) :: W
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,nt), intent(in) :: mem
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(out) :: kern

    integer(i4b) :: it,it1,io,igl,iphi
    type(meas) :: adj_load
    logical(lgt) :: load_flag
    character(len=256) :: file1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: delta_UVpa
    complex(dpc), dimension(ngl,nphi) :: ice_load,sltmp,ice_tmp,load_spat
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor,2) :: W_adj
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,2) :: mem_adj
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot,Wtmp
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot,memtmp

    UVp_adj = 0.0_dp
    W_adj = 0.0_dp
    mem_adj = 0.0_dp
    kern = 0.0_dp

    it = 1

    call time_point_interp(nt,dt,ice_data,ice_load)

    ! Calculate load at current time step
    allocate(adj_load%load_lm((lmax+1)**2))

    if (present(point_data)) then
       call adj_load_misfit_point(nt,it,dt,sl_fin(:,:,nt),UVp(:,:,nt),point_data, &
                                  adj_load,load_flag)
    else if (present(rsh_data)) then
       call adj_load_misfit_rsh(nt,it,dt,sl_fin(:,:,nt),UVp(:,:,nt),rsh_data, &
                                adj_load,load_flag)
    else 
       print *, "No data"
       return
    end if

    call fun_from_coefs(adj_load%load_lm,load_spat)


    open(newunit=io,file='adjload.1')
    do iphi = 1,nphi
       do igl = 1,ngl
          write(io,*) iphi,igl,real(load_spat(igl,iphi)),imag(load_spat(igl,iphi))
       end do
       write(io,*)
    end do
    close(io)

    call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1))

    call kern_visc(1,nt,dt,UVp(:,:,nt),W(:,:,nt),UVp_adj(:,:,1),W_adj(:,:,1), &
                   mem(:,:,:,:,:,:,nt),mem_adj(:,:,:,:,:,:,1),kern)

    do it = 2,nt

       print *, "adjoint time step = ",it," of ",nt

       ! Interpolate to find rate of change of load
       ! We want load at previous time step and time reversed
       it1 = time_point(nt-it+2,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+2,dt,ice_data(:,:,it1:it1+1),ice_load)
       end if

       UVpdot = 0.0_dp
       Wdot = 0.0_dp
       memdot = 0.0_dp

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_fin(:,:,nt-it+2),mem_adj(:,:,:,:,:,:,1), &
                                 UVp_adj(:,:,1),W_adj(:,:,1),memdot,UVpdot,Wdot)

       UVptmp = UVp_adj(:,:,1) + UVpdot*dt/2.0_dp
       Wtmp = W_adj(:,:,1) + Wdot*dt/2.0_dp
       memtmp = mem_adj(:,:,:,:,:,:,1) + memdot*dt/2.0_dp
       sltmp = 0.5_dp*(sl_fin(:,:,nt-it+2) + sl_fin(:,:,nt-it+1))

       it1 = time_point(nt-it+1,dt)
       if (it1 == nt_ice5g) then
          ice_tmp = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+1,dt,ice_data(:,:,it1:it1+1),ice_tmp)
       end if
       ice_tmp = 0.5_dp*(ice_tmp + ice_load)

       call time_derivs_full_adj(ice_tmp,sltmp,memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot)
       UVp_adj(:,:,2) = UVp_adj(:,:,1) + UVpdot*dt
       W_adj(:,:,2) = W_adj(:,:,1) + Wdot*dt
       mem_adj(:,:,:,:,:,:,2) = mem_adj(:,:,:,:,:,:,1) + memdot*dt

       ! Check to see if an elastic change occurs
       if (present(point_data)) then
          call adj_load_misfit_point(nt,it,dt,sl_fin(:,:,nt-it+1),UVp(:,:,nt-it+1),point_data, &
                                     adj_load,load_flag)
       else if (present(rsh_data)) then
          call adj_load_misfit_rsh(nt,it,dt,sl_fin(:,:,nt-it+1),UVp(:,:,nt-it+1),rsh_data, &
                                   adj_load,load_flag)
       else 
          print *, "No data"
          return
       end if
       
       if (load_flag) then
          call fun_from_coefs(adj_load%load_lm,load_spat)
          call string_cat_int('adjload.',it,file1)
          open(newunit=io,file=trim(file1))
          do iphi = 1,nphi
             do igl = 1,ngl
                write(io,*) iphi,igl,real(load_spat(igl,iphi)),imag(load_spat(igl,iphi))
             end do
             write(io,*)
          end do
          close(io)

          call stat_adj(ice_load,sl_fin(:,:,nt-it+1),adj_load,delta_UVpa(:,:))
          UVp_adj(:,:,2) = UVp_adj(:,:,2) + delta_UVpa(:,:)
       end if

       call kern_visc(it,nt,dt,UVp(:,:,nt-it+1),W(:,:,nt-it+1),UVp_adj(:,:,2),W_adj(:,:,2), &
                      mem(:,:,:,:,:,:,nt-it+1),mem_adj(:,:,:,:,:,:,2),kern)

       UVp_adj(:,:,1) = UVp_adj(:,:,2)
       W_adj(:,:,1) = W_adj(:,:,2)
       mem_adj(:,:,:,:,:,:,1) = mem_adj(:,:,:,:,:,:,2)
       
    end do


  end subroutine adj_kern_visc_misfit


  subroutine adj_kern_visc_misfit_rel(nt,dt,nout,ice_data,sl_fin,UVp,W,mem,kern, &
                                      point_data,rsh_data)

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_fin
    type(point_meas), dimension(:), intent(in), optional :: point_data
    type(rsh_meas), dimension(:), intent(in), optional :: rsh_data
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nout), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor,nout), intent(in) :: W
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,nout), intent(in) :: mem
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(out) :: kern

    integer(i4b) :: it,it1,io,igl,iphi,it_out,iout
    real(dp) :: rem,fac
    type(meas) :: adj_load
    logical(lgt) :: load_flag
    character(len=256) :: file1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: delta_UVpa
    complex(dpc), dimension(ngl,nphi) :: ice_load,sltmp,ice_tmp,load_spat
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor,2) :: W_adj
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,2) :: mem_adj
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension((lmax+1)**2,nglob_tor) :: Wdot,Wtmp
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec) :: memdot,memtmp

    UVp_adj = 0.0_dp
    W_adj = 0.0_dp
    mem_adj = 0.0_dp
    kern = 0.0_dp

    it_out = (nt-1)/(nout-1)

    it = 1

    call time_point_interp(nt,dt,ice_data,ice_load)

    ! Calculate load at current time step
    allocate(adj_load%load_lm((lmax+1)**2))

    if (present(point_data)) then
       call adj_load_misfit_point_rel(nt,it,dt,sl_fin,UVp,point_data, &
                                      adj_load,load_flag)
    else if (present(rsh_data)) then
       call adj_load_misfit_rsh(nt,it,dt,sl_fin(:,:,nt),UVp(:,:,nout),rsh_data, &
                                adj_load,load_flag)
    else 
       print *, "No data"
       return
    end if

    call fun_from_coefs(adj_load%load_lm,load_spat)


    open(newunit=io,file='adjload.1')
    do iphi = 1,nphi
       do igl = 1,ngl
          write(io,*) iphi,igl,real(load_spat(igl,iphi)),imag(load_spat(igl,iphi))
       end do
       write(io,*)
    end do
    close(io)

    call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1))

    call kern_visc(1,nt,dt,UVp(:,:,nout),W(:,:,nout),UVp_adj(:,:,1),W_adj(:,:,1), &
                   mem(:,:,:,:,:,:,nout),mem_adj(:,:,:,:,:,:,1),kern)

    iout = nout
    do it = 2,nt

       print *, "adjoint time step = ",it," of ",nt

       ! Interpolate to find rate of change of load
       ! We want load at previous time step and time reversed
       it1 = time_point(nt-it+2,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+2,dt,ice_data(:,:,it1:it1+1),ice_load)
       end if

       UVpdot = 0.0_dp
       Wdot = 0.0_dp
       memdot = 0.0_dp

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_fin(:,:,nt-it+2),mem_adj(:,:,:,:,:,:,1), &
                                 UVp_adj(:,:,1),W_adj(:,:,1),memdot,UVpdot,Wdot)

       UVptmp = UVp_adj(:,:,1) + UVpdot*dt/2.0_dp
       Wtmp = W_adj(:,:,1) + Wdot*dt/2.0_dp
       memtmp = mem_adj(:,:,:,:,:,:,1) + memdot*dt/2.0_dp
       sltmp = 0.5_dp*(sl_fin(:,:,nt-it+2) + sl_fin(:,:,nt-it+1))

       it1 = time_point(nt-it+1,dt)
       if (it1 == nt_ice5g) then
          ice_tmp = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+1,dt,ice_data(:,:,it1:it1+1),ice_tmp)
       end if
       ice_tmp = 0.5_dp*(ice_tmp + ice_load)

       call time_derivs_full_adj(ice_tmp,sltmp,memtmp,UVptmp,Wtmp,memdot,UVpdot,Wdot)
       UVp_adj(:,:,2) = UVp_adj(:,:,1) + UVpdot*dt
       W_adj(:,:,2) = W_adj(:,:,1) + Wdot*dt
       mem_adj(:,:,:,:,:,:,2) = mem_adj(:,:,:,:,:,:,1) + memdot*dt

       ! Check to see if an elastic change occurs
       if (present(point_data)) then
          call adj_load_misfit_point_rel(nt,it,dt,sl_fin,UVp,point_data,adj_load,load_flag)
       else if (present(rsh_data)) then
          ! UVp here needs changing!!!!!!!!!!
          call adj_load_misfit_rsh(nt,it,dt,sl_fin(:,:,nt-it+1),UVp(:,:,nt-it+1),rsh_data, &
                                   adj_load,load_flag)
       else 
          print *, "No data"
          return
       end if
       
       if (load_flag) then
          call fun_from_coefs(adj_load%load_lm,load_spat)
          call string_cat_int('adjload.',it,file1)
          open(newunit=io,file=trim(file1))
          do iphi = 1,nphi
             do igl = 1,ngl
                write(io,*) iphi,igl,real(load_spat(igl,iphi)),imag(load_spat(igl,iphi))
             end do
             write(io,*)
          end do
          close(io)

          call stat_adj(ice_load,sl_fin(:,:,nt-it+1),adj_load,delta_UVpa(:,:))
          UVp_adj(:,:,2) = UVp_adj(:,:,2) + delta_UVpa(:,:)
       end if

       rem = mod(it-1,it_out)
       if (rem == 0) then
          iout = iout - 1
          UVptmp = UVp(:,:,iout)
          Wtmp = W(:,:,iout)
          memtmp = mem(:,:,:,:,:,:,iout)
       else
          fac = (it_out - rem)/it_out
          UVptmp = UVp(:,:,iout-1) + fac*(UVp(:,:,iout) - UVp(:,:,iout-1))
          Wtmp = W(:,:,iout-1) + fac*(W(:,:,iout) - W(:,:,iout-1))
          memtmp = mem(:,:,:,:,:,:,iout-1) + fac*(mem(:,:,:,:,:,:,iout) - mem(:,:,:,:,:,:,iout-1))
       end if

       call kern_visc(it,nt,dt,UVptmp,Wtmp,UVp_adj(:,:,2),W_adj(:,:,2), &
                      memtmp,mem_adj(:,:,:,:,:,:,2),kern)

       UVp_adj(:,:,1) = UVp_adj(:,:,2)
       W_adj(:,:,1) = W_adj(:,:,2)
       mem_adj(:,:,:,:,:,:,1) = mem_adj(:,:,:,:,:,:,2)
       
    end do


  end subroutine adj_kern_visc_misfit_rel



  !=============================================!
  !      Routines for sensitivity kernels       !
  !=============================================!


  subroutine kern_visc(it,nt,dt,UVp,W,UVp_adj,W_adj,mem,mem_adj,kern)
    
    implicit none

    integer(i4b), intent(in) :: it,nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp,UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: W,W_adj
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: mem,mem_adj
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(inout) :: kern

    integer(i4b) :: l,ispec1,lmcur,m,ilayer,ispec,inode,ispec11,i,igl,iphi,ua,va,wa
    real(dp) :: rl,rmu,rsii,rr,zeta2,zeta2m2,zeta
    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV,dW,dUa,dVa,dWa
    complex(dpc), dimension(5,(lmax+1)**2) :: dcan,dcana
    complex(dpc), dimension(5,ngl,nphi) :: dspat,dspata
    complex(dpc), dimension(5,ngl,nphi) :: memspat,memspata
    ! ONLY BOTHERING WITH ITAU = 1 HERE

    call calculate_dUVW_full(UVp,W,dU,dV,dW)
    call calculate_dUVW_full(UVp_adj,W_adj,dUa,dVa,dWa)
    
    do ilayer = 1,nsect
       
       ! No terms in fluid or elastic regions
       if (sect_ind(3,ilayer) /= 2) cycle

       ! begin loop over the spectral elements
       !$OMP PARALLEL DO PRIVATE(ua,va,wa,ispec1,lmcur,zeta2,zeta,zeta2m2), &
       !$OMP& PRIVATE(rl,rr,rmu,rsii,l,m,i), &
       !$OMP& FIRSTPRIVATE(dcan,dcana,dspat,dspata,memspat,memspata), &
       !$OMP& SHARED(lmax,it,nt,dt,UVp,UVp_adj,W,W_adj,mem,mem_adj,kern,dU,dV,dW,dUa,dVa,dWa)
       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
                    
          do inode = 1,ngll

             dcan = 0.0_dp
             dcana = 0.0_dp

             dspat = 0.0_dp
             dspata = 0.0_dp
             
             memspat = 0.0_dp
             memspata = 0.0_dp

             rmu = mui_node(1,inode,ispec)   
             rr = r_node(inode,ispec)
             if ((ispec == 1) .and. (inode == 1)) then
                rr = r_node(inode+1,ispec)
             end if

             lmcur = 1
             do l = 1,lmax

                call get_depth(l,rl,ispec1)
                if (ispec1 > ispec) then
                   exit
                end if

                ua = gvn_ssg(inode,ispec,ispec1,2)
                va = ua + 1
                wa = gvn_tor(inode,ispec,ispec1)

                zeta2 = l*(l+1)
                zeta2m2 = zeta2 - 2.0_dp
                zeta = sqrt(zeta2)
             
                do m = -l,l
                   lmcur = lmcur + 1

                   dcan(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp(lmcur,va) - ii*W(lmcur,wa))
                   dcan(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                   + UVp(lmcur,ua) - ii*rr*dW(lmcur,inode,ispec) &
                                   + ii*W(lmcur,wa))
                   dcan(3,lmcur) = (rr*dU(lmcur,inode,ispec) - UVp(lmcur,ua) &
                                   + 0.5_dp*zeta2*UVp(lmcur,va)) &
                                   /(3.0_dp*rr)
                   dcan(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                   + UVp(lmcur,ua) + ii*rr*dW(lmcur,inode,ispec) &
                                   - ii*W(lmcur,wa))
                   dcan(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp(lmcur,va) + ii*W(lmcur,wa))

                   dcana(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp_adj(lmcur,va) - ii*W_adj(lmcur,wa))
                   dcana(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dVa(lmcur,inode,ispec) - UVp_adj(lmcur,va) &
                                   + UVp_adj(lmcur,ua) - ii*rr*dWa(lmcur,inode,ispec) &
                                   + ii*W_adj(lmcur,wa))
                   dcana(3,lmcur) = (rr*dUa(lmcur,inode,ispec) - UVp_adj(lmcur,ua) &
                                   + 0.5_dp*zeta2*UVp_adj(lmcur,va)) &
                                   /(3.0_dp*rr)
                   dcana(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dVa(lmcur,inode,ispec) - UVp_adj(lmcur,va) &
                                   + UVp_adj(lmcur,ua) + ii*rr*dWa(lmcur,inode,ispec) &
                                   - ii*W_adj(lmcur,wa))
                   dcana(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp_adj(lmcur,va) + ii*W_adj(lmcur,wa))


                end do

             end do

             ! Calculate spatial fields
             call fun_from_coefs_ten_tl(dcan,dspat)
             call fun_from_coefs_ten_tl(dcana,dspata)

             memspat = mem(:,:,:,1,inode,ispec)
             memspata = mem_adj(:,:,:,1,inode,ispec)

             do iphi = 1,nphi
                
                do igl = 1,ngl        

                   rsii = si_spat_node(igl,iphi,1,inode,ispec)

                   do i = 1,5

                      if ((it == 1) .or. (it == nt)) then

                         if ((i == 1) .or. (i == 5)) then

                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                                    + rmu*rsii*dt &
                                    *(dspat(i,igl,iphi)-memspat(i,igl,iphi)) &
                                    *(dspata(6-i,igl,iphi)-memspata(6-i,igl,iphi))

                         else if (i == 3) then

                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                          + 6.0_dp*rmu*rsii*dt &
                            *(dspat(i,igl,iphi)- memspat(i,igl,iphi)) &
                            *(dspata(i,igl,iphi)- memspata(i,igl,iphi))

                         else if ((i == 2) .or. (i == 4)) then
                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                               - 2.0_dp*rmu*rsii*dt &
                                 *(dspat(i,igl,iphi)- memspat(i,igl,iphi)) &
                                 *(dspata(6-i,igl,iphi)- memspata(6-i,igl,iphi))

                         end if

                      else

                         if ((i == 1) .or. (i == 5)) then

                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                                    + 2.0_dp*rmu*rsii*dt &
                                    *(dspat(i,igl,iphi)-memspat(i,igl,iphi)) &
                                    *(dspata(6-i,igl,iphi)-memspata(6-i,igl,iphi))

                         else if (i == 3) then

                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                          + 12.0_dp*rmu*rsii*dt &
                            *(dspat(i,igl,iphi)- memspat(i,igl,iphi)) &
                            *(dspata(i,igl,iphi)- memspata(i,igl,iphi))

                         else if ((i == 2) .or. (i == 4)) then
                            kern(igl,iphi,inode,ispec) = kern(igl,iphi,inode,ispec) &
                               - 4.0_dp*rmu*rsii*dt &
                                 *(dspat(i,igl,iphi)- memspat(i,igl,iphi)) &
                                 *(dspata(6-i,igl,iphi)- memspata(6-i,igl,iphi))

                         end if


                      end if


                   end do
                      
                end do
                   
             end do

          end do

       end do
       !$OMP END PARALLEL DO

    end do

    return

  end subroutine kern_visc


  subroutine adj_kern_ice_full(nt,dt,nout,ice_data,sl_fin,adj_load,ice_kern)

    implicit none

    integer(i4b), intent(in) :: nt,nout
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_fin
    type(meas), intent(in) :: adj_load
    complex(dpc), dimension(ngl,nphi,nout), intent(out) :: ice_kern

    integer(i4b) :: it,it1
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor,2) :: W_adj
    complex(dpc), dimension(ngl,nphi,2) :: Q_adj
    complex(dpc), dimension(2) :: K_adj
    complex(dpc), dimension(ngl,nphi) :: ice_load,ice_adj,sl_adj
    complex(dpc), dimension((lmax+1)**2) :: Us,ps
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec,2) :: mem

    ice_kern = 0.0_dp

    ! Get ice coverage
    call time_point_interp(nt,dt,ice_data,ice_load)

    call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1),Q_adj(:,:,1),K_adj(1))


    ps = 0.0_dp
    Us = 0.0_dp
    call surface_fields(UVp_adj(:,:,1),ps,Us)

    call adj_sl_ice(ice_load,sl_fin(:,:,nt),adj_load,Us,ps,Q_adj(:,:,1),K_adj(1),sl_adj,ice_adj)

    call kern_ice_it(ice_load,sl_fin(:,:,nt),Us,ps,sl_adj,ice_adj,ice_kern(:,:,nout))

    do it = 2,nt

       print *, "adjoint time step = ",it," of ",nt

       ! Interpolate to find rate of change of load
       ! We want load at previous time step and time reversed
       it1 = time_point(nt-it+2,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+2,dt,ice_data(:,:,it1:it1+1),ice_load)
       end if

       call update_full_rk_adj(dt,ice_load,sl_fin(:,:,nt-it+2),mem(:,:,:,:,:,:,1), &
                               UVp_adj(:,:,1),W_adj(:,:,1),mem(:,:,:,:,:,:,2),UVp_adj(:,:,2), &
                               W_adj(:,:,2),Q_adj(:,:,1),K_adj(1),Q_adj(:,:,2),K_adj(2))

       ps = 0.0_dp
       Us = 0.0_dp
       call surface_fields(UVp_adj(:,:,2),ps,Us)

       ! Calculate load at CURRENT reversed time
       it1 = time_point(nt-it+1,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+1,dt,ice_data(:,:,it1:it1+1),ice_load)
       end if

       call adj_sl_ice(ice_load,sl_fin(:,:,nt-it+1),adj_load,Us,ps,Q_adj(:,:,2),K_adj(2), &
                       sl_adj,ice_adj)

       ! Calculate kernel
       if (mod(it-1,(nt-1)/(nout-1)) == 0) then
          call kern_ice_it(ice_load,sl_fin(:,:,nt-it+1),Us,ps, &
                           sl_adj,ice_adj,ice_kern(:,:,nout - (it-1)*(nout-1)/(nt-1)))
       end if

       UVp_adj(:,:,1) = UVp_adj(:,:,2)
       W_adj(:,:,1) = W_adj(:,:,2)
       mem(:,:,:,:,:,:,1) = mem(:,:,:,:,:,:,2)
       Q_adj(:,:,1) = Q_adj(:,:,2)
       K_adj(1) = K_adj(2)

    end do


  end subroutine adj_kern_ice_full


  subroutine kern_ice_it(ice_load,sl_for,Us_adj,ps_adj,sl_adj,ice_adj,kern)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_load,sl_for
    complex(dpc), dimension((lmax+1)**2), intent(in) :: Us_adj,ps_adj
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_adj,ice_adj
    complex(dpc), dimension(ngl,nphi), intent(out) :: kern

    real(dp) :: area,gsurf,rsurf
    complex(dpc), dimension(ngl,nphi) :: ocean_func,Uadj,padj,CUp
    complex(dpc), dimension((lmax+1)**2) :: oflm,CUp_lm,sla_lm

    integer(i4b) :: io,iphi,igl
    real(dp) :: lat,long

    gsurf = grav_node(ngll,nspec)
    rsurf = r_node(ngll,nspec)

    ! Calculate ocean functions
    ocean_func = ocean_function(ice_load,sl_for)
    call coefs_from_fun(ocean_func,oflm)

    ! Calculate area
    area = sqrt(fourpi_d)*(rsurf**2)*oflm(1)

    ! Calculate spatial Us and ps
    call fun_from_coefs(Us_adj,Uadj)
    call fun_from_coefs(ps_adj,padj)

    ! Calculate CUp
    CUp = ocean_func*(gsurf*Uadj + padj)

    ! Calculate CUp_lm
    call coefs_from_fun(CUp,CUp_lm)

    ! Calculate sladj_lm
    call coefs_from_fun(sl_adj,sla_lm)

    kern = rice*(1.0_dp - ocean_func)*(gsurf*Uadj + padj) - rice*gsurf*ice_adj &
            - rice*(1.0_dp - ocean_func)*sqrt(fourpi_d)*(rsurf**2) &
              *(gsurf*sla_lm(1) + CUp_lm(1))/area
    

  end subroutine kern_ice_it



  ! INITIAL SEA LEVEL KERNEL

  subroutine surface_misfit(func1,func2,misfit)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: func1,func2
    complex(dpc), intent(out) :: misfit

    integer(i4b) :: io,igl,iphi
    complex(dpc), dimension(ngl,nphi) :: spatial_misfit

    spatial_misfit = 0.5_dp*(func1 - func2)**2.0_dp

    open(newunit=io,file='sldiff.out')
    do iphi = 1,nphi
       do igl=1,ngl
          write(io,*) iphi,igl,real(spatial_misfit(igl,iphi))
       end do
       write(io,*)
    end do
    close(io)

    call zero_coef_from_fun(spatial_misfit,misfit)

    misfit = misfit*sqrt(fourpi_d)*(r_node(ngll,nspec))**2.0_dp

  end subroutine surface_misfit


  subroutine real_sph_harm_misfit(func1,func2,misfit)
    ! func1 is real spherical harmonic coefficients
    ! func2 is complex spherical harmonic coefficients

    implicit none

    real(dp), dimension((lmax+1)**2), intent(in) :: func1
    complex(dpc), dimension((lmax+1)**2), intent(in) :: func2
    complex(dpc), intent(out) :: misfit

    integer(i4b) :: l,m,lmr,lmc

    misfit = 0.0_dp

    ! Loop over REAL l,m
    do l = 2,lmax
       do m = -l,l
          lmr = l**2 + l + m + 1
          if (m <= 0) then
             lmc = l**2 + l - m + 1
             misfit = misfit + 0.5_dp*(func1(lmr) - sqrt(2.0_dp)*real(func2(lmc)))**2
          else
             lmc = lmr
             misfit = misfit + 0.5_dp*(func1(lmr) + sqrt(2.0_dp)*imag(func2(lmc)))**2
          end if
       end do
    end do

  end subroutine real_sph_harm_misfit


  subroutine visc_inversion(t2,ice_data,sl_init,comb_data,file1)

    implicit none

    real(dp), intent(in) :: t2
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_init
    type(comb_meas), intent(in) :: comb_data
    character(len=256), intent(in) :: file1

    logical(lgt) :: negvisc
    integer(i4b) :: nloc,i,io,io2,ilayer,inode,ispec,igl,iphi,io3,itmax,ncoef,io4,nt,nout
    real(dp) :: misfit,rr,long,lat,tmin,dt
    character(len=256) :: visc_suff
    complex(dpc) :: alpha,beta
    complex(dpc), dimension(:,:,:), allocatable :: sl_out
    complex(dpc), dimension(:,:,:), allocatable :: UVp
    complex(dpc), dimension(:,:,:), allocatable :: W
    complex(dpc), dimension(:,:,:,:,:,:,:), allocatable :: mem
    complex(dpc), dimension(ngl,nphi,ngll,nspec) :: kern,pgrad,kerno,kernsm
    type(point_meas), dimension(:), allocatable :: point_data
    type(rsh_meas), dimension(:), allocatable :: rsh_data
    !type(rsh_meas) :: rsh_data

    call get_nt_nout(t2,dt,nt,nout)

    allocate(sl_out(ngl,nphi,nt))
    allocate(UVp((lmax+1)**2,nglob_ssg,nout))
    allocate(W((lmax+1)**2,nglob_tor,nout))
    allocate(mem(5,ngl,nphi,ntau,ngll,nspec,nout))

    if (comb_data%pres_point) then
       nloc = comb_data%n_point
       allocate(point_data(nloc))
       point_data = comb_data%point
    end if

    if (comb_data%pres_rsh) then
       ncoef = comb_data%n_rsh
       allocate(rsh_data(ncoef))
       rsh_data = comb_data%rsh
    end if

    ! Open misfit file
    open(newunit=io,file = trim(file1) // trim('misfit'))

    ! Initial forward calculation
    call forward_full_sl(nt,dt,nt,ice_data,sl_init,sl_out,UVp,W,mem,nout)

    if (comb_data%pres_point) then

       ! Calculate misfit - assuming all measurements are sea level
       call misfit_calc_point_rel(nt,dt,sl_out,UVp,point_data,misfit)

       write(io,*) 0, misfit*(r_norm**2)
       print *, misfit*(r_norm**2)

       ! Adjoint and kernel calculation
       call adj_kern_visc_misfit_rel(nt,dt,nout,ice_data,sl_out,UVp,W,mem,kern, &
                                     point_data=point_data)

    else if (comb_data%pres_rsh) then
       
       call misfit_calc_rsh(nt,dt,sl_out,UVp,rsh_data,misfit)

       write(io,*) 0, misfit*(pot_norm**2)

       ! Adjoint and kernel calculation
       call adj_kern_visc_misfit(nt,dt,nout,ice_data,sl_out,UVp,W,mem,kern,rsh_data=rsh_data)

    else
       print *, "No data supplied"
       return
    end if

    call smooth_kern(kern,kernsm)

    pgrad = -kernsm

    open(newunit=io3,file=trim(file1) // trim('tmin'))

    itmax = 10

    do i = 1,itmax

       ! Write out kernel
       call string_cat_int('kern.',i,visc_suff)
       open(newunit=io2,file = trim(file1) // trim(visc_suff))
       do ilayer = 1,nsect

          if (sect_ind(3,ilayer) /= 2) cycle

          do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
             do inode = 1,ngll
                rr = r_node(inode,ispec)*r_norm
                do iphi = 1,nphi
                   long = aphi(iphi)*rad2deg - 180.0_dp
                   do igl = 1,ngl
                      lat = (pio2_d - tgl(igl))*rad2deg
                      write(io2,*) rr,long,lat,kernsm(igl,iphi,inode,ispec), &
                                   pgrad(igl,iphi,inode,ispec)
                   end do
                   write(io2,*) rr, "ABC"
                   write(io2,*) "ABC", long
                end do

             end do
          end do

       end do
       close(io2)

       ! Calculate alpha
       if (comb_data%pres_point) then
          call alpha_calc(misfit,kern,pgrad,t2,ice_data,sl_init, &
                          point_data=point_data,alpha=alpha)
       else if (comb_data%pres_rsh) then
          call alpha_calc(misfit,kern,pgrad,t2,ice_data,sl_init, &
                          rsh_data=rsh_data,alpha=alpha)
       end if

       negvisc = .true.
       do while (negvisc)
          negvisc = .false.
          do ilayer = 1,nsect
             ! No terms in fluid or elastic regions
             if (sect_ind(3,ilayer) /= 2) cycle
             
             do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
                do inode = 1,ngll
                   do iphi = 1,nphi
                      do igl = 1,ngl
                         if (real(alpha)*real(pgrad(igl,iphi,inode,ispec)) >= 1.0_dp) then
                            alpha = alpha/2.0_dp
                            negvisc = .true.
                            exit
                         end if
                      end do
                      if (negvisc) exit
                   end do
                   if (negvisc) exit
                end do
                if (negvisc) exit
             end do
             if (negvisc) exit
          end do

       end do

       ! Update s by actually updating viscosity - changed my mind again
       !si_spat_node(:,:,1,:,:) = si_spat_node(:,:,1,:,:)/(1.0_dp + real(alpha)*real(pgrad))

       si_spat_node(:,:,1,:,:) = si_spat_node(:,:,1,:,:)*(1.0_dp - real(alpha)*real(pgrad))


       smax = 0.0_dp
       ! Write out viscosity
       call string_cat_int('visc.',i,visc_suff)
       open(newunit=io2,file = trim(file1) // trim(visc_suff))
       do ilayer = 1,nsect

          if (sect_ind(3,ilayer) /= 2) cycle

          do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
             do inode = 1,ngll
                rr = r_node(inode,ispec)*r_norm
                do iphi = 1,nphi
                   long = aphi(iphi)*rad2deg - 180.0_dp
                   do igl = 1,ngl
                      lat = (pio2_d - tgl(igl))*rad2deg
                      write(io2,*) rr,long,lat, &
                                   mu_node(inode,ispec)*visco_norm/real(si_spat_node(igl,iphi,1,inode,ispec))
                      if (real(si_spat_node(igl,iphi,1,inode,ispec)) > smax) then
                         smax = real(si_spat_node(igl,iphi,1,inode,ispec))
                      end if
                   end do
                   write(io2,*) rr, "ABC"
                   write(io2,*) "ABC", long
                end do

             end do
          end do

       end do
       close(io2)

       write(io3,*) i,(1.0_dp/smax)*t_norm*sec2yr

       call get_nt_nout(t2,dt,nt,nout)

       deallocate(sl_out,UVp,W,mem)
       allocate(sl_out(ngl,nphi,nt))
       allocate(UVp((lmax+1)**2,nglob_ssg,nout))
       allocate(W((lmax+1)**2,nglob_tor,nout))
       allocate(mem(5,ngl,nphi,ntau,ngll,nspec,nout))

       ! Forward calculation
       call forward_full_sl(nt,dt,nt,ice_data,sl_init,sl_out,UVp,W,mem,nout_opt=nout)

       ! Calculate misfit
       if (comb_data%pres_point) then
          call misfit_calc_point_rel(nt,dt,sl_out,UVp,point_data,misfit)
          write(io,*) i, misfit*(r_norm**2)
       else if (comb_data%pres_rsh) then       
          call misfit_calc_rsh(nt,dt,sl_out,UVp,rsh_data,misfit)
          write(io,*) i, misfit*(pot_norm**2)
       end if 
       
       if (i == itmax) exit

       ! Save old kernel
       kerno = kernsm

       ! Adjoint and kernel calculation
       if (comb_data%pres_point) then
          call adj_kern_visc_misfit_rel(nt,dt,nout,ice_data,sl_out,UVp,W,mem,kern, &
                                        point_data=point_data)
       else if (comb_data%pres_rsh) then
          call adj_kern_visc_misfit(nt,dt,nout,ice_data,sl_out,UVp,W,mem,kern, &
                                    rsh_data=rsh_data)
       end if

       call smooth_kern(kern,kernsm)

       ! Calculate beta
       call beta_calc(kernsm,kerno,beta)

       ! Update p
       pgrad = -kernsm + real(beta)*real(pgrad)

    end do

    close(io)
    close(io3)
    close(io4)

  end subroutine visc_inversion

  
  subroutine misfit_calc_point(nt,dt,sl_out,UVp,point_data,misfit)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2, nglob_ssg,nt), intent(in) :: UVp
    type(point_meas), dimension(:), intent(in) :: point_data
    real(dp), intent(out) :: misfit

    integer(i4b) :: iloc,itime,igl1,iphi1,it1,it2,iphi2,igl2
    real(dp) :: lat,long,xx,yy,time,tt
    complex(dpc) :: f111,f112,f121,f122,f211,f221,f212,f222,f_interp

    ! HAVEN'T ACTUALLY INCLUDED U AND P

    misfit = 0.0_dp

    do iloc = 1,size(point_data)
       lat = point_data(iloc)%lat
       long = point_data(iloc)%long

       call get_igl_iphi(lat,long,igl1,iphi1,xx,yy)
       iphi2 = iphi1 + 1
       igl2 = igl1 + 1

       do itime = 1,point_data(iloc)%nt
          time = point_data(iloc)%times(itime)

          call get_it(nt,dt,time,it1,tt)
          it2 = it1 + 1

          if (point_data(iloc)%type == 1) then
          
             f111 = sl_out(igl1,iphi1,it1)
             f112 = sl_out(igl2,iphi1,it1)
             f121 = sl_out(igl1,iphi2,it1)
             f122 = sl_out(igl2,iphi2,it1)
             f211 = sl_out(igl1,iphi1,it2)
             f212 = sl_out(igl2,iphi1,it2)
             f221 = sl_out(igl1,iphi2,it2)
             f222 = sl_out(igl2,iphi2,it2)

          end if

          f_interp = f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                      + (f122 - f112 - f121 + f111)*xx*yy &
                      + tt*(f211+ (f221 - f211)*xx + (f212 - f211)*yy &
                      + (f222 - f212 - f221 + f211)*xx*yy &
                      - (f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                      + (f122 - f112 - f121 + f111)*xx*yy))

          !print *, iloc,itime,real(f111)*r_norm,real(f211)*r_norm, &
          !         real(f_interp)*r_norm,tt,real(point_data(iloc)%values(itime)*r_norm)


          misfit = misfit + 0.5_dp*(f_interp - point_data(iloc)%values(itime))**2

       end do
    end do
       


  end subroutine misfit_calc_point


  subroutine misfit_calc_point_rel(nt,dt,sl_out,UVp,point_data,misfit)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2, nglob_ssg,nt), intent(in) :: UVp
    type(point_meas), dimension(:), intent(in) :: point_data
    real(dp), intent(out) :: misfit

    integer(i4b) :: iloc,itime,igl1,iphi1,it1,it2,iphi2,igl2,io
    real(dp) :: lat,long,xx,yy,time,tt
    complex(dpc) :: f111,f112,f121,f122,f211,f221,f212,f222,f_interp,fp_interp

    ! HAVEN'T ACTUALLY INCLUDED U AND P

    ! DATA IS RELATIVE SEA LEVEL
    ! Apart from last one, which is just present day sea level

    misfit = 0.0_dp

    !open(newunit=io,file='rsldiff.out')

    do iloc = 1,size(point_data)
       lat = point_data(iloc)%lat
       long = point_data(iloc)%long

       call get_igl_iphi(lat,long,igl1,iphi1,xx,yy)
       iphi2 = iphi1 + 1
       igl2 = igl1 + 1

       if (point_data(iloc)%type == 1) then

          f111 = sl_out(igl1,iphi1,nt)
          f112 = sl_out(igl2,iphi1,nt)
          f121 = sl_out(igl1,iphi2,nt)
          f122 = sl_out(igl2,iphi2,nt)

       end if

       fp_interp = f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                   + (f122 - f112 - f121 + f111)*xx*yy

       do itime = 1,point_data(iloc)%nt - 1
          time = point_data(iloc)%times(itime)

          call get_it(nt,dt,time,it1,tt)
          it2 = it1 + 1

          if (point_data(iloc)%type == 1) then
          
             f111 = sl_out(igl1,iphi1,it1)
             f112 = sl_out(igl2,iphi1,it1)
             f121 = sl_out(igl1,iphi2,it1)
             f122 = sl_out(igl2,iphi2,it1)
             f211 = sl_out(igl1,iphi1,it2)
             f212 = sl_out(igl2,iphi1,it2)
             f221 = sl_out(igl1,iphi2,it2)
             f222 = sl_out(igl2,iphi2,it2)

          end if

          f_interp = f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                      + (f122 - f112 - f121 + f111)*xx*yy &
                      + tt*(f211+ (f221 - f211)*xx + (f212 - f211)*yy &
                      + (f222 - f212 - f221 + f211)*xx*yy &
                      - (f111 + (f121 - f111)*xx + (f112 - f111)*yy &
                      + (f122 - f112 - f121 + f111)*xx*yy))

          misfit = misfit + 0.5_dp*(f_interp - fp_interp - point_data(iloc)%values(itime))**2

       end do

       ! Add present day misfit
       itime = point_data(iloc)%nt
       misfit = misfit + 0.5_dp*(fp_interp - point_data(iloc)%values(itime))**2

    end do


  end subroutine misfit_calc_point_rel


  subroutine misfit_calc_rsh(nt,dt,sl_out,UVp,rsh_data,misfit)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2, nglob_ssg,nt), intent(in) :: UVp
    type(rsh_meas), dimension(:), intent(in) :: rsh_data
    real(dp), intent(out) :: misfit

    integer(i4b) :: itime,it1,it2,l,m,lmr,lmc,ic,io
    real(dp) :: time,tt,f_interp,f1,f2

    misfit = 0.0_dp

    open(newunit=io,file='grace.diff')
    do ic = 1,size(rsh_data)

       do itime = 1,rsh_data(ic)%nt
          time = rsh_data(ic)%times(itime)

          call get_it(nt,dt,time,it1,tt)
          it2 = it1 + 1

          ! Loop over REAL l,m
          do l = 2,lmax
             do m = -l,l
                lmr = l**2 + l + m + 1
                if (m < 0) then
                   lmc = l**2 + l - m + 1
                   if (rsh_data(ic)%type == 3) then
                      f1 = sqrt(2.0_dp)*real(UVp(lmc,pa_surf(l),it1))
                      f2 = sqrt(2.0_dp)*real(UVp(lmc,pa_surf(l),it2))
                   end if
                else if (m == 0) then
                   lmc = lmr
                   if (rsh_data(ic)%type == 3) then
                      f1 = real(UVp(lmc,pa_surf(l),it1))
                      f2 = real(UVp(lmc,pa_surf(l),it2))
                   end if
                else
                   lmc = lmr
                   if (rsh_data(ic)%type == 3) then
                      f1 = -sqrt(2.0_dp)*imag(UVp(lmc,pa_surf(l),it1))
                      f2 = -sqrt(2.0_dp)*imag(UVp(lmc,pa_surf(l),it2))
                   end if                   
                end if
                f_interp = f1 + tt*(f2 - f1)
                misfit = misfit + 0.5_dp*(rsh_data(ic)%values(lmr,itime) - f_interp)**2
                write(io,*) l,m,rsh_data(ic)%values(lmr,itime),f_interp
             end do
          end do
       end do
    end do
    close(io)


  end subroutine misfit_calc_rsh


  subroutine adj_load_misfit_point(nt,it,dt,sl_out,UVp,point_data,adj_load,load_flag)

    implicit none

    integer(i4b), intent(in) :: nt,it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    type(point_meas), dimension(:), intent(in) :: point_data
    type(meas), intent(out) :: adj_load
    logical(lgt), intent(out) :: load_flag

    integer(i4b) :: iloc,itime,igl1,iphi1,igl2,iphi2,lm,io
    real(dp) :: tbp,lat,long,xx,yy,f11,f12,f21,f22,f_interp
    complex(dpc), dimension((lmax+1)**2) :: delta_al

    tbp = 21000.0_dp*yr2sec/t_norm - (nt-it)*dt

    allocate(adj_load%load_lm((lmax+1)**2))
    adj_load%load_lm = 0.0_dp
    load_flag = .false.

    do iloc = 1, size(point_data)
       lat = point_data(iloc)%lat
       long = point_data(iloc)%long

       ! Assume all the same type
       adj_load%type = point_data(iloc)%type

       do itime = 1, point_data(iloc)%nt
          if (abs(point_data(iloc)%times(itime) - tbp) < 1.0_dp*yr2sec/t_norm) then
             call point_load(lat,long,0,delta_al)

             ! Weight with difference between calculated and data
             call get_igl_iphi(lat,long,igl1,iphi1,xx,yy)
             igl2 = igl1 + 1
             iphi2 = iphi1 + 1

             if (point_data(iloc)%type == 1) then

                f11 = sl_out(igl1,iphi1)
                f12 = sl_out(igl2,iphi1)
                f21 = sl_out(igl1,iphi2)
                f22 = sl_out(igl2,iphi2)

             else
                print *, "NOT WRITTEN CORRECT BIT"
                return
             end if
             
             f_interp = f11 + (f21 - f11)*xx + (f12 - f11)*yy &
                         + (f22 - f12 - f21 + f11)*xx*yy

             if (real(f_interp - point_data(iloc)%values(itime)) == 0.0_dp) exit
             delta_al = delta_al*(f_interp - point_data(iloc)%values(itime))

             adj_load%load_lm = adj_load%load_lm + delta_al
             load_flag = .true.

             exit
          end if
          if (point_data(iloc)%times(itime) < tbp) exit
       end do
    end do
             
  end subroutine adj_load_misfit_point


  subroutine adj_load_misfit_point_rel(nt,it,dt,sl_out,UVp,point_data,adj_load,load_flag)

    implicit none

    integer(i4b), intent(in) :: nt,it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nt), intent(in) :: UVp
    type(point_meas), dimension(:), intent(in) :: point_data
    type(meas), intent(out) :: adj_load
    logical(lgt), intent(out) :: load_flag

    integer(i4b) :: iloc,itime,igl1,iphi1,igl2,iphi2,lm,io,ith
    real(dp) :: tbp,lat,long,xx,yy,f11,f12,f21,f22,f_interp,fp_interp
    complex(dpc), dimension((lmax+1)**2) :: delta_al

    allocate(adj_load%load_lm((lmax+1)**2))
    adj_load%load_lm = 0.0_dp
    load_flag = .false.

    tbp = 21000.0_dp*yr2sec/t_norm - (nt-it)*dt

    do iloc = 1, size(point_data)
       lat = point_data(iloc)%lat
       long = point_data(iloc)%long

       call point_load(lat,long,0,delta_al)

       ! Assume all the same type
       adj_load%type = point_data(iloc)%type

       call get_igl_iphi(lat,long,igl1,iphi1,xx,yy)
       igl2 = igl1 + 1
       iphi2 = iphi1 + 1

       if (point_data(iloc)%type == 1) then

          f11 = sl_out(igl1,iphi1,nt)
          f12 = sl_out(igl2,iphi1,nt)
          f21 = sl_out(igl1,iphi2,nt)
          f22 = sl_out(igl2,iphi2,nt)

       end if

       fp_interp = f11 + (f21 - f11)*xx + (f12 - f11)*yy &
                   + (f22 - f12 - f21 + f11)*xx*yy

       do itime = 1, point_data(iloc)%nt - 1
          if ((abs(point_data(iloc)%times(itime) - tbp) < 1.0_dp*yr2sec/t_norm) &
               .or. (it == 1)) then
             if (it == 1) then
                ith = floor((point_data(iloc)%times(itime) - 21000.0_dp*yr2sec/t_norm)/dt &
                            + 0.5_dp) + nt
             else
                ith = it
             end if

             if (point_data(iloc)%type == 1) then

                f11 = sl_out(igl1,iphi1,nt-ith+1)
                f12 = sl_out(igl2,iphi1,nt-ith+1)
                f21 = sl_out(igl1,iphi2,nt-ith+1)
                f22 = sl_out(igl2,iphi2,nt-ith+1)

                f_interp = f11 + (f21 - f11)*xx + (f12 - f11)*yy &
                           + (f22 - f12 - f21 + f11)*xx*yy

             else
                print *, "NOT WRITTEN CORRECT BIT"
                return
             end if
                
             if (it == 1) then
                adj_load%load_lm = adj_load%load_lm &
                                - delta_al*(f_interp - fp_interp - point_data(iloc)%values(itime))

             else
                adj_load%load_lm = adj_load%load_lm &
                                + delta_al*(f_interp - fp_interp - point_data(iloc)%values(itime))
             end if
             load_flag = .true.

          end if
          if (point_data(iloc)%times(itime) < tbp) exit
       end do
       ! Add present day if it == 1
       if (it == 1) then
          adj_load%load_lm = adj_load%load_lm &
                             + delta_al*(fp_interp - point_data(iloc)%values(point_data(iloc)%nt))
       end if

    end do

    !close(io)
       
  end subroutine adj_load_misfit_point_rel


  subroutine adj_load_misfit_rsh(nt,it,dt,sl_out,UVp,rsh_data,adj_load,load_flag)

    implicit none

    integer(i4b), intent(in) :: nt,it
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_out
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    type(rsh_meas), dimension(:), intent(in) :: rsh_data
    type(meas), intent(out) :: adj_load
    logical(lgt), intent(out) :: load_flag

    integer(i4b) :: iloc,itime,igl1,iphi1,igl2,iphi2,lmr,io,l,m,lmc,ic
    real(dp) :: tbp,lat,long,xx,yy,f11,f12,f21,f22,f_interp,r_surf,diff
    complex(dpc), dimension((lmax+1)**2) :: delta_al

    tbp = 21000.0_dp*yr2sec/t_norm - (nt-it)*dt

    adj_load%type = rsh_data(1)%type
    allocate(adj_load%load_lm((lmax+1)**2))
    adj_load%load_lm = 0.0_dp
    load_flag = .false.

    r_surf = r_node(ngll,nspec)

    do ic = 1,size(rsh_data)

       do itime = 1, rsh_data(ic)%nt
          if (abs(rsh_data(ic)%times(itime) - tbp) < 1.0_dp*yr2sec/t_norm) then

             do l = 2,lmax
                do m = -l,l
                   lmr = l**2 + l + m + 1
                   if (m < 0) then
                      lmc = l**2 + l - m + 1

                      if (rsh_data(ic)%type == 3) then
                         diff = sqrt(2.0_dp)*real(UVp(lmc,pa_surf(l))) &
                                - rsh_data(ic)%values(lmr,itime)
                      else
                         print *, "NOT WRITTEN CORRECT BIT"
                      end if

                      adj_load%load_lm(lmc) = adj_load%load_lm(lmc) + diff/sqrt(2.0_dp)

                      lmc = l**2 + l + m + 1
                      adj_load%load_lm(lmc) = adj_load%load_lm(lmc) + ((-1)**m)*diff/sqrt(2.0_dp)

                   else if (m == 0) then
                      lmc = lmr

                      if (rsh_data(ic)%type == 3) then
                         diff = real(UVp(lmc,pa_surf(l)) - rsh_data(ic)%values(lmr,itime))
                      else
                         print *, "NOT WRITTEN CORRECT BIT"
                      end if
                   
                      adj_load%load_lm(lmc) = adj_load%load_lm(lmc) + diff
                   
                   else
                      lmc = lmr
                      if (rsh_data(ic)%type == 3) then
                         diff = -sqrt(2.0_dp)*imag(UVp(lmc,pa_surf(l))) &
                                - rsh_data(ic)%values(lmr,itime)
                      else
                         print *, "NOT WRITTEN CORRECT BIT"
                      end if

                      adj_load%load_lm(lmc) = adj_load%load_lm(lmc) - ii*diff/sqrt(2.0_dp)

                      lmc = l**2 + l - m + 1
                      adj_load%load_lm(lmc) = adj_load%load_lm(lmc) &
                                              + ii*((-1)**m)*diff/sqrt(2.0_dp)

                   end if
                end do
             end do
             
             load_flag = .true.

             exit
          end if
          if (rsh_data(ic)%times(itime) < tbp) exit
       end do

    end do
             
  end subroutine adj_load_misfit_rsh


  subroutine alpha_calc(misfit,kern,pgrad,t2,ice_data,sl_init,point_data,rsh_data,alpha)

    implicit none

    real(dp), intent(in) :: misfit,t2
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(in) :: kern,pgrad
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_init
    type(point_meas), dimension(:), intent(in), optional :: point_data
    type(rsh_meas), dimension(:), intent(in), optional :: rsh_data
    complex(dpc), intent(out) :: alpha

    logical(lgt) :: negvisc
    integer(i4b) :: ilayer,inode,ispec,nloc,igl,iphi,nt,nout
    real(dp) :: misfit2,smax,rr,dt
    complex(dpc) :: kp,rint,nu
    complex(dpc), dimension(ngl,nphi) :: rprod
    complex(dpc), dimension(ngl,nphi,ntau,ngll,nspec) :: si_spat_node_save
    complex(dpc), dimension(:,:,:), allocatable :: sl_out
    complex(dpc), dimension(:,:,:), allocatable :: UVp

    kp = 0.0_dp

    ! Calculate <kern,p>
    do ilayer = 1,nsect
       ! No terms in fluid or elastic regions
       if (sect_ind(3,ilayer) /= 2) cycle

       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
          do inode = 1,ngll
             rr = r_node(inode,ispec)
             rprod = 0.0_dp
             rint = 0.0_dp

             rprod(:,:) = real(kern(:,:,inode,ispec))*real(pgrad(:,:,inode,ispec))
       
             call zero_coef_from_fun(rprod,rint)

             kp = kp + wgll(inode)*jac_element(ispec)*rint*sqrt(fourpi_d)*(rr**2.0_dp)

          end do

       end do

    end do

    nu = -2.0_dp*real(misfit)/real(kp)

    ! Update s by updating viscosity
    negvisc = .true.
    do while (negvisc)
       negvisc = .false.
       do ilayer = 1,nsect
          ! No terms in fluid or elastic regions
          if (sect_ind(3,ilayer) /= 2) cycle

          do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
             do inode = 1,ngll
                do iphi = 1,nphi
                   do igl = 1,ngl
                      if (real(nu)*real(pgrad(igl,iphi,inode,ispec)) >= 1.0_dp) then
                         nu = nu/2.0_dp
                         negvisc = .true.
                         exit
                      end if
                   end do
                   if (negvisc) exit
                end do
                if (negvisc) exit
             end do
             if (negvisc) exit
          end do
          if (negvisc) exit
       end do

    end do
    ! Save old viscosity
    si_spat_node_save = si_spat_node

    ! Update s 
    si_spat_node(:,:,1,:,:) = si_spat_node(:,:,1,:,:)*(1.0_dp - real(nu)*real(pgrad(:,:,:,:)))
    struct3d = .true.

    call calc_smax_3d

    call get_nt_nout(t2,dt,nt,nout)

    allocate(sl_out(ngl,nphi,nt))
    allocate(UVp((lmax+1)**2,nglob_ssg,nout))

    ! Initial forward calculation
    call forward_full_sl(nt,dt,nt,ice_data,sl_init,sl_out,UVp,nout_opt=nout)

    if (present(point_data)) then

       nloc = size(point_data)

       ! Calculate misfit
       call misfit_calc_point_rel(nt,dt,sl_out,UVp,point_data,misfit2)

    else if (present(rsh_data)) then

       call misfit_calc_rsh(nt,dt,sl_out,UVp,rsh_data,misfit2)

    end if

    ! Calculate alpha
    alpha = kp*(nu)**2/(2.0_dp*(misfit - misfit2 + kp*nu))

    print *, "kp, nu, misfit, misfit2, alpha = ",kp,nu,misfit,misfit2,alpha
    !print *, "smax alpha = ", (1.0_dp/smax)*t_norm*sec2yr

    ! Return to previous viscosity
    si_spat_node = si_spat_node_save


  end subroutine alpha_calc


  subroutine beta_calc(kern,kerno,beta)

    implicit none

    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(in) :: kern,kerno
    complex(dpc), intent(out) :: beta

    integer(i4b) :: ilayer,ispec,inode
    real(dp) :: rr
    complex(dpc) :: kkp1,kk,rint
    complex(dpc), dimension(ngl,nphi) :: rprod

    kkp1 = 0.0_dp
    kk = 0.0_dp

    ! Calculate <kern,p>
    do ilayer = 1,nsect
       ! No terms in fluid or elastic regions
       if (sect_ind(3,ilayer) /= 2) cycle

       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
          do inode = 1,ngll
             rr = r_node(inode,ispec)
             rprod = 0.0_dp
             rint = 0.0_dp

             ! Denominator
             rprod(:,:) = real(kerno(:,:,inode,ispec))*real(kerno(:,:,inode,ispec))
             call zero_coef_from_fun(rprod,rint)
             kk = kk + wgll(inode)*jac_element(ispec)*sqrt(fourpi_d)*(rr**2)*rint

             ! Numerator
             rprod(:,:) = real(kern(:,:,inode,ispec))*real(kern(:,:,inode,ispec) - kerno(:,:,inode,ispec))
             call zero_coef_from_fun(rprod,rint)
             kkp1 = kkp1 + wgll(inode)*jac_element(ispec)*sqrt(fourpi_d)*(rr**2)*rint

          end do

       end do

    end do

    beta = real(kkp1)/real(kk)


  end subroutine beta_calc


  subroutine get_nt_nout(t2,dt,nt,nout)

    real(dp), intent(in) :: t2
    real(dp), intent(out) :: dt
    integer(i4b), intent(out) :: nt,nout

    real(dp) :: tmin

    tmin = t_norm*sec2yr/smax

    if (tmin > 100.0_dp) then
       dt = 50.0_dp*yr2sec/t_norm
    else if (tmin > 50.0_dp) then
       dt = 25.0_dp*yr2sec/t_norm
    else if (tmin > 40.0_dp) then
       dt = 20.0_dp*yr2sec/t_norm
    else if (tmin > 20.0_dp) then
       dt = 10.0_dp*yr2sec/t_norm
    else if (tmin > 10.0_dp) then
       dt = 5.0_dp*yr2sec/t_norm
    else if (tmin > 5.0_dp) then
       dt = 2.5_dp*yr2sec/t_norm
    end if

    nt = (t2 + (0.5_dp*dt))/dt + 1
    print *, nt

    if (nt <= 421) then
       nout = nt
    else
       nout = (nt - 1)/10 + 1
    end if

  end subroutine get_nt_nout


  subroutine smooth_kern(kern,kernsm)

    implicit none

    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(in) :: kern
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(out) :: kernsm

    complex(dpc), dimension(ngl,nphi) :: kernr
    complex(dpc), dimension((lmax+1)**2) :: kernlm
    integer(i4b) :: ilayer,inode,ispec,l,m,lmcur

    kernsm = 0.0_dp

    do ilayer = 1,nsect
       ! No terms in fluid or elastic regions
       if (sect_ind(3,ilayer) /= 2) cycle
       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
          do inode = 1,ngll
             kernr = kern(:,:,inode,ispec)
             call coefs_from_fun(kernr,kernlm)
             lmcur = 0
             do l = 0,lmax
                do m = -l,l
                   lmcur = lmcur + 1
                   kernlm(lmcur) = kernlm(lmcur)*exp(-twopi_d*(l+1)/(lmax + 0.5_dp))
                end do
             end do
             call fun_from_coefs(kernlm,kernr)
             kernsm(:,:,inode,ispec) = kernr
          end do
       end do
    end do
             

  end subroutine smooth_kern


end module module_kern
