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
  
  subroutine adj_kern_visc_full(nt,dt,ice_data,sl_fin,adj_load,UVp,mem,kern)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi,nt), intent(in) :: sl_fin
    type(meas), intent(in) :: adj_load
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,nt), intent(in) :: UVp
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec,nt), intent(in) :: mem
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(out) :: kern

    integer(i4b) :: it,it1
    complex(dpc), dimension(ngl,nphi) :: ice_load,sltmp,ice_tmp
    complex(dpc), dimension((lmax+1)**2,nglob_ssg,2) :: UVp_adj
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec,2) :: mem_adj
    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVpdot,UVptmp
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec) :: memdot,memtmp

    UVp_adj = 0.0_dp
    mem_adj = 0.0_dp
    kern = 0.0_dp

    call time_point_interp(nt,dt,ice_data,ice_load)

    call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1))

    call kern_visc(1,nt,dt,UVp(:,:,nt),UVp_adj(:,:,1), &
                   mem(:,:,:,:,:,nt),mem_adj(:,:,:,:,:,1),kern)


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
       memdot = 0.0_dp

       ! Calculate time derivatives
       call time_derivs_full_adj(ice_load,sl_fin(:,:,nt-it+2),mem_adj(:,:,:,:,:,1), &
                                 UVp_adj(:,:,1),memdot,UVpdot)

       UVptmp = UVp_adj(:,:,1) + UVpdot*dt/2.0_dp
       memtmp = mem_adj(:,:,:,:,:,1) + memdot*dt/2.0_dp
       sltmp = 0.5_dp*(sl_fin(:,:,nt-it+2) + sl_fin(:,:,nt-it+1))

       it1 = time_point(nt-it+1,dt)
       if (it1 == nt_ice5g) then
          ice_load = ice_data(:,:,nt_ice5g)
       else
          call load_interp(nt-it+1,dt,ice_data(:,:,it1:it1+1),ice_tmp)
       end if
       ice_tmp = 0.5_dp*(ice_tmp + ice_load)

       call time_derivs_full_adj(ice_tmp,sltmp,memtmp,UVptmp,memdot,UVpdot)
       UVp_adj(:,:,2) = UVp_adj(:,:,1) + UVpdot*dt
       mem_adj(:,:,:,:,:,2) = mem_adj(:,:,:,:,:,1) + memdot*dt


       call kern_visc(it,nt,dt,UVp(:,:,nt-it+1),UVp_adj(:,:,2), &
                      mem(:,:,:,:,:,nt-it+1),mem_adj(:,:,:,:,:,2),kern)

       UVp_adj(:,:,1) = UVp_adj(:,:,2)
       mem_adj(:,:,:,:,:,1) = mem_adj(:,:,:,:,:,2)
       
    end do



  end subroutine adj_kern_visc_full



  subroutine adj_kern_init_sl(ice_fin,sl_fin,load,kern)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: ice_fin,sl_fin
    type(meas), intent(in) :: load
    complex(dpc), dimension(ngl,nphi), intent(out) :: kern

    complex(dpc), dimension((lmax+1)**2,nglob_ssg) :: UVp_adj
    complex(dpc), dimension(ngl,nphi) :: Q_adj
    complex(dpc) :: K_adj
    complex(dpc), dimension((lmax+1)**2) :: Us_adj,ps_adj
    complex(dpc), dimension(ngl,nphi) :: sl_adj,ice_adj

    kern = 0.0_dp

    call stat_adj(ice_fin,sl_fin,load,UVp_adj,Q_adj,K_adj)
    call surface_fields(UVp_adj,Us_adj,ps_adj)
    call adj_sl_ice(ice_fin,sl_fin,load,Us_adj,ps_adj,Q_adj,K_adj,sl_adj,ice_adj) 

    kern = roce*grav_node(ngll,nspec)*sl_adj


  end subroutine adj_kern_init_sl



  !=============================================!
  !      Routines for sensitivity kernels       !
  !=============================================!


  subroutine kern_visc(it,nt,dt,UVp,UVp_adj,mem,mem_adj,kern)
    
    implicit none

    integer(i4b), intent(in) :: it,nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp,UVp_adj
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec), intent(in) :: mem,mem_adj
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(inout) :: kern

    integer(i4b) :: l,ispec1,lmcur,m,ilayer,ispec,inode,ispec11,i,igl,iphi,ua,va
    real(dp) :: rl,rmu,rsii,rr,zeta2,zeta2m2,zeta
    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV,dUa,dVa
    complex(dpc), dimension(5,(lmax+1)**2) :: dcan,dcana,mcan,mcana
    complex(dpc), dimension(5,ngl,nphi) :: dspat,dspata
    complex(dpc), dimension(5,ngl,nphi) :: memspat,memspata
    ! ONLY BOTHERING WITH ITAU = 1 HERE

    call calculate_dUV_full(lmax,UVp,dU,dV)
    call calculate_dUV_full(lmax,UVp_adj,dUa,dVa)
    
    do ilayer = 1,nsect
       
       ! No terms in fluid or elastic regions
       if (sect_ind(3,ilayer) /= 2) cycle

       ! begin loop over the spectral elements
       !$OMP PARALLEL DO PRIVATE(ua,va,ispec1,lmcur,zeta2,zeta,zeta2m2), &
       !$OMP& PRIVATE(rl,rr,rmu,rsii,l,m,i), &
       !$OMP& FIRSTPRIVATE(dcan,dcana,dspat,dspata,mcan,mcana,memspat,memspata), &
       !$OMP& SHARED(lmax,it,nt,dt,UVp,UVp_adj,mem,mem_adj,kern,dU,dV,dUa,dVa)
       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
                    
          do inode = 1,ngll

             dcan = 0.0_dp
             dcana = 0.0_dp

             mcan = 0.0_dp
             mcana = 0.0_dp

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

                zeta2 = l*(l+1)
                zeta2m2 = zeta2 - 2.0_dp
                zeta = sqrt(zeta2)
             
                do m = -l,l
                   lmcur = lmcur + 1

                   dcan(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp(lmcur,va))
                   dcan(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                   + UVp(lmcur,ua))
                   dcan(3,lmcur) = (rr*dU(lmcur,inode,ispec) - UVp(lmcur,ua) &
                                   + 0.5_dp*zeta2*UVp(lmcur,va)) &
                                   /(3.0_dp*rr)
                   dcan(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                   + UVp(lmcur,ua))
                   dcan(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp(lmcur,va))

                   dcana(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp_adj(lmcur,va))
                   dcana(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dVa(lmcur,inode,ispec) - UVp_adj(lmcur,va) &
                                   + UVp_adj(lmcur,ua))
                   dcana(3,lmcur) = (rr*dUa(lmcur,inode,ispec) - UVp_adj(lmcur,ua) &
                                   + 0.5_dp*zeta2*UVp_adj(lmcur,va)) &
                                   /(3.0_dp*rr)
                   dcana(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *(rr*dVa(lmcur,inode,ispec) - UVp_adj(lmcur,va) &
                                   + UVp_adj(lmcur,ua))
                   dcana(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *(UVp_adj(lmcur,va))

                   mcan(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *mem(1,lmcur,1,inode,ispec)
                   mcan(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *mem(4,lmcur,1,inode,ispec)
                   mcan(3,lmcur) = mem(3,lmcur,1,inode,ispec)/(3.0_dp*rr)
                   mcan(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                   *mem(4,lmcur,1,inode,ispec)
                   mcan(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                   *mem(1,lmcur,1,inode,ispec)

                   mcana(1,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                    *mem_adj(1,lmcur,1,inode,ispec)
                   mcana(2,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                    *mem_adj(4,lmcur,1,inode,ispec)
                   mcana(3,lmcur) = mem_adj(3,lmcur,1,inode,ispec)/(3.0_dp*rr)
                   mcana(4,lmcur) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                    *mem_adj(4,lmcur,1,inode,ispec)
                   mcana(5,lmcur) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                    *mem_adj(1,lmcur,1,inode,ispec)


                end do

             end do

             ! Calculate spatial fields
             call fun_from_coefs_ten_tl(lmax,dcan,dspat)
             call fun_from_coefs_ten_tl(lmax,dcana,dspata)
             call fun_from_coefs_ten_tl(lmax,mcan,memspat)
             call fun_from_coefs_ten_tl(lmax,mcana,memspata)

             do iphi = 1,nphi
                
                do igl = 1,ngl        

                   rsii = si_node(1,inode,ispec)

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
    complex(dpc), dimension(ngl,nphi,2) :: Q_adj
    complex(dpc), dimension(2) :: K_adj
    complex(dpc), dimension(ngl,nphi) :: ice_load,ice_adj,sl_adj
    complex(dpc), dimension((lmax+1)**2) :: Us,ps
    complex(dpc), dimension(5,(lmax+1)**2,ntau,ngll,nspec,2) :: mem

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

       call update_full_rk_adj(dt,ice_load,sl_fin(:,:,nt-it+2),mem(:,:,:,:,:,1), &
                               UVp_adj(:,:,1),mem(:,:,:,:,:,2),UVp_adj(:,:,2), &
                               Q_adj(:,:,1),K_adj(1),Q_adj(:,:,2),K_adj(2))

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
       mem(:,:,:,:,:,1) = mem(:,:,:,:,:,2)
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
    call coefs_from_fun(lmax,ocean_func,oflm)

    ! Calculate area
    area = sqrt(fourpi_d)*(rsurf**2)*oflm(1)

    ! Calculate spatial Us and ps
    call fun_from_coefs(lmax,Us_adj,Uadj)
    call fun_from_coefs(lmax,ps_adj,padj)

    ! Calculate CUp
    CUp = ocean_func*(gsurf*Uadj + padj)

    ! Calculate CUp_lm
    call coefs_from_fun(lmax,CUp,CUp_lm)

    ! Calculate sladj_lm
    call coefs_from_fun(lmax,sl_adj,sla_lm)

    !open(newunit = io,file='iktest3.out')
    !do iphi = 1,nphi
    !   long = aphi(iphi)*rad2deg - 180.0_dp
    !   do igl = 1,ngl
    !      lat = (pio2_d - tgl(igl))*rad2deg
    !      write(io,*) long,lat,real(Uadj(igl,iphi)),real(padj(igl,iphi)), real(ice_adj(igl,iphi)), &
    !                  real(ocean_func(igl,iphi)),real(sl_adj(igl,iphi))
    !   end do
    !   write(io,*)
    !end do
    !close(io)

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

    call zero_coef_from_fun(lmax,spatial_misfit,misfit)

    misfit = misfit*sqrt(fourpi_d)*(r_node(ngll,nspec))**2.0_dp

  end subroutine surface_misfit

  subroutine init_sl_inversion(nt,dt,ice_data,sl_init_in,sl_fin_data,file1)

    implicit none

    integer(i4b), intent(in) :: nt
    real(dp), intent(in) :: dt
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_init_in
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_fin_data
    character(len=256), intent(in) :: file1

    integer(i4b) :: i,io,io2,igl,iphi
    real(dp) :: misfit,long,lat,xx,yy
    character(len=256) :: sl_suff
    complex(dpc) :: alpha,beta,misfitc
    type(meas) :: load
    complex(dpc), dimension(ngl,nphi) :: load_spat,sl_init
    complex(dpc), dimension(ngl,nphi,2) :: sl_out
    complex(dpc), dimension(ngl,nphi) :: kern,pgrad,kerno

    sl_init = sl_init_in

    ! Open misfit file
    open(newunit=io,file = trim(file1) // trim('.misfit'))

    ! Initial sea level
    call string_cat_int('.init.',0,sl_suff)
    open(newunit=io2,file = trim(file1) // trim(sl_suff))
    do iphi = 1,nphi
       long = aphi(iphi)*rad2deg - 180.0_dp
       do igl = 1,ngl
          lat = (pio2_d - tgl(igl))*rad2deg
          call get_xy(lat,long,xx,yy)
          write(io2,*) long,lat,xx,yy,real(sl_init(igl,iphi))*r_norm
       end do
       write(io2,*)
    end do
    close(io2)

    ! Initial forward calculation
    call forward_full_sl(nt,dt,2,ice_data,sl_init,sl_out)

    ! Final sea level
    call string_cat_int('.fin.',0,sl_suff)
    open(newunit=io2,file = trim(file1) // trim(sl_suff))
    do iphi = 1,nphi
       long = aphi(iphi)*rad2deg - 180.0_dp
       do igl = 1,ngl
          lat = (pio2_d - tgl(igl))*rad2deg
          call get_xy(lat,long,xx,yy)
          write(io2,*) long,lat,xx,yy,real(sl_out(igl,iphi,2))*r_norm
       end do
       write(io2,*)
    end do
    close(io2)

    ! Calculate misfit
    call surface_misfit(sl_fin_data,sl_out(:,:,2),misfitc)
    misfit = real(misfitc)

    write(io,*) 0,misfit*(r_norm**2)

    ! Set load
    load%type = 1
    allocate(load%load_lm((lmax+1)**2))
    load_spat = sl_out(:,:,2) - sl_fin_data(:,:)
    call coefs_from_fun(lmax,load_spat,load%load_lm)

    ! Adjoint and kernel calculation  
    call adj_kern_init_sl(ice_data(:,:,nt_ice5g),sl_out(:,:,2),load,kern)

    pgrad = -kern

    do i = 1,10

       ! Calculate alpha
       call alpha_calc_sl(misfit,kern,pgrad,nt,dt,ice_data,sl_init,sl_fin_data,alpha)

       ! Update initial sea level
       sl_init = sl_init + real(alpha)*real(pgrad)

       ! Initial sea level
       call string_cat_int('.init.',i,sl_suff)
       open(newunit=io2,file = trim(file1) // trim(sl_suff))
       do iphi = 1,nphi
          long = aphi(iphi)*rad2deg - 180.0_dp
          do igl = 1,ngl
             lat = (pio2_d - tgl(igl))*rad2deg
             call get_xy(lat,long,xx,yy)
             write(io2,*) long,lat,xx,yy,real(sl_init(igl,iphi))*r_norm
          end do
          write(io2,*)
       end do
       close(io2)

       ! Forward calculation
       call forward_full_sl(nt,dt,2,ice_data,sl_init,sl_out)

       ! Final sea level
       call string_cat_int('.fin.',i,sl_suff)
       open(newunit=io2,file = trim(file1) // trim(sl_suff))
       do iphi = 1,nphi
          long = aphi(iphi)*rad2deg - 180.0_dp
          do igl = 1,ngl
             lat = (pio2_d - tgl(igl))*rad2deg
             call get_xy(lat,long,xx,yy)
             write(io2,*) long,lat,xx,yy,real(sl_out(igl,iphi,2))*r_norm
          end do
          write(io2,*)
       end do
       close(io2)

       ! Calculate misfit
       call surface_misfit(sl_fin_data,sl_out(:,:,2),misfitc)
       misfit = real(misfitc)

       write(io,*) i, misfit*(r_norm**2)


       ! Save old kernel
       kerno = kern

       ! New load
       load_spat = sl_out(:,:,2) - sl_fin_data(:,:)
       call coefs_from_fun(lmax,load_spat,load%load_lm)

       ! Adjoint and kernel calculation  
       call adj_kern_init_sl(ice_data(:,:,nt_ice5g),sl_out(:,:,2),load,kern)

       ! Calculate beta
       call beta_calc_sl(kern,kerno,beta)

       ! Update p
       pgrad = -kern + real(beta)*real(pgrad)

    end do
    close(io)

  end subroutine init_sl_inversion


  subroutine alpha_calc_sl(misfit,kern,pgrad,nt,dt,ice_data,sl_init_in,sl_fin_data,alpha)

    implicit none

    integer(i4b) :: nt
    real(dp), intent(in) :: misfit,dt
    complex(dpc), dimension(ngl,nphi), intent(in) :: kern,pgrad
    complex(dpc), dimension(ngl,nphi,nt_ice5g), intent(in) :: ice_data
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_init_in
    complex(dpc), dimension(ngl,nphi), intent(in) :: sl_fin_data
    complex(dpc), intent(out) :: alpha

    integer(i4b) :: io
    real(dp) :: misfit2
    complex(dpc) :: nu,misfitc,kp
    complex(dpc), dimension(ngl,nphi) :: prod
    complex(dpc), dimension((lmax+1)**2) :: prod_lm
    complex(dpc), dimension(ngl,nphi) :: sl_init_tmp
    complex(dpc), dimension(ngl,nphi,2) :: sl_out


    ! Calculate kp
    prod = real(kern)*real(pgrad)
    prod_lm = 0.0_dp
    call coefs_from_fun(lmax,prod,prod_lm)
    kp = prod_lm(1)*sqrt(fourpi_d)*(r_node(ngll,nspec))**2

    nu = -2.0_dp*misfit/kp
    
    sl_init_tmp = sl_init_in + real(nu)*real(pgrad)

    call forward_full_sl(nt,dt,2,ice_data,sl_init_tmp,sl_out)

    call surface_misfit(sl_fin_data,sl_out(:,:,2),misfitc)
    misfit2 = real(misfitc)

    alpha = kp*(nu)**2/(2.0_dp*(misfit - misfit2 + kp*nu))

    open(newunit=io,file='alpha2.out')
    write(io,*) "misfit = ", misfit*(r_norm**2)
    write(io,*) "misfit2 = ", misfit2*(r_norm**2)
    write(io,*) "kp = ", kp
    write(io,*) "nu = ", nu*(r_norm**2)
    write(io,*) "alpha = ", alpha*(r_norm**2)
    close(io)

  end subroutine alpha_calc_sl



  subroutine beta_calc_sl(kern,kerno,beta)

    implicit none

    complex(dpc), dimension(ngl,nphi), intent(in) :: kern,kerno
    complex(dpc), intent(out) :: beta

    complex(dpc) :: kk,kkp1
    complex(dpc), dimension(ngl,nphi) :: prod
    complex(dpc), dimension((lmax+1)**2) :: prod_lm

    ! Numerator
    prod = real(kern)*real((kern - kerno))
    call coefs_from_fun(lmax,prod,prod_lm)
    kkp1 = prod_lm(1)*sqrt(fourpi_d)*(r_node(ngll,nspec)**2)

    ! Denominator
    prod = real(kerno)*real(kerno)
    call coefs_from_fun(lmax,prod,prod_lm)
    kk = prod_lm(1)*sqrt(fourpi_d)*(r_node(ngll,nspec)**2)

    beta = real(kkp1)/real(kk)


  end subroutine beta_calc_sl



end module module_kern
