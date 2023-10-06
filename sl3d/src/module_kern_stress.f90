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

  !=============================================!
  !      Routines for sensitivity kernels       !
  !=============================================!

  subroutine dev_stress_can(lmax,UVp,W,dU,dV,dW,d)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: W
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(in) :: dU,dV,dW
    complex(dpc), dimension(5,(lmax+1)**2,ngll,nspec), intent(inout) :: d

    integer(i4b) :: ispec,inode,ua,va,wa,m,ispec1,ilayer,lmcur,l
    real(dp) :: zeta2,zeta,zeta2m2,rr,rl

    ! A(1) = A(-,-)
    ! A(2) = A(-,0)
    ! A(3) = A(-,+) = 0.5*A(0,0)
    ! A(4) = A(0,+)
    ! A(5) = A(+,+)

    d = 0.0_dp

    do ilayer = 1,nsect

       ! No terms in fluid regions
       if (sect_ind(3,ilayer) == 3) cycle

    
       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)

          do inode = 1,ngll
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

                   d(1,lmcur,inode,ispec) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                            *(UVp(lmcur,va) - ii*W(lmcur,wa))
                   d(2,lmcur,inode,ispec) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                            *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                            + UVp(lmcur,ua) - ii*rr*dW(lmcur,inode,ispec) &
                                            + ii*W(lmcur,wa))
                   d(3,lmcur,inode,ispec) = (rr*dU(lmcur,inode,ispec) - UVp(lmcur,ua) &
                                            + 0.5_dp*zeta2*UVp(lmcur,va)) &
                                            /(3.0_dp*rr)
                   d(4,lmcur,inode,ispec) = zeta/(2.0_dp*sqrt(2.0_dp)*rr) &
                                            *(rr*dV(lmcur,inode,ispec) - UVp(lmcur,va) &
                                            + UVp(lmcur,ua) + ii*rr*dW(lmcur,inode,ispec) &
                                            - ii*W(lmcur,wa))
                   d(5,lmcur,inode,ispec) = zeta*sqrt(zeta2m2)/(2.0_dp*rr) &
                                            *(UVp(lmcur,va) + ii*W(lmcur,wa))


                end do

             end do
          
          end do

       end do

    end do
   

  end subroutine dev_stress_can


  subroutine kern_visc(it,nt,dt,lmax,UVp,W,UVp_adj,W_adj,mem,mem_adj,kern)
    
    implicit none

    integer(i4b), intent(in) :: it,nt,lmax
    real(dp), intent(in) :: dt
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp,UVp_adj
    complex(dpc), dimension((lmax+1)**2,nglob_tor), intent(in) :: W,W_adj
    complex(dpc), dimension(5,ngl,nphi,ntau,ngll,nspec), intent(in) :: mem,mem_adj
    complex(dpc), dimension(ngl,nphi,ngll,nspec), intent(inout) :: kern

    integer(i4b) :: l,ispec1,lmcur,m,ilayer,ispec,inode,ispec11,i,igl,iphi,ua,va,wa
    real(dp) :: rl,rmu,rsii,rr,zeta2,zeta2m2,zeta
    complex(dpc), dimension((lmax+1)**2,ngll,nspec) :: dU,dV,dW,dUa,dVa,dWa
    !complex(dpc), dimension(5,(lmax+1)**2,ngll,nspec) :: dcan,dcana
    complex(dpc), dimension(5,(lmax+1)**2) :: dcan,dcana
    complex(dpc), dimension(5,ngl,nphi) :: dspat,dspata
    complex(dpc), dimension(5,ngl,nphi) :: memspat,memspata
    ! ONLY BOTHERING WITH ITAU = 1 HERE

    call calculate_dUVW_full(lmax,UVp,W,dU,dV,dW)
    call calculate_dUVW_full(lmax,UVp_adj,W_adj,dUa,dVa,dWa)

    !call dev_stress_can(lmax,UVp,W,dU,dV,dW,dcan)
    !call dev_stress_can(lmax,UVp_adj,W_adj,dUa,dVa,dWa,dcana)

    
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
             call fun_from_coefs_ten_tl(lmax,dcan,dspat)
             call fun_from_coefs_ten_tl(lmax,dcana,dspata)

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


  subroutine kern_ice_it(lmax,ice_load,sl_for,Us_adj,ps_adj,sl_adj,ice_adj,kern)

    implicit none

    integer(i4b), intent(in) :: lmax
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

  subroutine surface_misfit(lmax,func1,func2,misfit)

    implicit none

    integer(i4b), intent(in) :: lmax
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

end module module_kern
