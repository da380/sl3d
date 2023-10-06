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


  subroutine calculate_dUV_all(lmax,UVp,dU,dV)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: lmax
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


  subroutine calculate_dUV_full(lmax,UVp,dU,dV)

    use nrtype

    implicit none

    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension((lmax+1)**2,nglob_ssg), intent(in) :: UVp
    complex(dpc), dimension((lmax+1)**2,ngll,nspec), intent(out) :: dU,dV

    integer(i4b) :: ispec,inode,jnode,ilayer,ispec11,l,m,ua,va,lmcur,ispec1
    real(dp) :: hp,jac,rl

    dU = 0.0_dp
    dV = 0.0_dp

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
                   else
                      ua = ua + 3
                   end if
                   va = ua + 1

                   do m = -l,l
                      lmcur = l**2 + 1 + l + m

                      dU(lmcur,inode,ispec) = dU(lmcur,inode,ispec) &
                                              + UVp(lmcur,ua)*hp/jac

                      dV(lmcur,inode,ispec) = dV(lmcur,inode,ispec) &
                                              + UVp(lmcur,va)*hp/jac

                   end do

                end do


             end do
          end do
       end do
    end do


  end subroutine calculate_dUV_full



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

  subroutine surface_fields_full(lmax,UVp,ps,Us)

    implicit none

    integer(i4b), intent(in) :: lmax
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


  subroutine surface_field(lmax,UVp,Us)

    implicit none

    integer(i4b), intent(in) :: lmax
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
