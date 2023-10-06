module module_sphexp

  use nrtype
  use SHTOOLS
  implicit none


  !=============================================================!
  !=============================================================!
  !                variables for use with SHTOOLS               !
  !=============================================================!
  integer(i4b), parameter :: maxdeg_SH = 1000
  integer(i4b), parameter :: maxgrid_SH = 2*(2*maxdeg_SH+1)
  integer(i4b), save :: ilat_SH
  integer(i4b), save :: ilon_SH
  integer(i4b), save :: nlat_SH
  integer(i4b), save :: nlon_SH
  integer(i4b), save :: nlat_SH_GL
  integer(i4b), save :: nlon_SH_GL 
  integer(i4b), save :: lmax_SH
  integer(i4b), save :: ncoef_SH
  integer(i4b), save :: ncoef_SH_real
  integer(i4b), save :: ndh_SH
  integer(i4b), dimension(:), allocatable, save :: ilm1_SH
  integer(i4b), dimension(:), allocatable, save :: ilm2_SH
  real(dp), dimension(maxgrid_SH), save :: lat_SH
  real(dp), dimension(maxgrid_SH), save :: lon_SH
  real(dp), dimension(maxgrid_SH), save :: th_SH
  real(dp), dimension(maxgrid_SH), save :: ph_SH
  real(dp), dimension(maxdeg_SH+1), save :: lat_SH_GL
  real(dp), dimension(2*maxdeg_SH+1), save :: lon_SH_GL
  real(dp), dimension(maxdeg_SH+1), save :: th_SH_GL
  real(dp), dimension(2*maxdeg_SH+1), save :: ph_SH_GL
  real(dp), dimension(maxgrid_SH,maxgrid_SH), save :: fgrid_SH
  real(dp), dimension(2, maxdeg_SH+1, maxdeg_SH+1), save :: flm_SH
  real(dp), dimension(maxdeg_SH+1,2*maxdeg_SH+1) :: fgrid_SH_GL
  real(dp), dimension(maxdeg_SH+1) :: zero_SH_GL
  real(dp), dimension(maxdeg_SH+1) :: w_SH_GL


  


  !---------------------------------------!
  !   variables for the spherical mesh    !
  !---------------------------------------!  
  integer(i4b), parameter :: npw_sph = 3
  integer(i4b), parameter :: ngll_sph = 5
  integer(i4b), save :: lmax_sph
  integer(i4b), save :: ncoef_sph
  integer(i4b), save :: nspec_sph
  integer(i4b), save :: nglob_sph
  integer(i4b), save :: nrec_sph
  integer(i4b), save :: nth_sph
  integer(i4b), save :: nph_sph
  real(dp), save :: jac_sph
  real(dp), save :: dth_sph
  real(dp), save :: dph_sph
  real(dp), dimension(:), allocatable, save :: tha_sph
  real(dp), dimension(:), allocatable, save :: pha_sph

  !-------------------------------------!
  !   variables for the reduced arrays  !
  !-------------------------------------!
  integer(i4b), save :: nthd_sph
  integer(i4b), save :: nphd_sph
  real(dp), dimension(:), allocatable, save :: thd_sph
  real(dp), dimension(:), allocatable, save :: phd_sph
  integer(i4b), dimension(:), allocatable, save :: indth_sph
  integer(i4b), dimension(:), allocatable, save :: indph_sph
  integer(i4b), dimension(:), allocatable, save :: rnkth_sph
  integer(i4b), dimension(:), allocatable, save :: rnkph_sph


  !-----------------------------------------------------!
  ! arrays for the spherical harmonics and exponentials !
  !-----------------------------------------------------!
  real(dp), dimension(:), allocatable, save :: sinth_sph
  real(dp), dimension(:,:), allocatable, save :: xlm_sph
  real(dp), dimension(:,:), allocatable, save :: xplm_sph
  real(dp), dimension(:,:), allocatable, save :: xclm_sph
  complex(dpc), dimension(:,:), allocatable, save :: eimp_sph
  complex(dpc), dimension(:), allocatable, save :: ylm_sph
  complex(dpc), dimension(:), allocatable, save :: dylmdt_sph
  complex(dpc), dimension(:), allocatable, save :: dylmdp_sph

  !---------------------------------------!
  ! arrays for the GLL points and weights !
  !---------------------------------------!
  real(dp), dimension(ngll_sph), save :: xigll_sph                        
  real(dp), dimension(ngll_sph), save :: wgll_sph                         
  real(dp), dimension(ngll_sph), save :: htmp_xi_sph
  real(dp), dimension(ngll_sph), save :: hptmp_xi_sph
  real(dp), dimension(ngll_sph), save :: htmp_eta_sph
  real(dp), dimension(ngll_sph), save :: hptmp_eta_sph



contains

  !================================================================!
  !================================================================!
  !                     wrappers to use SHTOOLS                    !
  !================================================================!
  !================================================================!

  
  subroutine set_SHTOOLS(lmax,ifgl)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: lmax
    logical(lgt), intent(in), optional :: ifgl
    integer(i4b) :: ilat,ilon,ilm1,ilm2,l
    real(dp) :: interval,lat,lon,th,ph,dlat,dlon
    
    lmax_SH = lmax
    ncoef_SH = (lmax_SH+1)*(lmax_SH+2)/2    
    ncoef_SH_real = (lmax+1)**2
    ndh_SH = 2*(lmax_SH+1)
    nlat_SH = ndh_SH
    nlon_SH = 2*ndh_SH
    dlat = 180.0_dp/nlat_SH
    dlon = dlat
    
    do ilat = 1,nlat_SH
       lat = 90.0-(ilat-1)*dlat
       th = (90.0_dp-lat)*deg2rad
       lat_SH(ilat) = lat
       th_SH(ilat) = th
    end do

    do ilon = 1,nlon_SH
       lon = (ilon-1)*dlon
       ph = lon*deg2rad
       lon_SH(ilon) = lon
       ph_SH(ilon) = ph
    end do

    allocate(ilm1_SH(0:lmax_SH))
    allocate(ilm2_SH(0:lmax_SH))

    ! set the ilm1 and ilm2 arrays
    ilm2 = 0
    do l = 0,lmax_SH
       ilm1 = ilm2+1
       ilm2 = ilm1+l
       ilm1_SH(l) = ilm1
       ilm2_SH(l) = ilm2
    end do


    ! set things up for the GL integrations if needed
    if(present(ifgl)) then
       if(ifgl) then
          
          ! get the latitude and longitude positions
          call GLQGridCoord(lat_SH_GL(1:lmax_SH+1),lon_SH_GL(1:2*lmax_SH+1), & 
                            lmax_SH,nlat_SH_GL,nlon_SH_GL)

          th_SH_GL = (90.0_dp-lat_SH_GL)*deg2rad
          ph_SH_GL = lon_SH_GL*deg2rad
          
          
          ! precompute the integration zeros and weights
          call  PreCompute(lmax_SH,zero_SH_GL,w_SH_GL,norm = 4)   


       end if
    end if


    return
  end subroutine set_SHTOOLS


  subroutine SHExpandGLQ_wrapper(flm,ifilt,lmax_calc)
    use nrtype
    use SHTOOLS
    implicit none
    complex(dpc), dimension(:), allocatable, intent(inout) :: flm
    real(dp), intent(in), optional :: ifilt
    integer(i4b), intent(in), optional :: lmax_calc
    
    integer(i4b) :: lmax_use,i,l,ilm1,ilm2,m
    real(dp) :: fac,ifilt_use


    if(present(ifilt)) then
       ifilt_use = ifilt
    else
       ifilt_use = 0.0_dp
    end if

    if(present(lmax_calc)) then
       lmax_use = lmax_calc
    else
       lmax_use = lmax_SH
    end if 
    
    
    if(allocated(flm)) then
       if(size(flm,1) /= ncoef_SH) then
          deallocate(flm)
          allocate(flm(ncoef_SH))
       end if
    else
       allocate(flm(ncoef_SH))
    end if

    call SHExpandGLQ(flm_SH,lmax_SH,fgrid_SH_GL(1:nlat_SH_GL, & 
                     1:nlon_SH_GL),w_SH_GL,ZERO = zero_SH_GL, NORM = 4)


    ! transform to complex coefficients
    i = 0
    do l = 0,lmax_SH
       do m = 0,l
          i = i+1
          if(m == 0) then
             flm(i) =  flm_SH(1,l+1,m+1)
          else
             flm(i) = (flm_SH(1,l+1,m+1)-ii*flm_SH(2,l+1,m+1))/sqrt2
          end if
       end do
    end do
    
    
    if(ifilt_use /= 0.0_dp) then
       ilm2 = 0
       do l = 0,lmax_SH
          ilm1 = ilm2+1
          ilm2 = ilm1+l
          fac = ifilt_use*twopi_d*(l+1)/(lmax_SH+0.5_dp)
          fac = exp(-fac)
          flm(ilm1:ilm2) = fac*flm(ilm1:ilm2)
       end do
    end if



    return
  end subroutine SHExpandGLQ_wrapper




  subroutine SHExpandDH_wrapper(flm,ifilt,lmax_calc)
    use nrtype
    use SHTOOLS
    implicit none    
    complex(dpc), dimension(:), allocatable, intent(inout) :: flm
    real(dp), intent(in), optional :: ifilt
    integer(i4b), intent(in), optional :: lmax_calc
    integer(i4b) :: n,lmax2,i,m,l,ilm1,ilm2,lmax_use
    real(dp) :: fac,ifilt_use

    if(present(ifilt)) then
       ifilt_use = ifilt
    else
       ifilt_use = 0.0_dp
    end if

    if(present(lmax_calc)) then
       lmax_use = lmax_calc
    else
       lmax_use = lmax_SH
    end if 

    if(allocated(flm)) then
       if(size(flm,1) /= ncoef_SH) then
          deallocate(flm)
          allocate(flm(ncoef_SH))
       end if
    else
       allocate(flm(ncoef_SH))
    end if

    n = nlat_SH
    call SHExpandDH(fgrid_SH,n,flm_SH,lmax2,lmax_calc=lmax_use,sampling=2,norm = 4)

      
    ! transform to complex coefficients
    i = 0
    do l = 0,lmax_SH
       do m = 0,l
          i = i+1
          if(m == 0) then
             flm(i) =  flm_SH(1,l+1,m+1)
          else
             flm(i) = (flm_SH(1,l+1,m+1)-ii*flm_SH(2,l+1,m+1))/sqrt2
          end if
       end do
    end do
    
    
    if(ifilt_use /= 0.0_dp) then
       ilm2 = 0
       do l = 0,lmax_SH
          ilm1 = ilm2+1
          ilm2 = ilm1+l
          fac = ifilt_use*twopi_d*(l+1)/(lmax_SH+0.5_dp)
          fac = exp(-fac)
          flm(ilm1:ilm2) = fac*flm(ilm1:ilm2)
       end do
    end if


    return
  end subroutine SHExpandDH_wrapper


  subroutine set_flm_SH(flm,lmax_calc)
    use nrtype
    implicit none
    complex(dpc), dimension(:), intent(in) :: flm
    integer(i4b), intent(in), optional :: lmax_calc
    integer(i4b) :: i,l,m,lmax_use
    
    if(present(lmax_calc)) then
       lmax_use = min(lmax_calc,lmax_SH)
    else
       lmax_use = lmax_SH
    end if


    ! initialize the coefficients
    flm_SH = 0.0_dp

    i = 0
    do l = 0,lmax_use
       do m = 0,l
          i = i+1
          if(m == 0) then
             flm_SH(1,l+1,m+1) = flm(i)
          else
             flm_SH(1,l+1,m+1) =  sqrt2*real(flm(i))
             flm_SH(2,l+1,m+1) = -sqrt2*aimag(flm(i))
          end if
       end do
    end do


    

    return
  end subroutine set_flm_SH


  subroutine MakeGridDH_wrapper(fgrid,lmax_calc)
    use nrtype
    use SHTOOLS
    implicit none
    real(dp), dimension(:,:), allocatable, intent(inout) :: fgrid
    integer(i4b), intent(in), optional :: lmax_calc
    integer(i4b) :: lmax_use


    if(present(lmax_calc)) then
       lmax_use = min(lmax_calc,lmax_SH)
    else
       lmax_use = lmax_SH
    end if


    if(allocated(fgrid)) then
       if(size(fgrid,1) /= nlat_SH .and.  & 
          size(fgrid,2) /= nlon_SH) then
          deallocate(fgrid)
          allocate(fgrid(nlat_SH,nlon_SH))
       end if
    else
       allocate(fgrid(nlat_SH,nlon_SH))
    end if
    
    call MakeGridDH(fgrid,ndh_SH,flm_SH,lmax_SH,lmax_calc=lmax_use,sampling=2,norm=4)


    
    return
  end subroutine MakeGridDH_wrapper


  subroutine real_to_complex_SH(lmax,flm_r,flm_c)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: lmax
    real(dp), dimension(:), intent(in) :: flm_r
    complex(dpc), dimension(:), intent(out) :: flm_c

    integer(i4b) :: ncoef_c,ncoef_r,l,m,ilm,ilm1,ilm2


    ncoef_c = (lmax+1)*(lmax+2)/2  
    ncoef_r = (lmax+1)*(lmax+1)

    ilm = 0
    ilm2 = 0
    do l = 0,lmax
       
       ilm1 = ilm2+1
       ilm2 = ilm1+2*l
       
       do m = 0,l
          ilm = ilm+1

          if(m == 0) then
             flm_c(ilm) = flm_r(ilm1+l)
          else
             flm_c(ilm) = (flm_r(ilm1+l-m) - ii*flm_r(ilm1+l+m))/sqrt2             
          end if


       end do       
    end do

    return
  end subroutine real_to_complex_SH


 subroutine complex_to_real_SH(lmax,flm_c,flm_r)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension(:), intent(in) :: flm_c
    real(dp), dimension(:), intent(out) :: flm_r

    integer(i4b) :: ncoef_c,ncoef_r,l,m,ilm,ilm1,ilm2


    ncoef_c = (lmax+1)*(lmax+2)/2  
    ncoef_r = (lmax+1)*(lmax+1)

    ilm = 0
    ilm2 = 0
    do l = 0,lmax
       
       ilm1 = ilm2+1
       ilm2 = ilm1+l
       
       do m = -l,l
          ilm = ilm+1
          if(m < 0) then
             flm_r(ilm) =  sqrt2*real(flm_c(ilm1-m))
          else if(m == 0) then
             flm_r(ilm) =  real(flm_c(ilm1))
          else if(m > 0) then
             flm_r(ilm) = -sqrt2*aimag(flm_c(ilm1+m))
          end if

       end do       
    end do

    return
  end subroutine complex_to_real_SH



  !================================================================!
  !================================================================!
  !                 my spherical harmonic routines                 !
  !================================================================!
  !================================================================!



  subroutine sphsum(th,ph,f_lm,f)
    ! this routine sums the spherical harmonic series
    ! for the real-valued function f with coefficients
    ! f_lm at the point (th,ph)
    ! 
    ! this routine first calls  get_ylm to compute the 
    ! necessary spherical harmonics, and then calls
    ! sphsum_sph to sum the series
    use nrtype
    implicit none
    real(dp), intent(in) :: th
    real(dp), intent(in) :: ph
    complex(dpc), dimension(:), intent(in) :: f_lm
    real(dp), intent(out) :: f
    
    ! compute the necessary spherical harmonics
    call get_ylm(th,ph)


    ! sum the series
    call sphsum_sph(f_lm,f)

    return
  end subroutine sphsum


  subroutine sphsum_sph(f_lm,f)
    ! routine sums the spherical harmonic series for the 
    ! real-valued function f with coeffcients f_{lm} 
    ! 
    ! This routine assumes that the spherical harmonics
    ! at the desired evaluation point have been evaluated
    ! and stored in the module-array ylm_sph
    use nrtype
    implicit none
    complex(dpc), dimension(:), intent(in) :: f_lm
    real(dp), intent(out) :: f
    
    integer(i4b) :: l,m,ilm
    
    
    f = 0.0_dp
    ilm = 0
    do l = 0,lmax_sph
       
       ! do the m = 0 term
       ilm = ilm+1
       f = f + real(f_lm(ilm)*ylm_sph(ilm))

       ! do the m > 0 terms
       do m = 1,l
          ilm = ilm+1
          f = f + 2.0_dp*real(f_lm(ilm)*ylm_sph(ilm))
       end do
       
    end do
    

    return
  end subroutine sphsum_sph


  subroutine sphsum_mesh(f_lm,f)
    ! this routine sums the spherical harmonic series
    ! for the real valued function f with coefficients 
    ! f_lm and returns the results at every node of the
    ! spherical mesh
    use nrtype
    implicit none
    complex(dpc), dimension(:), intent(in) :: f_lm
    real(dp), dimension(:), intent(out) :: f


    integer(i4b) :: irec,l,m,ilm,ithd,iphd,im


    do irec = 1,nrec_sph

       ithd = rnkth_sph(irec)
       iphd = rnkph_sph(irec)


       f(irec) = 0.0_dp
       ilm = 0
       do l = 0,lmax_sph

          ! do m = 0 term
          ilm = ilm+1
          f(irec) = f(irec) + real(f_lm(ilm)*xlm_sph(ithd,ilm))

          ! do m > 0 terms
          do m = 1,l
             im = m+1
             ilm = ilm+1
             f(irec) = f(irec) + 2.0_dp*real(f_lm(ilm) & 
                               * xlm_sph(ithd,ilm)     & 
                               *  eimp_sph(iphd,im))
          end do
          
          
       end do

    end do

    return
  end subroutine sphsum_mesh


  subroutine sphexp_cal(f_sph,f_lm)
    ! given a real valued function f tabulated 
    ! at the modes of the spherical mesh, this
    ! routine returns the spherical harmonic 
    ! expansion coefficients of the function    
    use nrtype
    implicit none
    real(dp), dimension(:), intent(in) :: f_sph
    complex(dpc), dimension(:), intent(out) :: f_lm

    integer(i4b) :: inode,jnode,ispec,irec,ithd,iphd,ilm,l,m

    ! initialize the coefficient vector
    f_lm = 0.0_dp


    ! perform integration over the grid
    irec = 0
    do ispec = 1,nspec_sph
       do inode = 1,ngll_sph
          do jnode = 1,ngll_sph
             irec = irec+1
             ithd = rnkth_sph(irec)
             iphd = rnkph_sph(irec)
             
             ilm = 0
             do l = 0,lmax_sph
                
                do m = 0,l

                   ilm = ilm+1
                   f_lm(ilm) = f_lm(ilm) + f_sph(irec)        & 
                                           *xlm_sph(ithd,ilm) & 
                                           *sinth_sph(ithd)   & 
                                           *conjg(eimp_sph(iphd,m+1))  &
                                           *wgll_sph(inode)   & 
                                           *wgll_sph(jnode)   & 
                                           *jac_sph

                end do

             end do
          end do
          
       end do
    end do

    
 
    return
  end subroutine sphexp_cal


  subroutine sphcoef_diff(f_lm)
    use nrtype
    implicit none
    complex(dpc), dimension(:), intent(out) :: f_lm
    
    integer(i4b) :: ilm,l,m
    real(dp) :: fac
    
    ilm = 0
    do l = 0,lmax_sph
       fac = twopi_d*(l+1)/(lmax_sph+0.5_dp)
       fac = exp(-fac)
       do m = 0,l
          ilm = ilm+1
          f_lm(ilm) = f_lm(ilm)*fac
       end do
    end do



    return
  end subroutine sphcoef_diff



  subroutine sph_polar_grid(lmax,verbin)
    ! given a value of lmax, this routine sets up
    ! a spherical mesh with the appropriate node 
    ! spacing for use in evaluating spherical harmonic
    ! transformations. The routine also calculates all
    ! the spherical harmonics upto degree lmax at each 
    ! node of the mesh
    use nrtype
    use module_util
    use module_function
    implicit none
    integer(i4b), intent(in) :: lmax
    integer(i4b), intent(in), optional :: verbin

    integer(i4b) :: ith,iph,irec,nth,nph,inode,ispec,jnode, & 
         l,m,ilm1,ilm2,ithd,iphd,verb
    real(dp) :: th1,th2,ph1,ph2,dth,dph,th11,th22,ph11,ph22, & 
         th,ph,sinth
    real(dp), dimension(lmax+1) :: xp_tmp,xc_tmp


    if(present(verbin)) then
       verb = verbin
    else
       verb = 0
    end if

    ! set some of the grid parameters
    lmax_sph = lmax
    ncoef_sph = (lmax_sph+1)*(lmax_sph+2)/2

    ! set the spacings for the element edges
    dth = twopi_d/((lmax_sph+0.5)*npw_sph)
    dph = dth



    ! compute the Jacobian
    jac_sph = 0.25_dp*dth*dph
    
    ! get the total number of elements
    th1 = 0.0_dp
    th2 = pi_d
    nth = (th2-th1)/dth+2
    dth = (th2-th1)/(nth-1)
    ph1 = 0.0_dp
    ph2 = twopi_d
    nph = (ph2-ph1)/dph
    dph = (ph2-ph1)/(nph-1)

    ! save some valus
    nth_sph = nth
    nph_sph = nph
    dth_sph = dth
    dph_sph = dph
    
    ! get the total number of elements 
    nspec_sph = (nth-1)*(nph-1)
    
    ! allocate the vector grid arrays
    nrec_sph = ngll_sph*ngll_sph*nspec_sph
    allocate(tha_sph(nrec_sph),pha_sph(nrec_sph))

    ! get the GLL points and weights
    call zwgljd(xigll_sph,wgll_sph,ngll_sph,0.0_dp,0.0_dp)
    if(mod(ngll_sph,2) /= 0) xigll_sph((ngll_sph-1)/2+1) = 0.0_dp    

    ! build up the grid
    if(verb == 1) print *, 'build the mesh points'
    ispec = 0
    irec  = 0
    do ith = 1,nth-1
       th11 = th1+(ith-1)*dth
       th22 = th11+dth
       do iph = 1,nph-1
          ph11 = ph1+(iph-1)*dph
          ph22 = ph11+dph
          ispec = ispec+1
          do inode = 1,ngll_sph
             th = th11+0.5_dp*(th22-th11)*(xigll_sph(inode)+1.0_dp)
             do jnode = 1,ngll_sph
                ph = ph11+0.5_dp*(ph22-ph11)*(xigll_sph(jnode)+1.0_dp)
                irec = irec+1
                tha_sph(irec) = th
                pha_sph(irec) = ph
             end do
          end do
       end do
    end do

    if(verb == 1) print *, ' nspec = ',nspec_sph
    if(verb == 1) print *, ' nrec  = ',nrec_sph

    ! construct the reduced th and ph arrays
    if(verb == 1) print *, ' indexing the reduced arrays'
    allocate(indth_sph(nrec_sph),rnkth_sph(nrec_sph))
    call reduce_array(nrec_sph,tha_sph,indth_sph,nthd_sph, & 
         thd_sph,rnkth_sph,dth/20.0_dp)
    allocate(indph_sph(nrec_sph),rnkph_sph(nrec_sph))
    call reduce_array(nrec_sph,pha_sph,indph_sph,nphd_sph, & 
         phd_sph,rnkph_sph,dph/20.0_dp)

    if(verb == 1) print *, ' nthd = ',nthd_sph
    if(verb == 1) print *, ' nphd = ',nphd_sph

    ! compute the legendre polynomials including sin(th) term
    if(verb == 1) print *, ' computing the Legendre functions'
    allocate(sinth_sph(nthd_sph))
    allocate(xlm_sph(nthd_sph,ncoef_sph))
    allocate(xplm_sph(nthd_sph,ncoef_sph))
    allocate(xclm_sph(nthd_sph,ncoef_sph))
    ilm2 = 0
    do l = 0,lmax_sph
       ilm1 = ilm2+1
       ilm2 = ilm1+l
       do ithd = 1,nthd_sph
          th = thd_sph(ithd)
          if(l == 0) sinth_sph(ithd) = sin(th)
          call legendre(th,l,l,  xlm_sph(ithd,ilm1:ilm2), & 
                                xplm_sph(ithd,ilm1:ilm2), & 
                                xclm_sph(ithd,ilm1:ilm2))          
       end do
    end do

    ! compute the complex exponentials
    if(verb == 1) print *, ' computing the complex exponentials'
    allocate(eimp_sph(nphd_sph,lmax_sph+1))
    do iphd = 1,nphd_sph
       ph = phd_sph(iphd)
       do m = 0,lmax_sph
          eimp_sph(iphd,m+1) = exp(ii*m*ph)
       end do
    end do
      

    ! allocate arrays for spherical harmonic interpolation
    allocate(ylm_sph(ncoef_sph))
    allocate(dylmdt_sph(ncoef_sph))
    allocate(dylmdp_sph(ncoef_sph))


    return
  end subroutine sph_polar_grid


  subroutine get_ylm_irec(irec,iopin)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: irec
    integer(i4b), intent(in), optional :: iopin

    integer(i4b) :: iop,ithd,iphd,l,m,ilm,im

    if(present(iopin)) then
       iop = iopin
    else
       iop = 0
    end if
    
    ithd = rnkth_sph(irec)
    iphd = rnkph_sph(irec)

    ilm = 0
    do l = 0,lmax_sph
       do m = 0,l
          im = m+1
          ilm = ilm+1
          ylm_sph(ilm) = xlm_sph(ithd,ilm)*eimp_sph(iphd,im)
          if(iop == 1) then
             dylmdt_sph(ilm) = xplm_sph(ithd,ilm)*eimp_sph(iphd,im)
             dylmdp_sph(ilm) = ii*m*xclm_sph(ithd,ilm)*eimp_sph(iphd,im)
          end if
       end do
    end do


    return
  end subroutine get_ylm_irec
  

  subroutine get_ylm(th,ph,iopin)
    use nrtype
    implicit none
    real(dp), intent(in) :: th
    real(dp), intent(in) :: ph
    integer(i4b), intent(in), optional :: iopin

    integer(i4b) :: ispec,irec,iop,ithd,iphd,l,m,ilm,im,inode,jnode

    real(dp) :: xi,eta

    if(present(iopin)) then
       iop = iopin
    else
       iop = 0
    end if
    
    ! find the point in the mesh
    call sph_mesh_finder(th,ph,ispec,xi,eta)

    ! get the Lagrange polynomials
    call lagrange_any(xi,ngll_sph,xigll_sph,htmp_xi_sph,hptmp_xi_sph)
    call lagrange_any(eta,ngll_sph,xigll_sph,htmp_eta_sph,hptmp_eta_sph)

    ilm = 0
    do l = 0,lmax_sph
        do m = 0,l
          im = m+1
          ilm = ilm+1
          ylm_sph(ilm) = 0.0_dp
          if(iop == 1) then
             dylmdt_sph(ilm) = 0.0_dp
             dylmdp_sph(ilm) = 0.0_dp
          end if
          irec = (ispec-1)*ngll_sph*ngll_sph
          do inode = 1,ngll_sph
             do jnode = 1,ngll_sph
                irec = irec+1
                ithd = rnkth_sph(irec)
                iphd = rnkph_sph(irec)
                ylm_sph(ilm) = ylm_sph(ilm) + xlm_sph(ithd,ilm)   & 
                                            * eimp_sph(iphd,im)   & 
                                            * htmp_xi_sph(inode)  & 
                                            * htmp_eta_sph(jnode)
                if(iop == 1) then
                   dylmdt_sph(ilm) = dylmdt_sph(ilm) + xplm_sph(ithd,ilm)  & 
                                                     * eimp_sph(iphd,im)   & 
                                                     * htmp_xi_sph(inode)  & 
                                                     * htmp_eta_sph(jnode)
                   dylmdp_sph(ilm) = dylmdp_sph(ilm) + ii*m                & 
                                                     * xclm_sph(ithd,ilm)  & 
                                                     * eimp_sph(iphd,im)   & 
                                                     * htmp_xi_sph(inode)  & 
                                                     * htmp_eta_sph(jnode)

                end if
                
             end do
          end do


       end do
    end do


    return
  end subroutine get_ylm




  subroutine sph_mesh_finder(th,ph,ispec,xi,eta)
    use nrtype
    implicit none
    real(dp), intent(in) :: th
    real(dp), intent(in) :: ph
    integer(i4b), intent(out) :: ispec
    real(dp), intent(out) :: xi 
    real(dp), intent(out) :: eta
    
    integer(i4b) :: ith,iph
    real(dp) :: th11,th22,ph11,ph22

    ith = floor(th/dth_sph)+1
    iph = floor(ph/dph_sph)+1
    

    if(ith == nth_sph) ith = ith-1
    if(iph == nph_sph) iph = iph-1

    ! work out bounding theta and phi
    th11 = (ith-1)*dth_sph
    th22 = th11+dth_sph
    ph11 = (iph-1)*dph_sph
    ph22 = ph11+dph_sph


    ! get the element number
    ispec = (ith-1)*(nph_sph-1) + iph 

    ! work out the xi and eta co-ordinates
    xi  = 2.0_dp*(th-th11)/(th22-th11)-1.0_dp
    eta = 2.0_dp*(ph-ph11)/(ph22-ph11)-1.0_dp


    return
  end subroutine sph_mesh_finder


  subroutine rotate_from_pole(lat,lon,alpha,beta,gamma)
      ! given the lat and lon of a point in the earth model, 
      ! this routine returns the Euler angles that rotates the North pole 
      ! to this point in the new co-ordinate system and aligns the 
      ! great circle between the point and the geographic north pole
      ! with the prime meridian of the new co-ordinate system
      use nrtype
      implicit none
      real(dp), intent(in) :: lat,lon
      real(dp), intent(out) :: alpha,beta,gamma


      alpha = pi_d
      beta  = (90.0_dp-lat)*deg2rad
      gamma = (180.0_dp - lon)*deg2rad

      
      return
    end subroutine rotate_from_pole


  subroutine rotvm0adj(alpha,beta,gamma,lmax,A,B)
    ! given the Euler angles alpha,beta, and gamma 
    ! as in Edmonds, this routine rotates the adjoint
    ! spherical harmonic coefficients A
    ! putting them in B
    !
    ! the input coefficients are only non-zero for m=0
    !
    ! see (C.283) of Dahlen & Tromp (1998) for the
    ! transformation formulae
    
    use nrtype
    use module_model
    implicit none
    real(dp), intent(in) :: alpha,beta,gamma
    integer(i4b), intent(in) :: lmax
    complex(dpc), dimension(lmax), intent(in) :: A
    complex(dpc), dimension(-lmax:lmax,lmax), intent(out) :: B
    
    integer(i4b) :: l,m,n
    real(dp), dimension(2*lmax+1,2*lmax+1) :: ldnm
    complex(dpc), dimension(2*lmax+1,2*lmax+1) :: bdnm
    
    ldnm = 0.0_dp
    bdnm = 0.0_dp
    B = 0.0_dp

    ! As input coefficients are only non-zero for m = 0
    n = 0
       
    do l = 1,lmax

       ! compute the Wigner d functions for the rotation matrix
       
       call rotmx2(l,l,beta,ldnm,2*l+1,2*l+1)
       
       ! put together the rotation matrix
       do m = -lmax,lmax
          if (m > l) exit
          if (abs(m) > l) cycle
          bdnm(l+m+1,l+n+1) = exp(ii*m*gamma)*ldnm(l+m+1,l+n+1) & 
                              *exp(ii*n*alpha)
       end do
    
       ! perform rotation of the coefficients
       do m = -lmax,lmax  
          if (m > l) exit
          if (abs(m) > l) cycle
          B(m,l) = bdnm(l+m+1,l+n+1)*A(l)
       end do

    end do

    B = conjg(B)
    
    return
  end subroutine rotvm0adj





end module module_sphexp
