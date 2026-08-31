program sea_level_3d

  use nrtype
  use module_util
  use module_function
  use module_gia
  use module_fourier
  use module_model
  use module_mat
  use module_sl
  use module_ice
  use module_config

  implicit none

  integer(i4b) :: narg,io,io2,ilayer

  logical(lgt) :: ltmp
  character(len=256) :: model_file,file1,icefile,file2

  real(dp) :: drmax,lambda_l,rr

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data,sea_level_sm
  complex(dpc), dimension(:,:), allocatable :: ice_loaddot,ice_load,ice_load2
  complex(dpc), dimension(:), allocatable :: sl_lm,water_lm,iceload_lm

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,it1,it1_old,icount
  real(dp) :: tmp,long,lat,xx,yy,age
  integer(i4b) :: iage

  integer(i4b) :: nt,it,l,ispec1,lmcur,ia,inode,ispec,m,ua,nout,pa,lmc
  real(dp) :: t1,t2,dt,rl,sltmp

  complex(dpc), dimension(:,:,:), allocatable :: sl_out
  complex(dpc), dimension(:,:,:), allocatable :: UVp_out

  !------------------------------------!
  !  set the normalization parameters  !
  !------------------------------------!

  call set_parameters

  !----------------------------------------------!
  !  the parameters of this run: these defaults  !
  !  may be overridden by a config file          !
  !----------------------------------------------!
  lmax_run     = 64
  time_step    = 5.0_dp      ! years
  output_every = 50.0_dp      ! years

  call read_config

  ! leave a record of what this run used, beside its results
  call write_config('config.out.nml')

  !--------------------------!
  !  read in the model file  !
  !--------------------------!
  open(newunit = io,file=trim(earth_model),action='read',form='formatted')
  call read_model(io)
  close(io)

  !---------------------------------!
  !      calculate mesh size        !
  !---------------------------------!

  ! set lmax
  lmax = lmax_run

  call set_sh_grid
  call calc_sph_harm

  call rob_projection

  lambda_l = twopi_d*r(nknot)/(lmax+0.5_dp)
  drmax = lambda_l/2.0_dp
  call mesh_model(drmax)

  open(newunit=io,file='model.out')
  do ispec = 1,nspec
     do inode = 1,ngll
        write(io,*) ispec, r_node(inode,ispec)*r_norm/1000.0_dp,rho_node(inode,ispec)*rho_norm, &
                    mu_node(inode,ispec)*con_norm, kappa_node(inode,ispec)*con_norm
     end do
  end do
  close(io)

  call global_points

  call visc_3d_from_vs
  !call visc_3d


  open(newunit=io,file='visc.s20.p01')
  do ilayer = 1,nsect
     if (sect_ind(3,ilayer) /= 2) cycle
     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
        do inode = 1,ngll

           rr = r_node(inode,ispec)*r_norm


           do iphi = 1,nphi
              long = aphi(iphi)*rad2deg - 180.0_dp

              do igl = 1,ngl
                 lat = (pio2_d - tgl(igl))*rad2deg
                 call get_xy(lat,long,xx,yy)

                 write(io,*) rr,long,lat,xx,yy, &
                      real(mui_node(1,inode,ispec)*visco_norm/si_spat_node(igl,iphi,1,inode,ispec))!, &
                      !imag(mui_node(1,inode,ispec)*visco_norm/si_spat_node(igl,iphi,1,inode,ispec))
              end do
              write(io,*) rr, "ABC"
              write(io,*) "ABC", long
           end do
             
        end do
     end do
  end do
  close(io)

!  stop


  !=========================================!
  !             calculate load              !
  !=========================================!

  allocate(ice_data(ngl,nphi,nt_ice))

  call set_ice_model(ice_data)    


  !==================================================!
  !          Calculate initial sea level             !
  !==================================================!

  ! set initial sea level to be negative topography
  allocate(sea_t(nlong_ice,nlat_ice))
  sea_t = 0.0_dp
  ! the surface altitude at 21 ka, as negative topography
  call read_ice_grid(1,'topography',sea_t)
  sea_t = -sea_t/r_norm

  allocate(sea_level_data(ngl,nphi))
  sea_level_data = 0.0_dp
  do iphi = 1,nphi
     long = aphi(iphi)*rad2deg - 180.0_dp
     do igl = 1,ngl
        lat = (pio2_d - tgl(igl))*rad2deg
        sea_level_data(igl,iphi) = interp(sea_t(:,:),lat,long)
        sea_level_data(igl,iphi) = sea_level_data(igl,iphi) + ice_data(igl,iphi,1)
     end do
  end do

  allocate(sl_lm((lmax+1)**2))
  sl_lm = 0.0_dp

  call coefs_from_fun(sea_level_data,sl_lm(:))

  
  !===========================================!
  !             GIA calculation               !
  !===========================================!

  ! set time stepping parameters
  t1 = 0.0_dp
!  t2 = 3000.0_dp*yr2sec/t_norm
  ! run over the whole of the chosen ice history
  t2 = run_time()*yr2sec/t_norm
  !dt = 5.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = nint(run_time()/time_step) + 1
!  nt = 601
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  !nout = t2*t_norm*sec2yr/1000.0_dp + 1
  nout = nint(run_time()/output_every) + 1
!  nout = 61
  !print *, nout
  !nout = nt

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_out(ngl,nphi,nout))
  !allocate(UVp_out((lmax+1)**2,nglob_ssg,nout))


  call forward_full_sl(nt,dt,nout,ice_data,sea_level_data,sl_out)

  !do it = 1,nout
  !   call string_cat_int('grace.64.l2.m22.',21-it/2,file1)
  !   if (mod(it,2) == 0) then
  !      file2 = trim(file1) // trim('.5')
  !   else
  !      file2 = trim(file1) // trim('.0')
  !   end if
  !   open(newunit=io,file=file2)
  !   do l = 2,lmax
  !      do m = -l,l
  !         if (m < 0) then
  !            lmc = l**2 + l - m + 1
  !            write(io,*) l,m,sqrt(2.0_dp)*real(UVp_out(lmc,pa_surf(l),it))*pot_norm
  !         else if (m==0) then
  !            lmc = l**2 + l + 1
  !            write(io,*) l,m,real(UVp_out(lmc,pa_surf(l),it))*pot_norm
  !         else
  !            lmc = l**2 + l + m + 1
  !            write(io,*) l,m,-sqrt(2.0_dp)*imag(UVp_out(lmc,pa_surf(l),it))*pot_norm
  !         end if
  !      end do
  !   end do
  !   close(io)
  !end do
  !stop
   
!  open(newunit=io,file='sl.64.20.p01.full')
  do it = 1,nout
     ! name each file by the age it holds, in ka, taken from the run's own
     ! time base rather than assumed, so the names stay right whichever ice
     ! model was chosen and however often the run writes output
     age = (ice_total_time() &
            - (it-1)*(t2-t1)*t_norm*sec2yr/max(nout-1,1))/1000.0_dp
     iage = nint(100.0_dp*age)          ! the age in hundredths of a ka, so
     write(file2,'(a,i0,a,i2.2)') &     ! that output written every ten years
          'sl.3d.',iage/100,'.',mod(iage,100)     ! or less often is distinct
     open(newunit=io,file=file2)
     do iphi = 1,nphi
        !long = aphi(iphi)*rad2deg
        long = aphi(iphi)*rad2deg - 180.0_dp
        do igl = 1,ngl
           lat = (pio2_d - tgl(igl))*rad2deg
           call get_xy(lat,long,xx,yy)
           write(io,*) it,long,lat,xx,yy,real(sl_out(igl,iphi,it))*r_norm, &
                      (real(sl_out(igl,iphi,it)) - real(sl_out(igl,iphi,1)))*r_norm
        end do
        write(io,*)
     end do
     write(io,*)
     close(io)
  end do
!  close(io)
  !stop


end program sea_level_3d
