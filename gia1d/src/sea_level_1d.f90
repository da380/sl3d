program sea_level_1d

  use nrtype
  use module_util
  use module_function
  use module_gia
  use module_fourier
  use module_model
  use module_mat
  use module_sl
  use module_ice

  implicit none

  integer(i4b) :: narg,io,io2

  logical(lgt) :: ltmp
  character(len=256) :: model_file,file1,icefile,file2

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data,sea_level_sm
  complex(dpc), dimension(:,:), allocatable :: ice_loaddot,ice_load,ice_load2
  complex(dpc), dimension(:), allocatable :: sl_lm,water_lm,iceload_lm

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,it1,it1_old,icount
  real(dp) :: tmp,long,lat,xx,yy,age
  integer(i4b) :: iage

  integer(i4b) :: nt,it,l,ispec1,lmcur,ia,inode,ispec,m,ua,nout,pa
  real(dp) :: t1,t2,dt,rl,sltmp

  complex(dpc), dimension(:,:,:), allocatable :: sl_out

  !------------------------------------!
  !  set the normalization parameters  !
  !------------------------------------!

  call set_parameters

  ! number of input arguments 
  narg = command_argument_count()

  ! if needed, display input info
  if(narg == 0) then
     print *, 'inputs: [model] {ice model: ice5g (default), ice6g, ice7g}'
     stop
  end if

  if(narg > 2) stop 'bad input'

  ! the ice model may be named as an optional second argument
  call init_ice_model

  !--------------------------!
  !  read in the model file  !
  !--------------------------!
  call get_command_argument(1,model_file)
  inquire(file=trim(model_file),exist=ltmp)
  if(.not.ltmp) stop 'can''t find model file'

  open(newunit = io,file=trim(model_file),action='read',form='formatted')
  call read_model(io)
  close(io)

  !---------------------------------!
  !      calculate mesh size        !
  !---------------------------------!

  ! set lmax
  lmax = 64

  call set_sh_grid
  call calc_sph_harm

  call rob_projection

  lambda_l = twopi_d*r(nknot)/(lmax+0.5_dp)
  drmax = lambda_l/2.0_dp
  call mesh_model(drmax)

  open(newunit=io,file='model.out')
  do ispec = 1,nspec
     write(io,*) ispec, r_node(1,ispec)*r_norm/1000.0_dp,rho_node(1,ispec)*rho_norm, &
                 mu_node(1,ispec)*con_norm, kappa_node(1,ispec)*con_norm
  end do
  close(io)

  ispec1 = 1
  call global_points

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

  call coefs_from_fun(lmax,sea_level_data,sl_lm(:))

  
  !===========================================!
  !             GIA calculation               !
  !===========================================!

  ! set time stepping parameters: run over the whole of the chosen ice
  ! history, in steps of 50 years, writing the sea level every 500
  t1 = 0.0_dp
  t2 = ice_total_time()*yr2sec/t_norm
  nt = nint(ice_total_time()/50.0_dp) + 1
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = nint(ice_total_time()/500.0_dp) + 1

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_out(ngl,nphi,nout))

  call forward_full_sl(nt,dt,nout,ice_data,sea_level_data,sl_out)
   

  do it = 1,nout
     ! name each file by the age it holds, in ka, taken from the run's own
     ! time base rather than assumed, so the names stay right whichever ice
     ! model was chosen and however often the run writes output
     age = (ice_total_time() &
            - (it-1)*(t2-t1)*t_norm*sec2yr/max(nout-1,1))/1000.0_dp
     iage = nint(100.0_dp*age)          ! the age in hundredths of a ka, so
     write(file2,'(a,i0,a,i2.2)') &     ! that output written every ten years
          'sl.1d.',iage/100,'.',mod(iage,100)     ! or less often is distinct
     open(newunit=io,file=file2)
     do iphi = 1,nphi
        !long = aphi(iphi)*rad2deg
        long = aphi(iphi)*rad2deg - 180.0_dp
        do igl = 1,ngl
           lat = (pio2_d - tgl(igl))*rad2deg
           !call get_xy(lat,long,xx,yy)
           write(io,*) long,lat,&!xx,yy,real(r_displ(igl,iphi))*r_norm, &
                       !real(grav(igl,iphi))*pot_norm, &
                       real(sl_out(igl,iphi,it))*r_norm, &
                      (real(sl_out(igl,iphi,it)) - real(sl_out(igl,iphi,1)))*r_norm
        end do
        write(io,*)
     end do
  end do
  close(io)


end program sea_level_1d
