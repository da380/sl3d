program sl_visc_inv

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
  use module_kern

  implicit none

  integer(i4b) :: narg,io,io2

  logical(lgt) :: ltmp
  character(len=256) :: model_file,file1,icefile,file2,grace_file

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data

  integer(i4b) :: ispec,ispec1,nt,it1,it,nout,ilayer,inode,nco,iphi,igl,ilat,ilong,ilonged

  real(dp) :: t1,t2,dt,rr,rl,tmp,lat,long

  type(rsh_meas), dimension(:), allocatable :: inv_data
  type(comb_meas) :: comb_data


  !------------------------------------!
  !  set the normalization parameters  !
  !------------------------------------!

  call set_parameters

  !----------------------------------------------!
  !  the parameters of this run: these defaults  !
  !  may be overridden by a config file          !
  !----------------------------------------------!
  lmax_run     = 64
  time_step    = 25.0_dp      ! years
  output_every = 25.0_dp      ! years
  run_length   = 5000.0_dp

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
        write(io,*) ispec, r_node(inode,ispec),r_node(inode,ispec)*r_norm/1000.0_dp,rho_node(inode,ispec)*rho_norm, &
                 mu_node(inode,ispec)*con_norm, kappa_node(inode,ispec)*con_norm
     end do
  end do
  close(io)

  ispec1 = 1
  call global_points

  call visc_3d_none

  !=========================================!
  !             calculate load              !
  !=========================================!

  allocate(ice_data(ngl,nphi,nt_ice))
  ice_data = 0.0_dp

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

  !===========================================!
  !             GIA calculation               !
  !===========================================!

  ! set time stepping parameters
  t1 = 0.0_dp
  t2 = run_time()*yr2sec/t_norm
  dt = time_step*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = nint(run_time()/time_step) + 1
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = nint(run_time()/output_every) + 1

  ! Construct A matrices
  call construct_system_matrix

  !===========================================!
  !                   Data                    !
  !===========================================!

  ! Number of sets of spherical harmonics
  nco = 1

  allocate(inv_data(nco))

  inv_data(1)%type = 3
  inv_data(1)%nt = 1
  allocate(inv_data(1)%times(inv_data(1)%nt))
  inv_data(1)%times(1) = 16000.0_dp
  allocate(inv_data(1)%values((lmax+1)**2,inv_data(1)%nt))

  grace_file = 'grace.64.l2.m22.'
  call extract_grace(grace_file,inv_data(1))

  inv_data(1)%times(:) = inv_data(1)%times(:)*yr2sec/t_norm
  inv_data(1)%values(:,:) = inv_data(1)%values(:,:)/pot_norm
  

  comb_data%n_point = 0
  comb_data%n_rsh = nco
  comb_data%pres_point = .false.
  comb_data%pres_rsh = .true.
  allocate(comb_data%rsh(comb_data%n_rsh))
  comb_data%rsh = inv_data

  !===========================================!
  !                 Inversion                 !
  !===========================================!

  file1 = 'grace.lim.16.'


  ! DA comment to remove compilation error
!  call visc_inversion(nt,dt,nt,ice_data,sea_level_data,comb_data,file1)

end program sl_visc_inv
