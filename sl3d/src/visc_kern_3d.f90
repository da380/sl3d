program visc_kern_3d

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
  character(len=256) :: model_file,file1,icefile,file2

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec,ispec1,nt,it1,it,nout,ilayer,inode,ua

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,load_lat,load_long,rr,rl

  complex(dpc), dimension(:,:,:), allocatable :: sl_fin
  complex(dpc), dimension(:,:,:), allocatable :: UVp,W
  complex(dpc), dimension(:,:,:,:,:,:,:), allocatable :: mem
  complex(dpc), dimension(:,:,:,:), allocatable :: kern

  type(meas) :: adj_load

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
  site_lat     = 57.0_dp
  site_long    = -110.0_dp
  run_length   = 18000.0_dp

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

  ispec1 = 1
  call global_points

  call visc_3d_from_vs

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
  !nout = nt

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_fin(ngl,nphi,nt))
  allocate(UVp((lmax+1)**2,nglob_ssg,nout))
  allocate(W((lmax+1)**2,nglob_tor,nout))
  allocate(mem(5,ngl,nphi,ntau,ngll,nspec,nout))

  sl_fin = 0.0_dp
  UVp = 0.0_dp
  W = 0.0_dp
  mem = 0.0_dp


  ! Forward calculation
  call forward_full_sl(nt,dt,nt,ice_data,sea_level_data,sl_fin,UVp,W,mem,nout)


  !============================================!
  !            Adjoint calculation             !
  !============================================!

  allocate(kern(ngl,nphi,ngll,nspec))
  kern = 0.0_dp

  ! Set load
  allocate(adj_load%load_lm((lmax+1)**2))
  adj_load%type = 1

  ! Not bulge
  !load_lat = 38.4719561826_dp
  !load_long = -99.84375_dp

  ! Bulge
  !load_lat = 43.254194665350937
  !load_long = -104.0625

  ! Tahiti - true
  !load_long = -149.426_dp
  !load_lat = -17.6509_dp

  ! North Canada
  load_long = site_long
  load_lat = site_lat
  !load_lat = 0.0_dp
  !load_long = 0.0_dp

  call point_load(load_lat,load_long,0,adj_load%load_lm)

  call adj_kern_visc_full(nt,dt,nout,ice_data,sl_fin,adj_load,UVp,W,mem,kern)


  open(newunit = io,file = 'vkncan.s20.p01.interp.18')
  do ilayer = 1,nsect
     
     ! No terms in fluid or elastic regions
     if (sect_ind(3,ilayer) /= 2) cycle

     ! begin loop over the spectral elements
     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
                    
        do inode = 1,ngll
           rr = r_node(inode,ispec)*r_norm
           call get_depth(1,rl,ispec1)
           ua = gvn_ssg(inode,ispec,ispec1,2)
           do iphi = 1,nphi
              long = aphi(iphi)*rad2deg - 180.0_dp
              do igl = 1,ngl
                 lat = (pio2_d - tgl(igl))*rad2deg
                 call get_xy(lat,long,xx,yy)
                 write(io,*) rr,long,lat,xx,yy,real(kern(igl,iphi,inode,ispec))
              end do
              write(io,*) rr, "ABC"
              write(io,*) "ABC", long
           end do
        end do
     end do
  end do
  close(io)


end program visc_kern_3d
