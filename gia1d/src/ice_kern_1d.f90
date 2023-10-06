program ice_kern_1d

  use nrtype
  use module_util
  use module_function
  use module_gia
  use module_fourier
  use module_model
  use module_mat
  use module_sl
  use module_kern

  implicit none

  integer(i4b) :: narg,io,io2

  logical(lgt) :: ltmp
  character(len=256) :: model_file,file1,icefile,file2

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data,ice_load

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec,ispec1,nt,it1,it,tmp1,tmp2,nout
  complex(dpc) :: tmp3

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,load_lat,load_long

  complex(dpc), dimension(:,:,:), allocatable :: sl_fin,ice_kern

  type(meas) :: adj_load

  !------------------------------------!
  !  set the normalization parameters  !
  !------------------------------------!

  call set_parameters

  ! number of input arguments 
  narg = command_argument_count()

  ! if needed, display input info
  if(narg == 0) then
     print *, 'inputs: [model]'
     stop
  end if

  if(narg /= 1) stop 'bad input'

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
  lmax = 128

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

  allocate(ice_data(ngl,nphi,nt_ice5g))
  ice_data = 0.0_dp

  call set_ice_model(ice_data)    


  !==================================================!
  !          Calculate initial sea level             !
  !==================================================!

  ! set initial sea level to be negative topography
  allocate(sea_t(nlong_ice5g,nlat_ice5g))
  sea_t = 0.0_dp
  open(newunit=io,file='/home/oc251/raid/data/ice5g/ascii/ice5g_v1.2_21.0k_1deg.orog.ascii')
  do ilat = 1,nlat_ice5g
     do ilong = 1,nlong_ice5g
        !ilonged = ilong
        ilonged = ilong + (nlong_ice5g/2)
        if (ilonged > nlong_ice5g) ilonged = ilonged - nlong_ice5g
        read(io,*) tmp
        sea_t(ilonged,ilat) = -1.0_dp*tmp/(r_norm)
     end do
  end do
  close(io)

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
  t2 = 21000.0_dp*yr2sec/t_norm
  dt = 50.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 421
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = 22

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_fin(ngl,nphi,nt))
  sl_fin = 0.0_dp

  ! Forward calculation
  call forward_full_sl(nt,dt,nt,ice_data,sea_level_data,sl_fin)


  !============================================!
  !            Adjoint calculation             !
  !============================================!

  allocate(ice_kern(ngl,nphi,nt))
  ice_kern = 0.0_dp

  ! Set load
  allocate(adj_load%load_lm((lmax+1)**2))
  adj_load%type = 1


  ! Bulge
  !load_lat = 38.4719561826_dp
  !load_long = -99.84375_dp

  ! Tahiti
  load_lat = -18.138970990239354_dp
  load_long = -149.0625_dp

  call point_load(load_lat,load_long,lmax,adj_load%load_lm)

  call adj_kern_ice_full(nt,dt,nout,ice_data,sl_fin,adj_load,ice_kern)

  allocate(ice_load(ngl,nphi))

  do it = 1,nout
     call time_point_interp(1 + (nt-1)*(it-1)/(nout-1),dt,ice_data,ice_load)
     call string_cat_int('iktah.128.1d.sm.proj.',it-1,file1)
     open(newunit=io,file=file1)
     do iphi = 1,nphi
        long = aphi(iphi)*rad2deg - 180.0_dp
        do igl = 1,ngl
           lat = (pio2_d - tgl(igl))*rad2deg
           call get_xy(lat,long,xx,yy)
           if (real(ice_load(igl,iphi)) == 0.0_dp) then
              write(io,*) long,lat,xx,yy,0.0_dp
           else
              write(io,*) long,lat,xx,yy,real(ice_kern(igl,iphi,it))
           end if
        end do
        write(io,*)
     end do
     close(io)
  end do


end program ice_kern_1d
