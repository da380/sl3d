program adj_test

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
  complex(dpc), dimension(:,:), allocatable :: sea_level_data

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec,ispec1,nt,it1,it,tmp1,tmp2
  complex(dpc) :: tmp3

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,load_lat,load_long

  complex(dpc), dimension(:,:,:), allocatable :: sl_fin,ice_kern
  complex(dpc), dimension(:,:,:), allocatable :: UVp_adj,Q_adj,W_adj
  complex(dpc), dimension(:,:,:,:,:,:,:), allocatable :: mem
  complex(dpc), dimension(:,:), allocatable :: ice_load,ice_adj,sl_adj
  complex(dpc), dimension(:), allocatable :: coefs,ps,Us
  complex(dpc), dimension(:), allocatable :: K_adj

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

  call visc_3d

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
  open(newunit=io,file='/raid/oc251/data/ice5g/ascii/ice5g_v1.2_21.0k_1deg.orog.ascii')
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
  t2 = 500.0_dp*yr2sec/t_norm
  dt = 50.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 11
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_fin(ngl,nphi,nt))
  sl_fin = 0.0_dp

  ! Forward calculation
  call forward_full_sl(nt,dt,nt,ice_data,sea_level_data,sl_fin)

  open(newunit=io,file='sltestnewi.1')
  do iphi = 1,nphi
     do igl = 1,ngl
        !read(io,*) tmp1,tmp2,tmp3
        !sl_t(igl,iphi,nt) = tmp3
       write(io,*) iphi,igl,real(sl_fin(igl,iphi,nt))-real(sl_fin(igl,iphi,1))
     end do
     write(io,*)
  end do
  close(io)
  stop


  !============================================!
  !            Adjoint calculation             !
  !============================================!

  allocate(UVp_adj((lmax+1)**2,nglob_ssg,2),Q_adj(ngl,nphi,2))
  allocate(W_adj((lmax+1)**2,nglob_tor,2))
  allocate(K_adj(2))
  allocate(ice_load(ngl,nphi))
  allocate(mem(5,ngl,nphi,ntau,ngll,nspec,2))
  allocate(Us((lmax+1)**2),ps((lmax+1)**2))
  allocate(sl_adj(ngl,nphi),ice_adj(ngl,nphi))
  allocate(ice_kern(ngl,nphi,nt))
  UVp_adj = 0.0_dp
  W_adj = 0.0_dp
  mem = 0.0_dp
  Q_adj = 0.0_dp
  K_adj = 0.0_dp
  ice_load = 0.0_dp
  Us = 0.0_dp
  ps = 0.0_dp
  sl_adj = 0.0_dp
  ice_adj = 0.0_dp
  ice_kern = 0.0_dp

  ! Set load
  allocate(adj_load%load_lm((lmax+1)**2))
  adj_load%type = 1

  load_lat = 0.0_dp
  load_long = 0.0_dp

  call point_load(load_lat,load_long,0,adj_load%load_lm)

  ! Get ice coverage
  call time_point_interp(nt,dt,ice_data,ice_load)

  call stat_adj(ice_load,sl_fin(:,:,nt),adj_load,UVp_adj(:,:,1),Q_adj(:,:,1),K_adj(1))


  ps = 0.0_dp
  Us = 0.0_dp
  call surface_fields(UVp_adj(:,:,1),ps,Us)

  call adj_sl_ice(ice_load,sl_fin(:,:,nt),adj_load,Us,ps,Q_adj(:,:,1),K_adj(1),sl_adj,ice_adj)

  call kern_ice_it(ice_load,sl_fin(:,:,nt), &
                   Us,ps,sl_adj,ice_adj,ice_kern(:,:,nt))

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

     call adj_sl_ice(ice_load,sl_fin(:,:,nt-it+1),adj_load,Us,ps,Q_adj(:,:,2),K_adj(2),sl_adj,ice_adj)

     ! Calculate kernel
     call kern_ice_it(ice_load,sl_fin(:,:,nt-it+1),Us,ps, &
                      sl_adj,ice_adj,ice_kern(:,:,nt-it+1))

     if (it == 21) then
        open(newunit=io,file='ikerntestnew.21')
        do iphi = 1,nphi
           do igl = 1,ngl
              write(io,*) iphi,igl,ice_kern(igl,iphi,nt-it+1)
           end do
           write(io,*)
        end do
        close(io)
        stop
     end if

     UVp_adj(:,:,1) = UVp_adj(:,:,2)
     W_adj(:,:,1) = W_adj(:,:,2)
     mem(:,:,:,:,:,:,1) = mem(:,:,:,:,:,:,2)
     Q_adj(:,:,1) = Q_adj(:,:,2)
     K_adj(1) = K_adj(2)

  end do


end program adj_test
