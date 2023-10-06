program sl_visc_inv

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
  character(len=256) :: model_file,file1,icefile,file2,sl_file

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec,ispec1,nt,it1,it,nout,ilayer,inode,ua,nloc

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,load_lat,load_long,rr,rl

  type(point_meas), dimension(:), allocatable :: inv_data
  type(comb_meas) :: comb_data


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
  t2 = 21000.0_dp*yr2sec/t_norm
  dt = 50.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 421
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = nt

  ! Construct A matrices
  call construct_system_matrix

  !===========================================!
  !                   Data                    !
  !===========================================!

  ! Number of locations
  nloc = 3

  allocate(inv_data(nloc))

  ! 1 = BG, 2 = Barb, 3 = Sunda, 4 = Richmond Gulf, 5 = James Bay
  ! 1 = BG
  inv_data(1)%type = 1
  inv_data(1)%nt = 10
  allocate(inv_data(1)%times(inv_data(1)%nt),inv_data(1)%values(inv_data(1)%nt))
  inv_data(1)%long = 127.96875000000000_dp
  inv_data(1)%lat = -10.992042206055032_dp
  inv_data(1)%times(1) = 21000.0_dp
  inv_data(1)%times(2) = 20500.0_dp
  inv_data(1)%times(3) = 20000.0_dp
  inv_data(1)%times(4) = 19500.0_dp
  inv_data(1)%times(5) = 19000.0_dp
  inv_data(1)%times(6) = 18500.0_dp
  inv_data(1)%times(7) = 18000.0_dp
  inv_data(1)%times(8) = 17500.0_dp
  inv_data(1)%times(9) = 17000.0_dp
  inv_data(1)%times(10) = 0.0_dp

  ! 2 = Barbados
  inv_data(2)%type = 1
  inv_data(2)%nt = 21
  allocate(inv_data(2)%times(inv_data(2)%nt),inv_data(2)%values(inv_data(2)%nt))
  inv_data(2)%long = -59.062500000000000_dp
  inv_data(2)%lat = 13.740049889134101_dp
  inv_data(2)%times(1) = 21000.0_dp
  inv_data(2)%times(2) = 20500.0_dp
  inv_data(2)%times(3) = 20000.0_dp
  inv_data(2)%times(4) = 19500.0_dp
  inv_data(2)%times(5) = 19000.0_dp
  inv_data(2)%times(6) = 18500.0_dp
  inv_data(2)%times(7) = 18000.0_dp
  inv_data(2)%times(8) = 15000.0_dp
  inv_data(2)%times(9) = 14500.0_dp
  inv_data(2)%times(10) = 14000.0_dp
  inv_data(2)%times(11) = 13500.0_dp
  inv_data(2)%times(12) = 13000.0_dp
  inv_data(2)%times(13) = 12500.0_dp
  inv_data(2)%times(14) = 12000.0_dp
  inv_data(2)%times(15) = 11500.0_dp
  inv_data(2)%times(16) = 11000.0_dp
  inv_data(2)%times(17) = 10500.0_dp
  inv_data(2)%times(18) = 10000.0_dp
  inv_data(2)%times(19) = 9500.0_dp
  inv_data(2)%times(20) = 9000.0_dp
  inv_data(2)%times(21) = 0.0_dp

  ! 3 = Sunda Shelf
  inv_data(3)%type = 1
  inv_data(3)%nt = 21
  allocate(inv_data(3)%times(inv_data(3)%nt),inv_data(3)%values(inv_data(3)%nt))
  inv_data(3)%long = 109.68750000000000_dp
  inv_data(3)%lat = 2.7480114865280294_dp
  inv_data(3)%times(1) = 21000.0_dp
  inv_data(3)%times(2) = 20500.0_dp
  inv_data(3)%times(3) = 20000.0_dp
  inv_data(3)%times(4) = 19500.0_dp
  inv_data(3)%times(5) = 19000.0_dp
  inv_data(3)%times(6) = 18500.0_dp
  inv_data(3)%times(7) = 18000.0_dp
  inv_data(3)%times(8) = 15000.0_dp
  inv_data(3)%times(9) = 14500.0_dp
  inv_data(3)%times(10) = 14000.0_dp
  inv_data(3)%times(11) = 13500.0_dp
  inv_data(3)%times(12) = 13000.0_dp
  inv_data(3)%times(13) = 12500.0_dp
  inv_data(3)%times(14) = 12000.0_dp
  inv_data(3)%times(15) = 11500.0_dp
  inv_data(3)%times(16) = 11000.0_dp
  inv_data(3)%times(17) = 10500.0_dp
  inv_data(3)%times(18) = 10000.0_dp
  inv_data(3)%times(19) = 9500.0_dp
  inv_data(3)%times(20) = 9000.0_dp
  inv_data(3)%times(21) = 0.0_dp

  ! 4 = Richmond Gulf
  !inv_data(4)%type = 1
  !inv_data(4)%nt = 5
  !allocate(inv_data(4)%times(inv_data(4)%nt),inv_data(4)%values(inv_data(4)%nt))
  !inv_data(4)%long = -75.93750000000000_dp
  !inv_data(4)%lat = 54.959452452978347_dp
  !inv_data(4)%times(1) = 7000.0_dp
  !inv_data(4)%times(2) = 6500.0_dp
  !inv_data(4)%times(3) = 6000.0_dp
  !inv_data(4)%times(4) = 4000.0_dp
  !inv_data(4)%times(5) = 0.0_dp

  ! 5 = James Bay
  !inv_data(5)%type = 1
  !inv_data(5)%nt = 10
  !allocate(inv_data(5)%times(inv_data(5)%nt),inv_data(5)%values(inv_data(5)%nt))
  !inv_data(5)%long =  -77.343750000000000_dp
  !inv_data(5)%lat =  52.211588148675084_dp
  !inv_data(5)%times(1) = 8000.0_dp
  !inv_data(5)%times(2) = 7500.0_dp
  !inv_data(5)%times(3) = 7000.0_dp
  !inv_data(5)%times(4) = 6500.0_dp
  !inv_data(5)%times(5) = 6000.0_dp
  !inv_data(5)%times(6) = 5500.0_dp
  !inv_data(5)%times(7) = 5000.0_dp
  !inv_data(5)%times(8) = 4500.0_dp
  !inv_data(5)%times(9) = 4000.0_dp
  !inv_data(5)%times(10) = 0.0_dp

  sl_file = 'sl.64.20.p01.'
  call extract_rsl(sl_file,inv_data)

  inv_data(1)%times(:) = inv_data(1)%times(:)*yr2sec/t_norm
  inv_data(1)%values(:) = inv_data(1)%values(:)/r_norm

  inv_data(2)%times(:) = inv_data(2)%times(:)*yr2sec/t_norm
  inv_data(2)%values(:) = inv_data(2)%values(:)/r_norm

  inv_data(3)%times(:) = inv_data(3)%times(:)*yr2sec/t_norm
  inv_data(3)%values(:) = inv_data(3)%values(:)/r_norm

  !inv_data(4)%times(:) = inv_data(4)%times(:)*yr2sec/t_norm
  !inv_data(4)%values(:) = inv_data(4)%values(:)/r_norm

  !inv_data(5)%times(:) = inv_data(5)%times(:)*yr2sec/t_norm
  !inv_data(5)%values(:) = inv_data(5)%values(:)/r_norm

  comb_data%n_point = nloc
  comb_data%n_rsh = 0
  comb_data%pres_point = .true.
  comb_data%pres_rsh = .false.
  allocate(comb_data%point(comb_data%n_point))
  comb_data%point = inv_data

  !===========================================!
  !                 Inversion                 !
  !===========================================!

  file1 = 'inv.check.rsl.'

  call visc_inversion(t2,ice_data,sea_level_data,comb_data,file1)

end program sl_visc_inv
