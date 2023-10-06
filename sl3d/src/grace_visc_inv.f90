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
  t2 = 5000.0_dp*yr2sec/t_norm
  dt = 25.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 201
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = nt

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
