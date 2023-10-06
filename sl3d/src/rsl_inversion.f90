program rsl_inversion

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
  character(len=256) :: model_file,file1,icefile,file2,sl_file,name

  real(dp) :: drmax,lambda_l

  real(dp), dimension(:,:), allocatable :: sea_t

  complex(dpc), dimension(:,:,:), allocatable :: ice_data
  complex(dpc), dimension(:,:), allocatable :: sea_level_data

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec,ispec1,nt,it1,it,nout,ilayer,inode,ua, &
                  nloc,iloc,ntimes,itime

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,load_lat,load_long,rr,rl,time

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

  !file1 = 'inv.check.rsl.visc.7'
  !call visc_3d_file(file1)
  call visc_3d_none()

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

  open(newunit=io,file='sl_loc2.dat')

  ! Number of locations
  read(io,*) nloc
  allocate(inv_data(nloc))

  open(newunit=io2,file='sl_times_rand.dat')

  do iloc = 1,nloc
     read(io,*) name,lat,long,ntimes
     inv_data(iloc)%type = 1
     inv_data(iloc)%nt = ntimes+1
     inv_data(iloc)%long = long
     inv_data(iloc)%lat = lat
     allocate(inv_data(iloc)%times(inv_data(iloc)%nt))
     allocate(inv_data(iloc)%values(inv_data(iloc)%nt))

     do itime = 1,inv_data(iloc)%nt-1
        read(io2,*) it,time
        if (it /= itime) then
           print *, iloc,itime,it
           stop
        end if
        inv_data(iloc)%times(itime) = time
     end do

     ! Present day
     inv_data(iloc)%times(inv_data(iloc)%nt) = 0.0_dp
  end do
  close(io)

  do iloc = 1,nloc

     inv_data(iloc)%times(:) = inv_data(iloc)%times(:)*yr2sec/t_norm

  end do


  sl_file = 'sl.64.20.p01.full'
  call extract_rsl_full(sl_file,nt,t2,inv_data)

  do iloc = 1,nloc

     inv_data(iloc)%values(:) = inv_data(iloc)%values(:)/r_norm

  end do


  comb_data%n_point = nloc
  comb_data%n_rsh = 0
  comb_data%pres_point = .true.
  comb_data%pres_rsh = .false.
  allocate(comb_data%point(comb_data%n_point))
  comb_data%point = inv_data

  !===========================================!
  !                 Inversion                 !
  !===========================================!

  file1 = 'inv.check.rsl2.'

  call visc_inversion(t2,ice_data,sea_level_data,comb_data,file1)

end program rsl_inversion
