program grace_visc_kern

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

  integer(i4b) :: ilat,ilong,ilonged,igl,iphi,ispec1,nt,it1,it,tmp1,tmp2,load_l,load_m, &
                  ilayer,inode,ispec
  complex(dpc) :: tmp3

  real(dp) :: tmp,long,lat,xx,yy,t1,t2,dt,rr

  complex(dpc), dimension(:,:,:), allocatable :: sl_fin,ice_kern
  complex(dpc), dimension(:,:,:), allocatable :: UVp
  complex(dpc), dimension(:,:,:,:,:,:), allocatable :: mem
  complex(dpc), dimension(:,:), allocatable :: ice_load
  complex(dpc), dimension(:,:,:,:), allocatable :: kern

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
     do inode = 1,ngll
        write(io,*) ispec, r_node(inode,ispec)*r_norm/1000.0_dp,rho_node(inode,ispec)*rho_norm, &
             mu_node(inode,ispec)*con_norm, kappa_node(inode,ispec)*con_norm
     end do
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
  t2 = 6000.0_dp*yr2sec/t_norm
  dt = 50.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 121
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  ! Construct A matrices
  call construct_system_matrix

  allocate(sl_fin(ngl,nphi,nt))
  allocate(UVp((lmax+1)**2,nglob_ssg,nt))
  allocate(mem(5,(lmax+1)**2,ntau,ngll,nspec,nt))
  sl_fin = 0.0_dp
  UVp = 0.0_dp
  mem = 0.0_dp

  ! Forward calculation
  call forward_full_sl(nt,dt,nt,ice_data,sea_level_data,sl_fin,UVp,mem)



  !============================================!
  !            Adjoint calculation             !
  !============================================!

  ! Set load
  allocate(adj_load%load_lm((lmax+1)**2))
  adj_load%type = 3

  load_l = 20
  load_m = 9

  call real_sph_load(load_l,load_m,adj_load%load_lm)

  allocate(kern(ngl,nphi,ngll,nspec))

  call adj_kern_visc_full(nt,dt,ice_data,sl_fin,adj_load,UVp,mem,kern)
  open(newunit = io,file = 'gracevk20964.6')

  do ilayer = 1,nsect
     
     ! No terms in fluid or elastic regions
     if (sect_ind(3,ilayer) /= 2) cycle

     ! begin loop over the spectral elements
     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
                    
        do inode = 1,ngll
           rr = r_node(inode,ispec)*r_norm
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

  stop


end program grace_visc_kern
