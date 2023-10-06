program kernel_int

  use nrtype
  use module_util
  use module_model
  use module_fourier

  implicit none

  integer(i4b) :: narg,io

  logical(lgt) :: ltmp

  character(len=256) :: model_file,str1

  integer(i4b) :: lmax,ilayer,ispec,inode,igl,iphi,int1,int2
  real(dp) :: drmax,lambda_l,real1,real2

  complex(dpc) :: sint
  real(dp) :: integral,rr,kern,rsi
  complex(dpc), dimension(:,:), allocatable :: spat

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

  call set_sh_grid(lmax)
  call calc_sph_harm(lmax)

  lambda_l = twopi_d*r(nknot)/(lmax+0.5_dp)
  drmax = lambda_l/2.0_dp
  call mesh_model(drmax)

  call global_points(lmax)

  call visc_3d(lmax)

  allocate(spat(ngl,nphi))

  integral = 0.0_dp

  open(newunit=io,file='vkbulge3dnew.21',action='read')
  do ilayer = 1,nsect
     if (sect_ind(3,ilayer) /= 2) cycle

     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
        do inode = 1,ngll
           spat = 0.0_dp
           rr = r_node(inode,ispec)
           rsi = si_node(1,inode,ispec)
           do iphi = 1,nphi
              do igl = 1,ngl
                 read(io,*) real1,int1,int2,kern
                 !spat(igl,iphi) = -kern*(si_spat_node(igl,iphi,1,inode,ispec)-rsi) &
                 !                     / si_spat_node(igl,iphi,1,inode,ispec)
                 spat(igl,iphi) = -kern*(si_spat_node(igl,iphi,1,inode,ispec)-rsi) &
                                      /rsi
              end do
              read(io,*) real1,str1
              read(io,*) str1,int1
           end do
           call zero_coef_from_fun(lmax,spat,sint)
           sint = sint*sqrt(fourpi_d)*rr**2
           integral = integral + wgll(inode)*jac_element(ispec)*real(sint)*r_norm
        end do
     end do
  end do
  close(io)

  print *, integral


end program kernel_int
