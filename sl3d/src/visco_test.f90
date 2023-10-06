program visco_test

  use nrtype
  use module_util
  use module_model
  use module_mat
  use module_fourier
  use module_gia

  implicit none

  integer(i4b) :: narg,io,io1

  logical(lgt) :: ltmp
  character(len=256) :: model_file,file1

  integer(i4b) :: ll,kl,ku,kd,ldab,info,ispec1,wa,inode,ispec,nrows, &
                  lmcur,igl,iphi,nt,nout,it,iout,i
  real(dp) :: drmax,lambda_l,rl,t1,t2,dt

  complex(dpc), dimension(:), allocatable :: load
  complex(dpc), dimension(:,:), allocatable :: load_spat,Us,ps,ten_test,ten_test2
  complex(dpc), dimension(:,:,:), allocatable :: U_spat,p_spat,ten_spat

  complex(dpc), dimension(:,:,:), allocatable :: UVp,W
  complex(dpc), dimension(:,:,:,:,:,:,:), allocatable :: mem
  

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

  lambda_l = twopi_d*r(nknot)/(lmax+0.5_dp)
  drmax = lambda_l/2.0_dp
  call mesh_model(drmax)

  call global_points

  call visc_3d

  !------------------------------------!
  !       Construct A matrices         !
  !------------------------------------!

  call construct_system_matrix

  !---------------------------------------!
  !               Set load                !
  !---------------------------------------!

  allocate(load((lmax+1)**2))
  load = 0.0_dp
  do ll = 1,lmax
     lmcur = ll**2 + ll + 1
     load(lmcur) = exp(-twopi_d*(ll+1)/(10 + 0.5_dp))*sqrt((2*ll + 1)/fourpi_d)
     !load(lmcur) = sqrt((2*ll + 1)/fourpi_d)
  end do

  ! Let's just assume time derivative of load is constant and equal to the above
  ! Make a reasonable size
  load = load*t_norm/(30.0_dp*365.0_dp*24.0_dp*60.0_dp*60.0_dp*r_norm)

  !allocate(load_spat(ngl,nphi))

  !call fun_from_coefs(lmax,load,load_spat)

  !open(newunit=io,file='load.out')
  !do iphi = 1,nphi
  !   do igl = 1,ngl
  !      write(io,*) iphi,igl,real(load_spat(igl,iphi))*r_norm/(t_norm)
  !   end do
  !end do
  !close(io)
  !stop

  !===========================================!
  !             GIA calculation               !
  !===========================================!

  ! set time stepping parameters
  t1 = 0.0_dp
  t2 = 5000.0_dp*yr2sec/t_norm
  dt = 50.0_dp*yr2sec/t_norm
  !nt = (t2-t1)/dt + 1
  nt = 101
  dt = (t2-t1)/(nt-1)
  print *, dt*t_norm*sec2yr

  nout = 6

  allocate(UVp((lmax+1)**2,nglob_ssg,2),W((lmax+1)**2,nglob_tor,2), &
           mem(5,ngl,nphi,ntau,ngll,nspec,2))


  !allocate(ten_test(5,(lmax+1)**2))
  !do ll = 1,lmax
  !   lmcur = ll**2 + ll + 1
  !   if (ll == 1) then
  !      ten_test(2:4,lmcur) = sqrt((2*ll + 1)/fourpi_d)
  !   else
  !      ten_test(1:5,lmcur) = sqrt((2*ll + 1)/fourpi_d)
  !   end if
  !end do
  !allocate(ten_spat(5,ngl,nphi))
  !allocate(ten_test2(5,(lmax+1)**2))
  !call fun_from_coefs_ten_tl(lmax,ten_test,ten_spat)
  !call coefs_from_fun_ten_tl(lmax,ten_spat,ten_test2)

  !call fun_from_coefs(lmax,ten_test(3,:),ten_spat(3,:,:))

  !open(newunit=io,file='tenspat.out')
  !do iphi = 1,nphi
  !   do igl = 1,ngl
  !      write(io,*) iphi,igl,ten_spat(3,igl,iphi)
  !   end do
  !   write(io,*)
  !end do
  !close(io)


  !open(newunit=io,file='tentest.out')
  !do ll = 1 ,lmax
  !   do i = 1,5
  !      write(io,*) ll, i, real(ten_test(i,ll**2+ll+1)), real(ten_test2(i,ll**2+ll+1))
  !   end do
  !   write(io,*)
  !end do
  !close(io)
  !stop



  UVp = 0.0_dp
  W = 0.0_dp
  mem = 0.0_dp

  allocate(ps((lmax+1)**2,nout),Us((lmax+1)**2,nout))
  allocate(p_spat(ngl,nphi,nout),U_spat(ngl,nphi,nout))
  ps = 0.0_dp
  Us = 0.0_dp
  p_spat = 0.0_dp
  U_spat = 0.0_dp

  iout = 1

  do it = 2,nt

     print *, it

     if (it <= 201) then

        call update_full_rk(dt,mem(:,:,:,:,:,:,1),UVp(:,:,1),W(:,:,1),mem(:,:,:,:,:,:,2), &
                            UVp(:,:,2),W(:,:,2),load)
     else

        call update_full_rk(dt,mem(:,:,:,:,:,:,1),UVp(:,:,1),W(:,:,1),mem(:,:,:,:,:,:,2), &
                            UVp(:,:,2),W(:,:,2))
     end if

     if (mod(it,20) == 1) then
        iout = iout + 1

        call surface_fields_full(UVp(:,:,2),ps(:,iout),Us(:,iout))
        call fun_from_coefs(Us(:,iout),U_spat(:,:,iout))
        call fun_from_coefs(ps(:,iout),p_spat(:,:,iout))

     end if

     mem(:,:,:,:,:,:,1) = mem(:,:,:,:,:,:,2)
     UVp(:,:,1) = UVp(:,:,2)
     W(:,:,1) = W(:,:,2)

  end do

  open(newunit=io1,file='visctdispl.out')
  do iout = 1,nout
     write(io1,*) (iout - 1)*1000.0_dp, real(U_spat(64,64,iout))*r_norm, &
                  real(p_spat(64,64,iout))*pot_norm

     call string_cat_int('viscdispl.',iout-1,file1)
     open(newunit=io,file=file1)
     do iphi = 1,nphi
        do igl = 1,ngl
           write(io,*) iphi,igl,real(U_spat(igl,iphi,iout)*r_norm), &
                       real(p_spat(igl,iphi,iout)*pot_norm)
        end do
        write(io,*)
     end do
     close(io) 
  end do
  close(io1)


end program visco_test
