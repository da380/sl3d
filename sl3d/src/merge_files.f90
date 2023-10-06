program merge_files

  use nrtype
  use module_util
  use module_fourier
  use module_gia

  implicit none

  integer(i4b) :: iphi,igl,io1,io2,io3,io4,it
  integer(i4b), dimension(10) :: io
  real(dp) :: lat,long,xx,yy,tmp1,tmp2,tmp3,tmp4,tmp5,var1,var2,var3
  real(dp), dimension(10) :: var
  character(len=256) :: file1,file2

  lmax = 64

  call set_sh_grid
  call rob_projection

  !do it = 1,43
  !   call string_cat_int('sea_level_check25.',21-it/2,file1)
  !   if (mod(it,2) == 0) then
  !      file2 = trim(file1) // trim('.5')
  !   else
  !      file2 = trim(file1) // trim('.0')
  !   end if     
  !   open(newunit=io1,file=file2,action='read')

  !   call string_cat_int('jasl25.',21-it/2,file1)
  !   if (mod(it,2) == 0) then
  !      file2 = trim(file1) // trim('.5')
  !   else
  !      file2 = trim(file1) // trim('.0')
  !   end if
  !   open(newunit=io2,file=file2,action='read')

  !   call string_cat_int('sl_both_25_diff.',21-it/2,file1)
  !   if (mod(it,2) == 0) then
  !      file2 = trim(file1) // trim('.5')
  !   else
  !      file2 = trim(file1) // trim('.0')
  !   end if
  !   open(newunit=io3,file=file2)

  !open(newunit=io1,file='sl.64.pert.200.21.0',action='read')
  !open(newunit=io2,file= 'sl.inv.200.init.10',action='read')
  !open(newunit=io3,file='sl.inv.200.init.comp.10')

  !open(newunit=io1,file='inv.test.visc.1.968',action='read')
  !open(newunit=io2,file= 'inv.test.visc.2.968',action='read')
  !open(newunit=io3,file='inv.test.visc.comp.968')

  !open(newunit=io1,file = 'vktah.1d.18.1756',action='read')
  open(newunit=io2,file = 'vkbulge.s20.p01.interp.9.635',action='read')
  open(newunit=io3,file = 'vkbulge.s20.p01.xy.9.635')

  do iphi = 1,nphi
     do igl = 1,ngl
        !read(io1,*) long,lat,xx,yy,var1
        !read(io2,*) long,lat,xx,yy,var2
        !if ((abs(var1) < 0.002_dp) .or. (abs(var2) < 0.002_dp)) then
        !   var1 = 0.0_dp
        !   var2 = 0.0_dp
        !end if
        !write(io3,*) long,lat,xx,yy,var1,var2
        !read(io1,*) tmp1,long,lat,var2
        read(io2,*) tmp2,long,lat,var1
        call get_xy(lat,long,xx,yy)
        write(io3,*) long,lat,xx,yy,var1!,var2
     end do
     !read(io1,*)
     read(io2,*)
     write(io3,*)
  end do
  !close(io1)
  close(io2)
  close(io3)
  !end do


  stop

  open(newunit=io1,file='inv.test.grace.2.visc.all.2150')

  do it = 1,10
     call string_cat_int('inv.test.grace.2.visc.',it,file1)
     file1 = trim(file1) // trim('.2150')
     open(newunit=io(it),file=file1,action='read')
  end do

  do iphi = 1,nphi
     do igl = 1,ngl
        do it = 1,10
           read(io(it),*) tmp1,long,lat,var1
           var(it) = var1
        end do
        write(io1,*) long,lat,var(:)
     end do
     do it = 1,10
        read(io(it),*)
     end do
     write(io1,*)
  end do

  close(io1)
  do it = 1,10
     close(io(it))
  end do
        




end program merge_files
