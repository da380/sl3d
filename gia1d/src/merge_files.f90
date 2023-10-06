program merge_files

  use nrtype
  use module_util
  use module_fourier
  use module_gia

  implicit none

  integer(i4b) :: lmax,iphi,igl,io1,io2,io3,io4,it
  real(dp) :: lat,long,xx,yy,tmp1,tmp2,tmp3,tmp4,tmp5,var1,var2,var3
  character(len=256) :: file1,file2

  lmax = 64

  call set_sh_grid(lmax)
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

  open(newunit=io1,file='sl.64.20.p015.l.0.0',action='read')
  open(newunit=io2,file= 'sl.64.40.0.0.0',action='read')
  open(newunit=io3,file='sl.64.20.p015.lcomp.0.0')

     do iphi = 1,nphi
        do igl = 1,ngl
           read(io1,*) long,lat,tmp1,var1
           read(io2,*) long,lat,tmp2,var2
           write(io3,*) long,lat,var1,var2
        end do
        read(io1,*)
        read(io2,*)
        write(io3,*)
     end do
     close(io1)
     close(io2)
     close(io3)
  !end do




end program merge_files
