program read_rob

  use nrtype
  use module_util
  use module_fourier
  use module_sl

  implicit none

  character(len=256) :: name
  integer(i4b) :: io,igl,iphi,io2,nloc,iloc,num
  real(dp) :: tmp,long,lat,xx,yy,visc,var1,var2,var3,var4

  lmax = 64

  call set_sh_grid
  call calc_sph_harm

  call rob_projection

  open(newunit=io,file='sl.64.20.p01.21.0',action='read')
  open(newunit=io2,file='sl.64.20.p01.xy.21.0',action='write')

  do iphi = 1,nphi
     do igl = 1,ngl
        read(io,*) long,lat,var1,var2
        call get_xy(lat,long,xx,yy)
        write(io2,*) long,lat,xx,yy,var1,var2
     end do
     read(io,*)
     write(io2,*)
  end do
  close(io)
  close(io2)

  !open(newunit=io,file='sl_loc2.dat',action='read')
  !open(newunit=io2,file='sl_loc3.dat',action='write')
  !read(io,*) nloc
  !do iloc = 1,nloc
  !   read(io,*) name,lat,long,num
  !   call get_xy(lat,long,xx,yy)
  !   write(io2,*) name,long,lat,xx,yy,num
  !end do
  !close(io)
  !close(io2)
  


end program read_rob
