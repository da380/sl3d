program read_kl

  use nrtype
  use module_util
  use module_fourier
  use module_gia

  implicit none

  integer(i4b) :: it,nlat,nlong,ilat,ilong,io,ilat1,ilat2,ilong1,ilong2,igl,iphi
  real(dp) :: long,lat,tmp1,tmp2,int1d,int3d,xx,yy,f11,f12,f21,f22
  real(dp), dimension(:,:), allocatable :: sl1d,sl3d

  lmax = 256

  call set_sh_grid

  nlat = 181
  nlong = 361

  allocate(sl1d(nlat,nlong),sl3d(nlat,nlong))
  sl1d = 0.0_dp
  sl3d = 0.0_dp

  open(newunit=io,file='p55_ICE5G_kl.0.0',action='read')
  do ilong = 1,nlong
     do ilat = 1,nlat
        read(io,*) long,lat,tmp1,tmp2
        sl1d(ilat,ilong) = tmp1
        sl3d(ilat,ilong) = tmp2
     end do
     read(io,*)
  end do
  close(io)

  open(newunit=io,file='p55_ICE5G_kl_ed.0.0')
  do iphi = 1,nphi
     long = aphi(iphi)*rad2deg - 180.0_dp
     ilong1 = int(floor(long)) + 181
     ilong2 = ilong1 + 1

     xx = long - floor(long)

     do igl = 1,ngl

        lat = (pio2_d - tgl(igl))*rad2deg
        ilat1 = int(floor(lat)) + 91
        ilat2 = ilat1 + 1

        yy = lat - floor(lat)

        f11 = sl1d(ilat1,ilong1)
        f12 = sl1d(ilat1,ilong2)
        f21 = sl1d(ilat2,ilong1)
        f22 = sl1d(ilat2,ilong2)

        int1d = f11 + xx*(f12 - f11) + yy*(f21 - f11) + xx*yy*(f11 + f22 - f12 - f21)

        f11 = sl3d(ilat1,ilong1)
        f12 = sl3d(ilat1,ilong2)
        f21 = sl3d(ilat2,ilong1)
        f22 = sl3d(ilat2,ilong2)

        int3d = f11 + xx*(f12 - f11) + yy*(f21 - f11) + xx*yy*(f11 + f22 - f12 - f21)

        write(io,*) long,lat,int1d,int3d

     end do
     write(io,*)
  end do
  close(io)
        


end program read_kl
