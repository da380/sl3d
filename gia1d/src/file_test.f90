program file_test

  use nrtype
  use module_util

  implicit none

  integer(i4b) :: io,ilat,ilong,nlat,nlong
  real(dp) :: long,lat,data

  nlat = 256
  nlong = 512

  open(newunit=io,file='p55_ICE5G_ja.21.0',action='read')
  do ilong = 1,nlong
     do ilat = 1,nlat
        read(io,*) long,lat,data
        print *, data
     end do
     read(io,*)
  end do
  close(io)
     


end program file_test
