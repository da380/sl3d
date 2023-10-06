program read_kl_2

  use nrtype
  use module_util
  
  implicit none

  integer(i4b) :: io,igl,iphi,io2,ngl,nphi
  real(dp) :: tmp

  ngl = 256
  nphi = 512

  open(newunit=io,file='n_200',action='read')
  open(newunit=io2,file='n_200.out')
  do igl = 1,ngl
     do iphi = 1,nphi
        read(io,*) tmp
        write(io2,*) iphi,(ngl+1-igl),tmp
     end do
     write(io2,*)
  end do
  close(io)
  close(io2)



end program read_kl_2
