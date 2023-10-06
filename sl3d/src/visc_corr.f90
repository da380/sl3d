program visc_corr

  use nrtype
  use module_fourier
  use module_sl

  implicit none

  logical(lgt) :: ltmp
  character(len=256) :: model_file
  integer(i4b) :: io,io2,io3,ilayer,ispec,inode,narg,igl,iphi
  real(dp) :: rr,tmp1,tmp2,tmp3,tmp4,visc,corr,drmax,lambda_l,visct2,visci2,long,lat,xx,yy
  complex(dpc) :: visct_0,visci_0,viscp_0
  complex(dpc), dimension(:,:), allocatable :: visct,visci,viscp
  !real(dp), dimension(:,:), allocatable :: vinti,vintt,vintp

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

  !allocate(vinti(ngl,nphi),vintt(ngl,nphi),vintp(ngl,nphi))
  !vinti = 0.0_dp
  !vintt = 0.0_dp
  !vintp = 0.0_dp

  allocate(visci(ngl,nphi),visct(ngl,nphi),viscp(ngl,nphi))

  !open(newunit=io,file='visc.s20.p01',action='read')
  !open(newunit=io2,file='inv.check.rsl2.visc.10',action='read')
  !do ilayer = 1,nsect
     ! No terms in fluid or elastic regions
  !   if (sect_ind(3,ilayer) /= 2) cycle

  !   do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
  !      do inode = 1,ngll
  !         do iphi = 1,nphi
  !            do igl = 1,ngl
  !               read(io,*) rr,tmp1,tmp2,tmp3,tmp4,visc
  !               if (ispec <= 43) then
  !                  visct2 = log10(visc) - log10(5e21)
  !               else
  !                  visct2 = log10(visc) - log10(5e20)
  !               end if
  !               read(io2,*) rr,tmp1,tmp2,visc
  !               if (ispec <= 43) then
  !                  visci2 = log10(visc) - log10(5e21)
  !               else
  !                  visci2 = log10(visc) - log10(5e20)
  !               end if
  !               vinti(igl,iphi) = vinti(igl,iphi) + (visci2**2)*wgll(inode)*jac_element(ispec)
  !               vintt(igl,iphi) = vintt(igl,iphi) + (visct2**2)*wgll(inode)*jac_element(ispec)
  !               vintp(igl,iphi) = vintp(igl,iphi) + visci2*visct2*wgll(inode)*jac_element(ispec)
  !            end do
  !            read(io,*)
  !            read(io,*)
  !            read(io2,*)
  !            read(io2,*)
  !         end do
  !      end do
  !   end do
  !end do
  !close(io)
  !close(io2)

  !open(newunit=io3,file='visc.inv.corr.surf',action='write')
  !do iphi = 1,nphi
  !   long = aphi(iphi)*rad2deg - 180.0_dp
  !   do igl = 1,ngl
  !      lat = (pio2_d - tgl(igl))*rad2deg
  !      call get_xy(lat,long,xx,yy)
  !      write(io3,*) long,lat,xx,yy,sqrt(vintt(igl,iphi)),sqrt(vinti(igl,iphi)), &
  !                   vintp(igl,iphi),vintp(igl,iphi)/(sqrt(vintt(igl,iphi)*vinti(igl,iphi)))
  !   end do
  !end do
  !close(io3)

  open(newunit=io,file='visc.s20.p01',action='read')
  open(newunit=io2,file='inv.check.rsl2.visc.10',action='read')
  open(newunit=io3,file='visc.inv.corr2',action='write')
  do ilayer = 1,nsect
     ! No terms in fluid or elastic regions
     if (sect_ind(3,ilayer) /= 2) cycle

     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
        do inode = 1,ngll
           do iphi = 1,nphi
              do igl = 1,ngl
                 read(io,*) rr,tmp1,tmp2,tmp3,tmp4,visc
                 if (ispec <= 43) then
                    visct(igl,iphi) = log10(visc) - 21.699
                 else
                    visct(igl,iphi) = log10(visc) - log10(5e20)
                 end if
                 read(io2,*) rr,tmp1,tmp2,visc
                 if (ispec <= 43) then
                    visci(igl,iphi) = log10(visc) - 21.699
                 else
                    visci(igl,iphi) = log10(visc) - log10(5e20)
                 end if
                 viscp(igl,iphi) = visct(igl,iphi)*visci(igl,iphi)
              end do
              read(io,*)
              read(io,*)
              read(io2,*)
              read(io2,*)
           end do
           visct = visct**2
           visci = visci**2
           call zero_coef_from_fun(visct,visct_0)
           call zero_coef_from_fun(visci,visci_0)
           call zero_coef_from_fun(viscp,viscp_0)
           visct_0 = sqrt(real(visct_0))
           visci_0 = sqrt(real(visci_0))
           write(io3,*) rr, real(visct_0),real(visci_0),real(viscp_0),real(viscp_0/(real(visct_0*visci_0)))
        end do
     end do
  end do
  close(io)
  close(io2)
  close(io3)

end program visc_corr
