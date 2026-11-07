program visc_from_vs

  !==========================================================================!
  !  Writes the 3D viscosity that visc_3d_from_vs builds from S20RTS, node   !
  !  by node, so that another code's version of the same chain can be held   !
  !  to it.  Set up as visc_test is: the earth model and lmax_run from the   !
  !  command line or a config file, the mesh from lmax.                      !
  !                                                                          !
  !  Output, visc_from_vs.out, one line per node of the viscoelastic region  !
  !  and point of the grid:                                                  !
  !                                                                          !
  !    radius (m)  colatitude (rad)  longitude (rad)  dlnvs  eta / eta_1D    !
  !                                                                          !
  !  where the longitude is the node's geographic one in [-pi, pi), the     !
  !  codes' long = aphi - 180 degrees, dlnvs is the S20RTS value the chain   !
  !  applied at that node (visc_3d_from_vs reads its synthesis across the    !
  !  same half-turn), and eta / eta_1D = si_node / si_spat_node, the inverse !
  !  relaxation time being proportional to 1 / eta.                          !
  !==========================================================================!

  use nrtype
  use module_util
  use module_function
  use module_fourier
  use module_model
  use module_config

  implicit none

  integer(i4b) :: io,ilayer,ispec,inode,igl,iphi,iphih
  real(dp) :: drmax,lambda_l,rr
  real(sp) :: rx
  complex(dpc), dimension(:,:), allocatable :: dlnvs

  call set_parameters

  lmax_run = 64
  call read_config

  open(newunit = io,file=trim(earth_model),action='read',form='formatted')
  call read_model(io)
  close(io)

  lmax = lmax_run
  call set_sh_grid
  call calc_sph_harm

  lambda_l = twopi_d*r(nknot)/(lmax+0.5_dp)
  drmax = lambda_l/2.0_dp
  call mesh_model(drmax)

  call visc_3d_from_vs

  allocate(dlnvs(ngl,nphi))
  open(newunit=io,file='visc_from_vs.out',action='write')
  do ilayer = 1,nsect
     if (sect_ind(3,ilayer) /= 2) cycle
     do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
        do inode = 1,ngll
           rr = r_node(inode,ispec)
           ! as visc_3d_from_vs maps the region onto the splines' [-1, 1]
           rx = (2.0_dp*rr - r_node(1,sect_ind(1,ilayer)) - r_moho) &
                /(r_moho - r_node(1,sect_ind(1,ilayer)))
           call get_dvsvs(rx,dlnvs)
           do iphi = 1,nphi
              if (iphi <= 2*lmax) then
                 iphih = iphi + 2*lmax
              else
                 iphih = iphi - 2*lmax
              end if
              do igl = 1,ngl
                 write(io,'(es24.16,4es25.16)') rr*r_norm, tgl(igl), aphi(iphi) - pi_d, &
                      real(dlnvs(igl,iphih)), &
                      real(si_node(1,inode,ispec))/real(si_spat_node(igl,iphi,1,inode,ispec))
              end do
           end do
        end do
     end do
  end do
  close(io)

end program visc_from_vs
