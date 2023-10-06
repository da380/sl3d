module module_model

  
  use nrtype
  use module_fourier
  use module_spline
  implicit none


  !===============================================================!
  !           normalization parameters for the model              !
  !===============================================================!
  real(dp), parameter :: r_norm = 6.371e+6_dp                     !
  real(dp), parameter :: rho_norm = 5.515e+3_dp                   ! 
  real(dp), parameter :: big = 1.e+40_dp                          !
  real(dp), save      :: vel_norm                                 !
  real(dp), save      :: acl_norm                                 !
  real(dp), save      :: fre_norm                                 !
  real(dp), save      :: t_norm                                   !
  real(dp), save      :: con_norm                                 !
  real(dp), save      :: grav_norm                                !
  real(dp), save      :: moment_norm                              !
  real(dp), save      :: pot_norm                                 !
  real(dp), save      :: pibigg                                   !
  real(dp), save      :: rot_rate                                 !
  real(dp), save      :: rotfac                                   !
  real(dp), save      :: mass_norm                                !
  real(dp), save      :: visco_norm                               !
  real(dp), save      :: tmin                                     !
  real(dp), save      :: rice                                     !
  real(dp), save      :: roce                                     !
  real(dp), save      :: rcru                                     !
  !---------------------------------------------------------------!



  !===============================================================!
  !            arrays holding the elastic model properties        !
  !===============================================================!
  ! The parameters are:                                           !
  ! nlayer == the number of distinct layers in the model, that    !
  !           is, the number of regions in the model separated by !
  !           parameter discontinities.                           !
  ! nknot  == the total number of radial knots in the model.      !
  ! layer_index == an array containing the botom and top indices  !
  !                of a given layer in the model. The dimensions  !
  !                of the array are nlayer-by-3, and the storage  !
  !                format is such that:                           !
  !                layer_index(n,1) == bottom index of nth layer, !
  !                layer_index(n,2) == top index of nth layer.    !
  ! ntau   == number of Maxwell elements in the time-dependent    !
  !           shear modulus                                       !
  ! r      == array of dimension nknot containing the position of !
  !           each radial knot in the model.                      !
  ! rho    == array of dimension nknot containing the density at  !
  !           each radial knot in the model.                      !
  ! rho_cs == array of dimension knot containing the cubic spline !
  !           coefficient of density at each radial knot in the   !
  !           model.                                              !
  ! kappa   == array of dimension nknot containing the bulk       !
  !            modulus at each radial knot in the model.          !
  ! kappa_cs == array of dimension knot containing the cubic      !
  !             spline coefficient of the bulk modulus at each    !
  !             radial knot in the  model.                        !
  ! mu       == array of dimension nknot containing the shear     !
  !             modulus at each radial knot in the model.         !
  ! mu_cs    == array of dimension knot containing the cubic      !
  !             spline coefficient of the shear modulus at each   !
  !             radial knot in the  model.                        !
  !---------------------------------------------------------------!
  logical(lgt), save                          :: ocean            !
  integer(i4b), save                          :: nlayer           !
  integer(i4b), save                          :: nknot            !
  integer(i4b), dimension(:,:), &                                 !
       allocatable , save                     :: layer_index      !
  real(dp), save                              :: mass             !
  real(dp), save                              :: moment           !
  integer(i4b), save                          :: ntau             !
  real(dp), save                              :: tau_min          !
  real(dp), save                              :: tau_max          !
  real(dp), save                              :: tau1             !
  real(dp), dimension(:), allocatable, save   :: r                !
  real(dp), dimension(:), allocatable, save   :: rho              !
  real(dp), dimension(:), allocatable, save   :: rho_cs           !
  real(dp), dimension(:), allocatable, save   :: kappa            !
  real(dp), dimension(:), allocatable, save   :: kappa_cs         !
  real(dp), dimension(:), allocatable, save   :: mu               !
  real(dp), dimension(:,:), allocatable, save :: mui              !
  real(dp), dimension(:), allocatable, save   :: mu_cs            !
  real(dp), dimension(:,:), allocatable, save :: mui_cs           !
  real(dp), dimension(:), allocatable, save   :: visco            !
  real(dp), dimension(:,:), allocatable, save :: viscoi           !
  real(dp), dimension(:), allocatable, save   :: visco_cs         !
  real(dp), dimension(:,:), allocatable, save :: viscoi_cs        !
  !---------------------------------------------------------------!


  !===============================================================!
  ! merged radius array                                           !
  !===============================================================!
  real(dp), dimension(:), allocatable :: rm


  !===============================================================!
  ! number of nodes points per wavelength in the mesh             !
  !===============================================================!
  integer(i4b), parameter :: npw = 12                             !
  !---------------------------------------------------------------!

  !---------------------------------------------------------------!
  ! layer_spec == n by 2 dimensional array, where n is the number !
  !               of layers separated by discontinuities. Stores  !
  !               spectral element at TOP of layer and whether    !
  !               layer is solid (1) or fluid (2).                !                        
  ! ispec_ind does this
  !---------------------------------------------------------------!



  !===============================================================!
  ! Mesh parameters and variables:                                !
  !===============================================================!
  integer(i4b), parameter :: ngll = 5                             !
  integer(i4b), save :: nspec                                     !
  integer(i4b), save :: nglob                                     !
  integer(i4b), save :: nglob_ssg                                 !
  integer(i4b), save :: nglob_tor                                 !
  integer(i4b), save :: nglob_ssg_0                               !
  integer(i4b), save :: ispec_nic                                 !
  integer(i4b), save :: ispec_noc                                 !
  integer(i4b), save :: ispec_ddp                                 !
  integer(i4b), save :: ispec_670                                 !
  integer(i4b), save :: ispec_moho                                !
  integer(i4b), save :: ispec_nsl                                 !
  integer(i4b), save :: ispec_lith                                !
  integer(i4b), save :: nspec_mantle                              !
  integer(i4b), save :: nsect                                     !
  integer(i4b), dimension(:), allocatable, save :: pa_surf        !
  integer(i4b), dimension(:,:), allocatable, save :: ibool        !
  integer(i4b), dimension(:,:), allocatable, save :: ispec_ind    !
  integer(i4b), dimension(:,:), allocatable, save :: sect_ind     !
  real(dp), dimension(ngll), save :: xigll                        !
  real(dp), dimension(ngll), save :: wgll                         !
  real(dp), dimension(:), allocatable, save :: jac_element        !
  real(dp), dimension(:), allocatable, save :: ijac_element       !
  real(dp), dimension(ngll,ngll), save  :: hprime                 !
  real(dp), dimension(:,:), allocatable, save :: r_node           !
  real(dp), dimension(:,:), allocatable, save :: rho_node         !
  real(dp), dimension(:,:), allocatable, save :: drho_node        !
  real(dp), dimension(:,:), allocatable, save :: mu_node          !
  real(dp), dimension(:,:,:), allocatable, save :: mui_node       !
  real(dp), dimension(:,:), allocatable, save :: kappa_node       !   
  real(dp), dimension(:,:), allocatable, save :: grav_node        !
  real(dp), dimension(:,:), allocatable, save :: ellip_node       !
  real(dp), dimension(:,:), allocatable, save :: pres_node        !
  real(dp), dimension(:,:), allocatable, save :: phi_node         !
  real(dp), dimension(:,:), allocatable, save :: rho20_node       !
  real(dp), dimension(:,:), allocatable, save :: phi20_node       !
  real(dp), dimension(:,:), allocatable, save :: pres20_node      !
  real(dp), dimension(:,:), allocatable, save :: visco_node       !
  real(dp), dimension(:,:,:), allocatable, save :: viscoi_node    !
  real(dp), dimension(:,:,:), allocatable, save :: si_node        !
  !---------------------------------------------------------------!

  
contains

  


  
  subroutine read_model(io,isostatic_in)                                     
    !===========================================================!
    ! given a unit number to an open model file, this           !
    ! routine reads in the elastic model file and then puts the !
    ! model parameters into the required form for the           !
    ! program.                                                  !
    !===========================================================!
    use nrtype
    use module_spline
    implicit none
    

    ! Inputs: 
    integer(i4b), intent(in) :: io
    logical(lgt), intent(in), optional :: isostatic_in
    

    ! local variables: 
    character(len=10) :: ctmp
    integer(i4b) :: i,j,itau
    logical(lgt) :: isostatic

    if(present(isostatic_in)) then
       isostatic = isostatic_in
    else
       isostatic = .false.
    end if
    

    
    ! read in the header to the model file
    read(io,*) ctmp      
    read(io,*) nknot,ntau
 

    ! allocate the arrays needed for the model
    call allocate_model(nknot)
    
    ! read in the model file
    do i=1,nknot
       !read(io,*) r(i),rho(i),kappa(i),mu(i),mui(:,i),visco(i),viscoi(:,i)
       !mui(:,i) = mui(:,i)*mu(i)
       !viscoi(:,i) = viscoi(:,i)*visco(i)
       read(io,*) r(i), rho(i), kappa(i), mu(i), visco(i)
       mui(1,i) = mu(i)
       viscoi(1,i) = visco(i)
    end do

    where(visco /= 0.0_dp) visco = 10.0_dp**visco
    where(viscoi /= 0.0_dp) viscoi = 10.0_dp**viscoi

    r=r/r_norm
    rho = rho/rho_norm
    kappa = kappa*1e9_dp/con_norm
    mu = mu*1e9_dp/con_norm  
    mui = mui*1e9_dp/con_norm 
    visco = visco/visco_norm
    viscoi = viscoi/visco_norm


    ! determine the number of layers in the model
    nlayer=1
    do i=2,nknot
       if(r(i) == r(i-1)) nlayer=nlayer+1
    end do
    allocate(layer_index(nlayer,2))
    layer_index(1,1)=1
    layer_index(nlayer,2)=nknot
    j=0
    do i=2,nknot
       if(r(i) == r(i-1)) then
          j=j+1
          layer_index(j,2)=i-1
          layer_index(j+1,1)=i
       end if
    end do
    
    print *, nknot,ntau

    ! compute some of the cubic spline coefficients 
    do i=1,nlayer
       call spline(r(layer_index(i,1):layer_index(i,2)), & 
            rho(layer_index(i,1):layer_index(i,2)), & 
            big,big,rho_cs(layer_index(i,1):layer_index(i,2)))
       call spline(r(layer_index(i,1):layer_index(i,2)), & 
            kappa(layer_index(i,1):layer_index(i,2)), & 
            big,big,kappa_cs(layer_index(i,1):layer_index(i,2)))
       call spline(r(layer_index(i,1):layer_index(i,2)), & 
            mu(layer_index(i,1):layer_index(i,2)), & 
            big,big,mu_cs(layer_index(i,1):layer_index(i,2)))
       call spline(r(layer_index(i,1):layer_index(i,2)), & 
            visco(layer_index(i,1):layer_index(i,2)), & 
            big,big,visco_cs(layer_index(i,1):layer_index(i,2)))
       do itau = 1,ntau
          call spline(r(layer_index(i,1):layer_index(i,2)), & 
               mui(itau,layer_index(i,1):layer_index(i,2)), & 
               big,big,mui_cs(itau,layer_index(i,1):layer_index(i,2)))
          call spline(r(layer_index(i,1):layer_index(i,2)), & 
               viscoi(itau,layer_index(i,1):layer_index(i,2)), & 
               big,big,viscoi_cs(itau,layer_index(i,1):layer_index(i,2)))
       end do
    end do




    return
  end subroutine read_model





  


  subroutine mesh_model(drmax)


    use nrtype
    use module_spline
    implicit none
    
    
    real(dp), intent(in) :: drmax
    
    

    integer(i4b) :: i,j,k,ispec,inode,j1,j2,jnode,itau,ilayer,isect
    integer(i4b), dimension(:), allocatable :: nspec_tmp
    
    real(dp) :: r1,r2,r11,r22,dr, & 
                rmu,rkappa,rvisco,rrho,tau
    real(dp), dimension(ntau) :: rmui,rviscoi
    
    ! external routines used:
    real(dp), external :: lagrange_deriv_GLL
    
    
    
    
    !=======================================================!
    !     get the GLL points and Lagrange derivatives       !
    !=======================================================!
    
    ! get the GLL points and weights
    call zwgljd(xigll,wgll,ngll,0.0_dp,0.0_dp)
    if(mod(ngll,2) /= 0) xigll((ngll-1)/2+1) = 0.0_dp

    
    ! get the derivatives of the Lagrange polynomials at 
    ! the GLL points; recall that  hprime(i,j)=h'_{j}(xigll_{i}) 
    do j=1,ngll
       do i=1,ngll
          hprime(i,j) = lagrange_deriv_GLL(j-1,i-1,xigll,ngll)
       end do
    end do




    !=========================================================!
    ! work out the number of spectral elements that should be !
    ! in the mesh                                             !
    !=========================================================!
    
    allocate(nspec_tmp(nlayer))     

    ! loop over each layer in the model determining how
    ! many spectral elements it should contain
    nspec = 0
    nspec_tmp = 0

    do i = 1,nlayer
       j1 = layer_index(i,1)
       j2 = layer_index(i,2)
       r1 = r(j1)
       r2 = r(j2)
       nspec_tmp(i) = max(1,floor(npw*abs(r2-r1)/(drmax*ngll)))
       if (r1 > 0.894) nspec_tmp(i) = nspec_tmp(i)*5.0_dp
       nspec = nspec + nspec_tmp(i)
    end do

    ! allocate the mesh arrays
    if(allocated(r_node)) then
       call allocate_mesh(ngll,nspec,1)
       call allocate_mesh(ngll,nspec)
    else
       call allocate_mesh(ngll,nspec)
    end if

    !========================================================!
    ! build up the mesh                                      !
    !========================================================!   

    ! allocate array to store the element numbers of 
    ! the elements at the top and bottom of each of the layers
    ! and whether layer is solid or fluid
    allocate(ispec_ind(3,nlayer))

    tau_max = 0.0_dp
    tau_min = 10.0_dp**(20.0_dp)
    
    ispec = 0
    do i = 1,nlayer
       j1  = layer_index(i,1)
       j2  = layer_index(i,2)
       r11 = r(j1)
       r22 = r(j2)
       dr  = (r22-r11)/nspec_tmp(i) 
       do k = 1,nspec_tmp(i)
          ispec = ispec+1
          r1 = r11+(k-1)*dr
          r2 = r1+dr
          jac_element(ispec)  = 0.5_dp*(r2-r1)
          ijac_element(ispec) = 1.0_dp/jac_element(ispec) 
          
          if(k == 1) then
             ispec_ind(1,i) = ispec
          end if
          if(k == nspec_tmp(i)) then
             ispec_ind(2,i) = ispec
          end if
          
          do inode = 1,ngll
             
             r_node(inode,ispec) = r1+0.5_dp*(r2-r1) & 
                  *(xigll(inode)+1.0_dp)
             j = find_radius_index(r_node(inode,ispec))
             do 
                if(j >= j1) exit
                j = j+1
             end do
             do 
                if(j <= j2) exit
                j = j-1
             end do
             
             ! get reference elastic constants
             rrho       = splint_dis(r,rho,rho_cs,r_node(inode,ispec),j)
             rmu        = splint_dis(r,mu,mu_cs,r_node(inode,ispec),j)
             rkappa     = splint_dis(r,kappa,kappa_cs,r_node(inode,ispec),j)
             rvisco     = splint_dis(r,visco,visco_cs,r_node(inode,ispec),j)
             do itau = 1,ntau
                rmui(itau)    = splint_dis(r,mui(itau,:),mui_cs(itau,:),r_node(inode,ispec),j)
                rviscoi(itau) = splint_dis(r,viscoi(itau,:),viscoi_cs(itau,:), &
                                           r_node(inode,ispec),j)
             end do

             ! figure out whether layer is elastic, viscoelastic or fluid
             ! 1 = elastic
             ! 2 = viscoelastic
             ! 3 = fluid
             if (k == 1) then
                if (rmu == 0.0_dp) then
                   ! fluid
                   ispec_ind(3,i) = 3
                else
                   if (rvisco == 0.0_dp) then
                      ! elastic
                      ispec_ind(3,i) = 1
                   else
                      ! viscoelastic
                      ispec_ind(3,i) = 2
                   end if
                end if
             end if


             ! store the nodal values of the model parameters
             rho_node(inode,ispec)      = rrho
             mu_node(inode,ispec)       = rmu
             mui_node(:,inode,ispec)    = rmui
             kappa_node(inode,ispec)    = rkappa
             visco_node(inode,ispec)    = rvisco    
             viscoi_node(:,inode,ispec) = rviscoi
             do itau = 1,ntau
                if (viscoi_node(itau,inode,ispec) == 0.0_dp) then 
                   si_node(itau,inode,ispec) = 0.0_dp
                else
                   si_node(itau,inode,ispec) = mui_node(itau,inode,ispec) &
                                          /viscoi_node(itau,inode,ispec)
                end if
             end do
          end do
       end do
    end do

    ! array similar to ispec_ind, BUT only counts as a different layer if 
    ! there is a change from elastic to viscoelastic etc.
    
    ! calculate number of sections
    nsect = 1
    do ilayer = 2,nlayer
       if (ispec_ind(3,ilayer) /= ispec_ind(3,ilayer-1)) nsect = nsect + 1
    end do
    
    allocate(sect_ind(3,nsect))
    sect_ind(:,1) = ispec_ind(:,1)
    sect_ind(2,nsect) = nspec
    isect = 1
    do ilayer = 2,nlayer
       if (ispec_ind(3,ilayer) /= ispec_ind(3,ilayer-1)) then
          isect = isect + 1
          sect_ind(1,isect) = ispec_ind(1,ilayer)
          sect_ind(3,isect) = ispec_ind(3,ilayer)
          sect_ind(2,isect-1) = ispec_ind(2,ilayer-1)
       end if
    end do


    !====================================================!
    ! calculate the connectivity array for the mesh      !
    !====================================================!
    
    k = 1
    do ispec = 1,nspec
       do inode = 1,ngll
          if(inode /= 1) k = k+1
          ibool(inode,ispec) = k
       end do
    end do
    
    ! number of global nodes
    nglob = maxval(ibool)
    
    
    ! calculate gravity and ellipicity
    call grav_cal


    ! calculate the derivative of the density
    do ispec = 1,nspec
       do inode = 1,ngll
          drho_node(inode,ispec) = 0.0_dp
          do jnode = 1,ngll
             drho_node(inode,ispec) = drho_node(inode,ispec)      & 
                  +rho_node(jnode,ispec) &
                  *hprime(inode,jnode)   & 
                  *ijac_element(ispec)   
          end do
       end do
    end do

    ! find the shortest and longest Maxwell times for the model
    do ilayer = 1,nsect
       if (sect_ind(3,ilayer) /= 2) cycle
       do ispec = sect_ind(1,ilayer),sect_ind(2,ilayer)
          do inode = 1,ngll
             do itau = 1,ntau
                if (si_node(itau,inode,ispec) /= 0.0_dp) then
                   tau = 1.0_dp/si_node(itau,inode,ispec)
                   if(ispec == ispec_noc+1 .and. inode == 1) then
                      tau_min = tau
                      tau_max = tau
                   else if(tau < tau_min) then
                      tau_min = tau
                   else if (tau > tau_max) then
                      tau_max = tau
                   end if
                end if
             end do
          end do
       end do
    end do
 

    if(tau_min == 0.0_dp) tau_min = 100000.0_dp*yr2sec/t_norm    

    return
  end subroutine mesh_model


  subroutine grav_cal
    !==========================================================!
    ! This subroutine computes the, mass, moment of inertia,   !
    ! pressure,  gravity and ellipticity in the earth model    !
    !==========================================================!     
    
    ! modules used
    use nrtype
    use nrutil

    
    implicit none
    
    integer(i4b) :: ispec,iglob,inode,jnode,info,ldab
    integer(i4b), dimension(:), allocatable :: ipivot
    
    real(dp) :: psi
    real(dp), dimension(ngll) :: fin_tmp,fout_tmp
    real(dp), dimension(:,:), allocatable :: f,drho,amat
    
    
    
    
    ! compute the mass of the earth model
    mass = 0.0_dp
    do ispec = 1,nspec
       do inode = 1,ngll
          mass = mass + rho_node(inode,ispec)    & 
                      *wgll(inode)               & 
                      *r_node(inode,ispec)**2    & 
                      *jac_element(ispec)
       end do
    end do
    mass = fourpi_d*mass
    

    ! compute the moment of inertia of the earth model
    moment  = 0.0_dp
    do ispec = 1,nspec
       do inode = 1,ngll
          moment = moment + rho_node(inode,ispec)    & 
                            *wgll(inode)             &
                            *r_node(inode,ispec)**4  &
                            *jac_element(ispec)
       end do
    end do
    moment = (8.0_dp*pi_d/3.0_dp)*moment
    
    ! allocate arrays for the linear systems
    ldab = 3*(ngll-1)+1
    allocate(f(nglob,1))
    allocate(amat(ldab,nglob))
    allocate(ipivot(nglob))

    ! calculate gravitational potential in the reference model
    call force_potr(f(:,1))
    call amat_potr(amat)
    
    ! perform an LU-decomposition on A
    call dgbtrf(nglob,nglob,ngll-1,ngll-1,amat, & 
         ldab,ipivot,info)
    
    ! perform the backsubstitution
    call dgbtrs('N', nglob, ngll-1, ngll-1, 1, amat, &
         ldab,ipivot,f,nglob,info)
    
    
    ! store the reference gravitational potential
    ! and gravitational acceleration            
    
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)
          phi_node(inode,ispec) = f(iglob,1)
       end do
       do inode = 1,ngll
          grav_node(inode,ispec) = 0.0_dp
          do jnode = 1,ngll
             grav_node(inode,ispec) = grav_node(inode,ispec) & 
                                     +phi_node(jnode,ispec)  & 
                                     *hprime(inode,jnode)    & 
                                     *ijac_element(ispec)
          end do
       end do
    end do



    ! compute pressure in the reference earth model
    
    
    do ispec = 1,nspec
       ! integrand for the main pressure term
       fin_tmp(1:ngll) = rho_node(1:ngll,ispec)*grav_node(1:ngll,ispec)
       ! add the contribution from the spherical part of the rotational
       ! potential - should probably shift this to the hetrogeneous 
       ! model section
       fin_tmp(1:ngll) = fin_tmp(1:ngll) -(2.0_dp/3.0_dp)*rot_rate**2 & 
                                          *r_node(1:ngll,ispec)
       call integrate_between_nodes(r_node(1:ngll,ispec),fin_tmp,fout_tmp)
       
       if(ispec == 1) then
          pres_node(1:ngll,ispec) = fout_tmp(1:ngll)
       else
          pres_node(1:ngll,ispec)  =  pres_node(ngll,ispec-1) & 
                                     +fout_tmp(1:ngll)
       end if
       
    end do
    
    pres_node = pres_node(ngll,nspec) - pres_node
      

    
    

    ! compute the radial derivative of the density
    allocate(drho(ngll,nspec))
    
    do ispec = 1,nspec
       do inode = 1,ngll
          drho(inode,ispec) = 0.0_dp
          do jnode = 1,ngll
             drho(inode,ispec) = drho(inode,ispec)      & 
                  +rho_node(jnode,ispec) &
                  *hprime(inode,jnode)   & 
                  *ijac_element(ispec)   
          end do
       end do
    end do
    
    ! assemble the equations for the ellipticity calculation
    call amat_ellip(amat,drho)
    call force_ellip(f(:,1),drho)
    
    ! perform an LU-decomposition on A
    call dgbtrf(nglob,nglob,ngll-1,ngll-1,amat, & 
         ldab,ipivot,info)
    
    ! perform the backsubstitution
    call dgbtrs('N', nglob, ngll-1, ngll-1, 1, amat, &
         ldab,ipivot,f,nglob,info)
    
    ! extract the various parameters
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)
          psi =  rotfac*r_node(inode,ispec)**2
          phi20_node(inode,ispec) = f(iglob,1)+psi
          if(inode /= 1 .or. ispec /= 1) then
             ellip_node(inode,ispec) = (phi20_node(inode,ispec)) & 
                                       /(r_node(inode,ispec)        & 
                                       *grav_node(inode,ispec))
             rho20_node(inode,ispec)  = drho(inode,ispec)       & 
                                       *(phi20_node(inode,ispec)) & 
                                       /grav_node(inode,ispec)
          end if

          pres20_node(inode,ispec) = -rho_node(inode,ispec) & 
                                    *(phi20_node(inode,ispec))

           

       end do
    end do
    
    ellip_node(1,1) = ellip_node(2,1)
    rho20_node(1,1) = rho20_node(2,1)

    ! include the scaling factor for the ellipticity
    !  E_DT = (3/4)sqrt(5/pi)E_ME
    ellip_node = 1.5_dp*sqrt(5.0_dp/fourpi_d)*ellip_node
      


    return
  end subroutine grav_cal



  subroutine amat_potr(amat)
    ! returns system matrix for the reference gravitational potential
    use nrtype
    implicit none
    
    real(dp), dimension(:,:), intent(out) :: amat
    
    integer(i4b) :: inode,jnode,knode, &
         iglob,jglob,ispec,ku,kl,kglob,j1,j2
    

    ! lower and upper bandwidths  for the matrix
    kl = ngll-1
    ku = ngll-1
    
    amat = 0.0_dp

      

    ! calculation of the matrix
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)  
          do jnode = 1,ngll
             jglob = ibool(jnode,ispec)
             ! only do calculations for the lower triangle
             if(jnode > inode) cycle
             kglob = kl+ku+1+iglob-jglob
             do knode = 1,ngll
                amat(kglob,jglob) = amat(kglob,jglob) & 
                     +hprime(knode,inode)    &
                     *hprime(knode,jnode)    &
                     *r_node(knode,ispec)**2 & 
                     *wgll(knode)*ijac_element(ispec)
             end do
          end do
       end do
    end do
    
    ! add in the surface boundary term
    amat(kl+ku+1,nglob) = amat(kl+ku+1,nglob)  & 
         + r_node(ngll,nspec)
    


    ! fill in the upper triangle using the  
    ! symmetry of A
    do iglob = 1,nglob-1
       j1 = iglob+1
       j2 = min(nglob,iglob+ku)
       do jglob = j1,j2
          amat(kl+ku+1+iglob-jglob,jglob)   & 
               = amat(kl+ku+1+jglob-iglob,iglob)
       end do
    end do
    
    return
  end subroutine amat_potr


  
  subroutine amat_ellip(amat,drho)
    ! returns system matrix for the ellipticity 
    ! potential calculations
    use nrtype
    implicit none
    

    real(dp), dimension(:,:), intent(out) :: amat
    real(dp), dimension(:,:), intent(in) :: drho
    
    integer(i4b) :: i,inode,jnode,knode, &
         iglob,jglob,ispec,ku,kl,kglob,j1,j2,l
    real(dp) :: a1,a2
    
    
    l = 2
    
    ! lower and upper bandwidths  for the matrix
    kl = ngll-1
    ku = ngll-1
    
    amat = 0.0_dp
    
      
    !==============================!
    !  calculation of the matrix   !
    !==============================!
    
    !-------------------------!
    ! add the integral terms  !
    !-------------------------!
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)  
          kglob = kl+ku+1
          if(ispec == 1 .and. inode == 1) then
             a1 = drho(inode+1,ispec)/grav_node(inode+1,ispec)
          else
             a1 = drho(inode,ispec)/grav_node(inode,ispec)
          end if
          amat(kglob,iglob) = amat(kglob,iglob) + l*(l+1) & 
               *wgll(inode)*jac_element(ispec)
          amat(kglob,iglob) = amat(kglob,iglob) + 4.0_dp*pibigg*a1 & 
               *r_node(inode,ispec)**2 & 
               *wgll(inode)*jac_element(ispec)
          do jnode = 1,ngll
             jglob = ibool(jnode,ispec)
             if(jnode > inode) cycle
             kglob = kl+ku+1+iglob-jglob
             do knode = 1,ngll
                amat(kglob,jglob) = amat(kglob,jglob) & 
                     +hprime(knode,inode)    &
                     *hprime(knode,jnode)    &
                     *r_node(knode,ispec)**2 & 
                     *wgll(knode)*ijac_element(ispec)
             end do
          end do
       end do
    end do
    
    !-------------------------!
    ! add the boundary terms  !
    !-------------------------!
    
    ! add in the surface boundary term
    amat(kl+ku+1,nglob) = amat(kl+ku+1,nglob)  & 
         + (l+1)*r_node(ngll,nspec)
    
    
    ! add in the discontinuity terms
    do i = 1,nsect
         
       if(i < nsect) then
          ispec = sect_ind(2,i)
          iglob = ibool(ngll,ispec)
          a2 = (rho_node(1,ispec+1)-rho_node(ngll,ispec)) & 
               /grav_node(ngll,ispec)
          amat(kl+ku+1,iglob) = amat(kl+ku+1,iglob) & 
               +4.0_dp*pibigg*a2    & 
               *r_node(ngll,ispec)**2
       else
          a2 = -rho_node(ngll,nspec)/grav_node(ngll,nspec)
          amat(kl+ku+1,nglob) = amat(kl+ku+1,nglob) & 
               +4.0_dp*pibigg*a2    & 
               *r_node(ngll,nspec)**2
          
       end if
       
    end do


    ! fill in the upper triangle using the  
    ! symmetry of A
    do iglob = 1,nglob-1
       j1 = iglob+1
       j2 = min(nglob,iglob+ku)
       do jglob = j1,j2
          amat(kl+ku+1+iglob-jglob,jglob)   & 
               = amat(kl+ku+1+jglob-iglob,iglob)
       end do
    end do
    
    return
  end subroutine amat_ellip



   


  subroutine force_potr(f)
    ! returns force vector for the reference gravitational 
    ! potential calculation
    use nrtype
    implicit none
    
    real(dp), dimension(:), intent(out) :: f
    
    integer(i4b) :: inode,iglob,ispec
    
    f(1:nglob) = 0.0_dp
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)
          f(iglob) = f(iglob) + rho_node(inode,ispec) & 
               *r_node(inode,ispec)*r_node(inode,ispec) & 
               *wgll(inode)*jac_element(ispec)
       end do
    end do
    
    ! recall the normalization is such that
    ! pi*G = 1
    f = -4.0_dp*f
    
    return
  end subroutine force_potr






  subroutine force_ellip(f,drho)
    ! returns the force vector for the ellipticity
    ! potential calculations
    use nrtype
    implicit none
    
    real(dp), dimension(:), intent(out) :: f
    real(dp), dimension(:,:), intent(in) :: drho

    integer(i4b) :: inode,iglob,ispec,ilayer
    
    real(dp) :: a1,a2,psi
    
    f(1:nglob) = 0.0_dp
    
    
    ! add the integral term
    do ispec = 1,nspec
       do inode = 1,ngll
          iglob = ibool(inode,ispec)
          if(iglob == 1) then
             a1 = 0.0_dp
          else
             a1 = drho(inode,ispec)/grav_node(inode,ispec)
          end if
          psi = rotfac*r_node(inode,ispec)**2
          f(iglob) = f(iglob) - 4.0_dp*pibigg*a1*psi     & 
               *r_node(inode,ispec)**2             & 
               *wgll(inode)*jac_element(ispec)
       end do
    end do
    
    
    ! add in the boundary terms
    
    do ilayer = 1,nsect
       
       if(ilayer < nsect) then
          ispec = sect_ind(2,ilayer)
          iglob = ibool(ngll,ispec)
          a2 = (rho_node(1,ispec+1)-rho_node(ngll,ispec)) & 
               /grav_node(ngll,ispec)
          psi = rotfac*r_node(ngll,ispec)**2
          f(iglob) = f(iglob) -4.0_dp*pibigg*a2*psi & 
               *r_node(ngll,ispec)**2
       else
          a2 = -rho_node(ngll,nspec)/grav_node(ngll,nspec)
          psi = rotfac*r_node(ngll,nspec)**2
          f(nglob) = f(nglob) -4.0_dp*pibigg*a2*psi & 
               *r_node(ngll,nspec)**2
       end if
    end do

    return
  end subroutine force_ellip





  subroutine allocate_model(nknot_in,isign)
    !===========================================================!
    ! given a value for nknot this routine allocates            !
    ! the various arrays of the model. If the optional          !
    ! variable isign is set equal to 1 then the routine         !
    ! deallocates these same arrays.                            !
    !===========================================================!
    use nrtype
    implicit none
    integer(i4b), intent(in) :: nknot_in
    integer(i4b), intent(in), optional :: isign
    logical(lgt) :: ltmp
    ltmp=present(isign)
    if(ltmp) then
       if(isign /= 1) ltmp=.false.
    end if
    if(.not.ltmp) then
       allocate(r(nknot_in),        & 
            rho(nknot_in),          & 
            rho_cs(nknot_in),       & 
            kappa(nknot_in),        &   
            kappa_cs(nknot_in),     & 
            mu(nknot_in),           & 
            mui(ntau,nknot_in),     & 
            mu_cs(nknot_in),        & 
            mui_cs(ntau,nknot_in),  & 
            visco(nknot_in),        & 
            viscoi(ntau,nknot_in),  & 
            visco_cs(nknot_in),     & 
            viscoi_cs(ntau,nknot_in)) 
            
    else
       deallocate(r,       & 
            rho,           & 
            rho_cs,        &
            kappa,         & 
            kappa_cs,      & 
            mu,            &
            mui,           &
            mu_cs,         &
            mui_cs,        &
            visco,         &
            viscoi,        &
            visco_cs,      &
            viscoi_cs,     & 
            layer_index)
    end if
    return
  end subroutine allocate_model





  subroutine allocate_mesh(ngll,nspec,isign)
    use nrtype
    implicit none
    integer(i4b), intent(in) :: ngll,nspec
    integer(i4b), intent(in), optional :: isign
    logical(lgt) :: ltmp
    ltmp=present(isign)
    if(ltmp) then
       if(isign /= 1) ltmp=.false.
    end if
    if(.not.ltmp) then
       allocate(jac_element(nspec),   &
            ijac_element(nspec),      &
            r_node(ngll,nspec),       &
            rho_node(ngll,nspec),     &
            drho_node(ngll,nspec),    &
            grav_node(ngll,nspec),    &
            ellip_node(ngll,nspec),   & 
            pres_node(ngll,nspec),    &
            phi_node(ngll,nspec),     & 
            rho20_node(ngll,nspec),   & 
            phi20_node(ngll,nspec),   & 
            pres20_node(ngll,nspec),  & 
            mu_node(ngll,nspec),      & 
            mui_node(ntau,ngll,nspec),&
            si_node(ntau,ngll,nspec), &
            kappa_node(ngll,nspec),   &
            visco_node(ngll,nspec),      &
            viscoi_node(ntau,ngll,nspec),&
            ibool(ngll,nspec))
    else
       deallocate(jac_element,   &
            ijac_element,  &
            r_node,        &
            rho_node,      &
            drho_node,      &
            grav_node,     &
            ellip_node,    & 
            pres_node, & 
            phi_node,      &
            sect_ind,     &
            rho20_node,     & 
            phi20_node,     & 
            pres20_node,    & 
            mu_node,       & 
            kappa_node,    &
            visco_node,    &
            ibool)
    end if
    
    
    return
  end subroutine allocate_mesh

  

  
  
  



  function find_radius_index(rr,up)
    use nrtype
    implicit none
    integer(i4b) :: find_radius_index
    real(dp), intent(in) :: rr
    integer(i4b), intent(in), optional :: up
    integer(i4b) :: iup,i
    if(present(up)) then
       iup=up
    else
       iup=1
    end if
    if(iup /= 1 .and. iup /= -1) stop 'bad input to fri'
    if(rr < 0.0_dp .or. rr > 1.0_dp) then
       find_radius_index=0
       return
    end if
    do i=1,nknot-1
       if(iup == 1) then
          if(rr == r(i) .and. rr == r(i+1)) then
             find_radius_index=i+1
             exit
          end if
          if(rr >= r(i) .and. rr < r(i+1)) then
             find_radius_index=i
             exit
          end if
       else
          if(rr == r(i) .and. rr == r(i+1)) then
             find_radius_index=i
             exit
          end if
          if(rr > r(i) .and. rr <= r(i+1)) then
             find_radius_index=i
             exit
          end if
       end if
    end do
    return
  end function find_radius_index





  subroutine integrate_between_nodes(rr,fin,fout)
    use nrtype
    implicit none
    
    real(dp), dimension(ngll), intent(in) :: rr,fin
    real(dp), dimension(ngll), intent(out) :: fout
    
    integer(i4b) :: inode,jnode,knode
    real(dp) :: rtmp,ftmp,xi,jac_tmp,int_tmp
    real(dp), dimension(ngll) :: htmp,hptmp
    
    interface lagrange_any
       
       subroutine lagrange_any(xi,NGLL,xigll,h,hprime)
         use nrtype
         implicit none         
         integer(i4b), intent(in) :: NGLL
         real(dp), intent(in) :: xi
         real(dp), dimension(NGLL), intent(in) :: xigll
         real(dp), dimension (NGLL), intent(out) :: h,hprime
       end subroutine lagrange_any
    end interface
    


    fout(1) = 0.0_dp
    do inode = 2,ngll
       int_tmp = 0.0_dp
       do jnode = 1,ngll
          rtmp = rr(inode-1)+0.5_dp*(xigll(jnode)+1.0_dp) & 
               *(rr(inode)-rr(inode-1))
          xi = 2.0_dp*(rtmp-rr(1))/(rr(ngll)-rr(1))-1.0_dp
          call lagrange_any(xi,ngll,xigll,htmp,hptmp)
          ftmp = 0.0_dp
          do knode = 1,ngll
             ftmp = ftmp + fin(knode)*htmp(knode)
          end do
          jac_tmp = 0.5_dp*(rr(inode)-rr(inode-1))
          int_tmp = int_tmp + ftmp*wgll(jnode)*jac_tmp
       end do
       fout(inode) = fout(inode-1)+int_tmp
    end do
    
    return
  end subroutine integrate_between_nodes



    subroutine mesh_finder(rtmp,ispec,xi)
      use nrtype
      real(dp), intent(in) :: rtmp
      integer(i4b), intent(out) :: ispec
      real(dp), intent(out) :: xi
      integer(i4b) :: ijspec,jspec
      real(dp) :: r1,r2


      if(rtmp < 0.0_dp .or. rtmp > r_node(ngll,nspec)) then         
         stop 'point not in the model'
      end if
      ispec = 0
      do jspec = nspec,1,-1
         r1 = r_node(1,jspec)
         r2 = r_node(ngll,jspec)
         if(rtmp >= r1 .and. rtmp <= r2) then
            ispec = jspec
            exit
         end if
      end do
      if(ispec == 0) stop 'problem in mesh_finder'
      r1 = r_node(1,ispec)
      r2 = r_node(ngll,ispec)
      xi = 2.0_dp*(rtmp-r1)/(r2-r1)-1.0_dp
      return
    end subroutine mesh_finder
    

    subroutine rad_mesh(r1,r2,dr,nr,ra)
      use nrtype
      implicit none
      real(dp), intent(in) :: r1
      real(dp), intent(in) :: r2
      real(dp), intent(in) :: dr
      integer(i4b), intent(out) :: nr
      real(dp), dimension(:), allocatable, intent(inout) :: ra



      integer(i4b) :: i,j,ir,ir1,ir2,nrl
      real(dp) :: r11,r22,drr

      real(dp), parameter :: rlil = 1.0_dp/r_norm

      
      ! locate bottom of the mesh
      ir1 = find_radius_index(r1,1)


      ! locate top of the mesh
      ir2 = find_radius_index(r2,-1)
      

      nr = 0
      do i = ir1,ir2         
         if(i == ir1) then
            if(i /= ir2) then
               r11 = r1
               r22 = r(i+1)
            else
               r11 = r1
               r22 = r2
            end if
         else if(i > ir1 .and. i < ir2) then
            r11 = r(i)
            r22 = r(i+1)
         else
            r11 = r(i)
            r22 = r2
         end if       
         if(r11 == r22) then
            nrl = 2
         else
            nrl = (r22-r11)/dr+2
         end if
         nr = nr+nrl
      end do

      allocate(ra(nr))
      

      ir = 0
      do i = ir1,ir2         
         if(i == ir1) then
            if(i /= ir2) then
               r11 = r1
               r22 = r(i+1)
            else
               r11 = r1
               r22 = r2
            end if
         else if(i > ir1 .and. i < ir2) then
            r11 = r(i)
            r22 = r(i+1)
         else
            r11 = r(i)
            r22 = r2
         end if
         if(r11 == r22) then
            nrl = 2
         else
            nrl = (r22-r11)/dr+2
         end if
         drr = (r22-r11)/(nrl-1)
         do j = 1,nrl
            ir = ir+1
            ra(ir) = r11+(j-1)*drr
         end do
      end do


      do ir = 2,nr-1
         if(ra(ir) == ra(ir+1)) then
            ra(ir) = ra(ir)-rlil
            ra(ir+1) = ra(ir+1)+rlil           
         end if
      end do



      return
    end subroutine rad_mesh




  subroutine get_depth(l,rl,ispec1)
    ! given a value of l, this routine returns the deepest
    ! radius that needs to be used in the calculations
    use nrtype
    implicit none

    integer(i4b), intent(in) :: l
    real(dp), intent(out) :: rl
    integer(i4b), intent(out) :: ispec1

    real(dp) :: rn,xi

    ! set upper limit to the sea floor
    rn = r(nknot)

    ! set the lower limit using 'rule of thumb'
    ! was 10 instead of 20
    rl = rn*(1.0_dp-20.0_dp/(l+0.5_dp))

    !rl = 0.0_dp

    if(rl < 0.0_dp) rl = 0.0_dp

    ! locate the appropriate spectral element
    call mesh_finder(rl,ispec1,xi)
    !ispec1 = 33
    
  end subroutine get_depth

  !subroutine get_depth2(l,rl)
    ! given a value of l, this routine returns the deepest
    ! radius that needs to be used in the calculations
  !  use nrtype
  !  implicit none

  !  integer(i4b), intent(in) :: l
  !  real(dp), intent(out) :: rl
    
  !  integer(i4b) :: iglob1,iglob2
  !  real(dp) :: rn,xi

    ! set upper limit to the sea floor
  !  rn = r(nsl)


    ! set the lower limit using 'rule of thumb'
  !  rl = rn*(1.0_dp-10.0_dp/(l+0.5_dp))

  !  if(rl < 0.0_dp) rl = 0.0_dp
    
    
  !end subroutine get_depth2


  function get_layer(ispec)

    ! return layer than element is in

    use nrtype

    implicit none

    integer(i4b) :: get_layer
    integer(i4b), intent(in) :: ispec

    integer(i4b) :: ilayer
    
    get_layer = 0
    do ilayer = 1,nsect
       if (ispec <= sect_ind(2,ilayer)) then
          get_layer = ilayer
          exit
       end if
    end do

  end function get_layer

  
  subroutine global_points
    
    use nrtype
    
    implicit none

    integer(i4b) :: l,ispec1
    real(dp) :: rl

    ! Calculate number of spherical and toroidal points
    nglob_ssg = gvn_ssg(ngll,nspec,1,3)
    nglob_tor = gvn_tor(ngll,nspec,1)
    nglob_ssg_0 = gvn_ssg_0(ngll,nspec,2)

    ! Calculate surface point for each l
    allocate(pa_surf(0:lmax))
    do l = 0,lmax
       call get_depth(l,rl,ispec1)
       pa_surf(l) = gvn_ssg(ngll,nspec,ispec1,1)
    end do

  end subroutine global_points


  function gvn_tor(inode,ispec,ispec1)
    use nrtype
    implicit none
    integer(i4b) :: gvn_tor
    integer(i4b), intent(in) :: inode
    integer(i4b), intent(in) :: ispec
    integer(i4b), intent(in) :: ispec1

    integer(i4b) :: spec_layer,spec1_layer,ilayer
    
    gvn_tor = 0

    if(ispec < ispec1) then
       print *, "ispec is less than ispec1."
       gvn_tor = 0
       return
    end if

    spec_layer = get_layer(ispec)

    ! gvn_tor = 0 if layer is fluid
    if (sect_ind(3,spec_layer) == 3) then
       print *, "W is undefined in fluid region."
       gvn_tor = 0
       return
    end if

    spec1_layer = get_layer(ispec1)
    
    ! if ispec1 and ispec are same layer
    if (spec_layer == spec1_layer) then
       gvn_tor = ibool(inode,ispec) - ibool(1,ispec1) + 1
       return
    end if

    ! Add points from layer ispec1 is in, 
    ! provided it's not fluid
    if (sect_ind(3,spec1_layer) /= 3) then
       gvn_tor = ibool(ngll,sect_ind(2,spec1_layer)) - ibool(1,ispec1)
    end if

    ! Add points from layers between ispec1 and ispec
    do ilayer = spec1_layer + 1, spec_layer - 1
       if (sect_ind(3,ilayer) == 3) then
          gvn_tor = gvn_tor + 1
       else
          gvn_tor = gvn_tor + ibool(ngll,sect_ind(2,ilayer)) &
                            - ibool(1,sect_ind(1,ilayer))
       end if
    end do

    ! Add points from layer ispec is in
    gvn_tor = gvn_tor + ibool(inode,ispec) - ibool(1,sect_ind(1,spec_layer)) + 1

    return
  end function gvn_tor


  function gvn_ssg(inode,ispec,ispec1,ivar)
    ! ivar == 1 ==> phi
    ! ivar == 2 ==> U
    ! ivar == 3 ==> V
    use nrtype
    implicit none
    integer(i4b) :: gvn_ssg
    integer(i4b), intent(in) :: inode
    integer(i4b), intent(in) :: ispec
    integer(i4b), intent(in) :: ispec1
    integer(i4b), intent(in) :: ivar

    integer(i4b) :: spec_layer,spec1_layer,ilayer

    if(ispec < ispec1) then
       !print *, "ispec is less than ispec1."
       gvn_ssg = 0
       return
    end if

    spec_layer = get_layer(ispec)

    ! gvn_ssg = 0 if layer is fluid and ivar = 2,3
    if ((sect_ind(3,spec_layer) == 3) .and. ((ivar == 2) .or. (ivar == 3))) then
       !print *, "U and V are undefined in fluid region."
       gvn_ssg = 0
       return
    end if

    spec1_layer = get_layer(ispec1)
    
    ! if ispec1 and ispec are same layer
    if (spec_layer == spec1_layer) then
       if (sect_ind(3,spec_layer) == 3) then
          ! layer is fluid
          gvn_ssg = ibool(inode,ispec) - ibool(1,ispec1) + 1
          return
       else
          ! layer is solid
          gvn_ssg = 3*(ibool(inode,ispec) - ibool(1,ispec1)) + ivar
          return
       end if
    end if

    ! Add points from layer ispec1 is in to point below discontinuity
    if (sect_ind(3,spec1_layer) == 3) then
       gvn_ssg = ibool(ngll,sect_ind(2,spec1_layer)) - ibool(1,ispec1)
    else
       gvn_ssg = 3*(ibool(ngll,sect_ind(2,spec1_layer)) - ibool(1,ispec1))
    end if

    ! Add points from layers between ispec1 and ispec
    do ilayer = spec1_layer + 1, spec_layer - 1
       if (sect_ind(3,ilayer) == 3) then
          gvn_ssg = gvn_ssg + ibool(ngll,sect_ind(2,ilayer)) &
                            - ibool(1,sect_ind(1,ilayer))
          if (sect_ind(3,ilayer-1) /= 3) then
             gvn_ssg = gvn_ssg + 2
          end if
       else
          gvn_ssg = gvn_ssg + 3*(ibool(ngll,sect_ind(2,ilayer)) &
                                 - ibool(1,sect_ind(1,ilayer)))
       end if
    end do

    ! Add points from layer ispec is in
    if (sect_ind(3,spec_layer) == 3) then
       gvn_ssg = gvn_ssg + ibool(inode,ispec) &
                         - ibool(1,sect_ind(1,spec_layer)) + 1
       if (sect_ind(3,spec_layer-1) /= 3) then
          gvn_ssg = gvn_ssg + 2
       end if
    else
       gvn_ssg = gvn_ssg + 3*(ibool(inode,ispec) &
                         - ibool(1,sect_ind(1,spec_layer))) + ivar
    end if
   

    return
  end function gvn_ssg



  ! degree 0 and so there is no V and ispec1 = 1
  function gvn_ssg_0(inode,ispec,ivar)
    ! ivar == 1 ==> phi
    ! ivar == 2 ==> U
    use nrtype
    implicit none
    integer(i4b) :: gvn_ssg_0
    integer(i4b), intent(in) :: inode
    integer(i4b), intent(in) :: ispec
    integer(i4b), intent(in) :: ivar

    integer(i4b) :: spec_layer,ilayer

    gvn_ssg_0 = 0

    spec_layer = get_layer(ispec)
    
    ! Add points from layers below ispec
    do ilayer = 1, spec_layer - 1
       gvn_ssg_0 = gvn_ssg_0 + 2*(ibool(ngll,sect_ind(2,ilayer)) &
                                 - ibool(1,sect_ind(1,ilayer)))
    end do

    ! Add points from layer ispec is in
    gvn_ssg_0 = gvn_ssg_0 + 2*(ibool(inode,ispec) &
                               - ibool(1,sect_ind(1,spec_layer))) + ivar
   

    return
  end function gvn_ssg_0




  subroutine set_parameters
    !===========================================================!
    ! this routine computes some of the normalization           !
    ! parameters used in the main program - it should           !
    ! always be the first routine called.                       !
    !===========================================================!
    use nrtype
    implicit none



    vel_norm=r_norm*sqrt(pi_d*bigg*rho_norm)
    acl_norm=pi_d*bigg*rho_norm*r_norm
    fre_norm=vel_norm/r_norm
    t_norm=1/fre_norm
    con_norm=vel_norm**2*rho_norm
    grav_norm=pi_d*bigg
    moment_norm=r_norm**5*rho_norm**2*grav_norm
    pibigg=pi_d*bigg/grav_norm
    rot_rate = 7.292115e-05_dp/fre_norm 
    rotfac = sqrt(fourpi_d/5.0_dp)*rot_rate*rot_rate/3.0_dp
    pot_norm = acl_norm*r_norm
    mass_norm = rho_norm*r_norm**3
    visco_norm = con_norm*t_norm
    tmin = 0.5_dp*yr2sec/t_norm
    rice = rice_ref/rho_norm
    roce = roce_ref/rho_norm
    rcru = rcru_ref/rho_norm
    return
  end subroutine set_parameters



end module module_model
