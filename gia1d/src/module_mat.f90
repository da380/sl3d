module module_mat


  use nrtype
  implicit none

  complex(dpc), dimension(:,:,:), allocatable, save :: aa_ssg
  complex(dpc), dimension(:,:,:), allocatable, save :: aa_tor

  ! l = 1
  complex(dpc), dimension(:,:), allocatable, save :: aa_ssg_1
  integer(i4b), dimension(:), allocatable, save :: ipiv_ssg_1
  complex(dpc), dimension(:,:), allocatable, save :: aa_tor_1
  integer(i4b), dimension(:), allocatable, save :: ipiv_tor_1
 

contains


  subroutine construct_system_matrix
    ! The A matrices for all ls are calculated and decomposed

    use nrtype
    use module_model

    implicit none

    integer(i4b) :: l,kd,kl,ku,ldab,nrows,info,ispec1,ilayer,inode,ispec
    real(dp) :: rl

    ! l = 1
    l = 1
    call get_depth(l,rl,ispec1)

    ! Spheroidal
    kl = 3*(ngll-1) + 2
    ku = 3*(ngll-1) + 2
    ldab = 2*kl + ku + 1

    nrows = nglob_ssg

    allocate(aa_ssg_1(ldab,nrows),ipiv_ssg_1(nrows))
    aa_ssg_1 = 0.0_dp
    ipiv_ssg_1 = 0.0_dp
    call amat_ssg(ispec1,1,aa_ssg_1(1:ldab,1:nrows))

    call zgbtrf(nrows,nrows,kl,ku,aa_ssg_1(1:ldab,1:nrows),ldab,ipiv_ssg_1(1:nrows),info)

    ! Toroidal
    kl = ngll - 1
    ku = ngll - 1
    ldab = 2*kl + ku + 1

    nrows = nglob_tor

    allocate(aa_tor_1(ldab,nrows),ipiv_tor_1(nrows))
    aa_tor_1 = 0.0_dp
    ipiv_tor_1 = 0.0_dp
    call amat_tor(ispec1,1,aa_tor_1)

    call zgbtrf(nrows,nrows,kl,ku,aa_tor_1(1:ldab,1:nrows),ldab,ipiv_tor_1(1:nrows),info)

    ! Spheroidal
    kd = 3*(ngll-1) + 2
    ldab = kd + 1
    nrows = nglob_ssg
    allocate(aa_ssg(ldab,nrows,lmax))
    do l = 2,lmax

       call get_depth(l,rl,ispec1)

       nrows = gvn_ssg(ngll,nspec,ispec1,3)    

       ! Calculate spheroidal matrix

       ! Calculate amat
       call amat_ssg_sym(ispec1,l,aa_ssg(:,1:nrows,l))
    
       ! Decompose amat
       call zpbtrf('U',nrows,kd,aa_ssg(:,1:nrows,l),ldab,info)

    end do

    ! Toroidal
    nrows = nglob_tor
    kd = ngll - 1
    ldab = kd + 1
    allocate(aa_tor(ldab,nrows,lmax))
    aa_tor = 0.0_dp
    do l = 2,lmax

       call get_depth(l,rl,ispec1)   
       nrows = gvn_tor(ngll,nspec,ispec1)

       ! Calculate toroidal matrix

       ! Calculate amat
       call amat_tor_sym(ispec1,l,aa_tor(:,1:nrows,l))
    
       ! Decompose amat
       call zpbtrf('U',nrows,kd,aa_tor(:,1:nrows,l),ldab,info)

    end do
   

  end subroutine construct_system_matrix


  !==========================================================!
  !                 system matrix routines                   !
  !==========================================================!

    
  subroutine amat_tor(ispec1,l,aa)
    use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: ispec1
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(:,:), intent(out) :: aa

    integer(i4b) :: kl,ku,ldab,ia,ja,ka, & 
                    ispec,inode,jnode,knode,   & 
                    ispec11,j1,j2,ilayer
    real(dp) :: zeta2,zeta,zeta2m2,rr,rrho,rmu,tmp1,tmp2
    integer(i4b), dimension(3) :: ispec_lay
      
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zeta = sqrt(zeta2)

    kl = ngll-1
    ku = ngll-1
    ldab = 2*kl+ku+1

    ! initialize the matrix
    aa = 0.0_dp

    ia = 1
    do ilayer = 1,nsect
       ispec_lay = sect_ind(:,ilayer)
       
       if (ispec1 > ispec_lay(2)) cycle

       ! No terms in fluid regions
       if (ispec_lay(3) == 3) then
          if (ispec1 < ispec_lay(1)) then
             ia = ia + 1
          end if
          cycle
       end if

       ispec11 = max(ispec1,ispec_lay(1))

       ! begin loop over the spectral elements
       do ispec = ispec11,ispec_lay(2)

          ! begin loop over the internal nodes       
          do inode = 1,ngll

             if (inode /= 1) ia = ia + 1

             ! get row number 
             !ia2 = gvn_tor(inode,ispec,ispec1)
             ! get column number
             ja = ia

          
             ! get some parameters
             rr = r_node(inode,ispec)
             rmu = mu_node(inode,ispec)
             rrho = rho_node(inode,ispec)

             ! add the diagonal term to the matrix
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + zeta2*zeta2m2*rmu*wgll(inode) & 
                                     *jac_element(ispec)
          
             ! begin inner loop over nodes
             do jnode = 1,ngll

                ! only do the first triangle
                !if(jnode > inode) exit
                
                ! get column number 
                ja = gvn_tor(jnode,ispec,ispec1)

                ka = kl+ku+1+ia-ja
                
                ! start integration loop over nodes
                do knode = 1,ngll

                   ! get some parameters
                   rr = r_node(knode,ispec)
                   rmu = mu_node(knode,ispec)
                   
                   tmp1 = rr*hprime(knode,inode) & 
                          *ijac_element(ispec)
                   if(inode == knode) tmp1 = tmp1-1.0_dp
                
                   tmp2 = rr*hprime(knode,jnode) & 
                          *ijac_element(ispec)
                   if(jnode == knode) tmp2 = tmp2-1.0_dp
                
                   aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1*tmp2 & 
                                           *wgll(knode)*jac_element(ispec)
       
                end do
                ! end integration loop over nodes

             end do
             ! end inner loop over nodes

          end do
          ! end outer loop over nodes

       end do
       ! end loop over spectral elements

    end do
    

    ! fill in the upper triangle using symmetry
    !do ia = 1,nglob_tor-1
    !   j1 = ia+1
    !   j2 = min(nglob_tor,ia+ku)
    !   do ja = j1,j2
    !      aa(kl+ku+1+ia-ja,ja)= aa(kl+ku+1+ja-ia,ia)
    !   end do
    !end do


    return
  end subroutine amat_tor


  subroutine amat_tor_sym(ispec1,l,aa)
    use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: ispec1
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(:,:), intent(out) :: aa

    integer(i4b) :: kd,ldab,ia,ja,ka, & 
                    ispec,inode,jnode,knode,   & 
                    ispec11,j1,j2,ilayer
    real(dp) :: zeta2,zeta,zeta2m2,rr,rrho,rmu,tmp1,tmp2
    integer(i4b), dimension(3) :: ispec_lay
      
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zeta = sqrt(zeta2)

    kd = ngll-1
    ldab = kd + 1

    ! initialize the matrix
    aa = 0.0_dp

    ia = 1
    do ilayer = 1,nsect
       ispec_lay = sect_ind(:,ilayer)
       
       if (ispec1 > ispec_lay(2)) cycle

       ! No terms in fluid regions
       if (ispec_lay(3) == 3) then
          if (ispec1 < ispec_lay(1)) then
             ia = ia + 1
          end if
          cycle
       end if

       ispec11 = max(ispec1,ispec_lay(1))

       ! begin loop over the spectral elements
       do ispec = ispec11,ispec_lay(2)

          ! begin loop over the internal nodes       
          do inode = 1,ngll

             if (inode /= 1) ia = ia + 1

             ! get row number 
             !ia2 = gvn_tor(inode,ispec,ispec1)
             ! get column number
             ja = ia

          
             ! get some parameters
             rr = r_node(inode,ispec)
             rmu = mu_node(inode,ispec)
             rrho = rho_node(inode,ispec)

             ! add the diagonal term to the matrix
             ka = kd+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + zeta2*zeta2m2*rmu*wgll(inode) & 
                                     *jac_element(ispec)
          
             ! begin inner loop over nodes
             do jnode = inode,ngll
                
                ! get column number 
                ja = gvn_tor(jnode,ispec,ispec1)

                ka = kd+1+ia-ja
                
                ! start integration loop over nodes
                do knode = 1,ngll

                   ! get some parameters
                   rr = r_node(knode,ispec)
                   rmu = mu_node(knode,ispec)
                   
                   tmp1 = rr*hprime(knode,inode) & 
                          *ijac_element(ispec)
                   if(inode == knode) tmp1 = tmp1-1.0_dp
                
                   tmp2 = rr*hprime(knode,jnode) & 
                          *ijac_element(ispec)
                   if(jnode == knode) tmp2 = tmp2-1.0_dp
                
                   aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1*tmp2 & 
                                           *wgll(knode)*jac_element(ispec)
       
                end do
                ! end integration loop over nodes

             end do
             ! end inner loop over nodes

          end do
          ! end outer loop over nodes

       end do
       ! end loop over spectral elements

    end do
    

    return
  end subroutine amat_tor_sym






  subroutine amat_ssg(ispec1,l,aa)
    use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: ispec1
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(:,:), intent(out) :: aa


    integer(i4b) :: kl,ku,ldab, & 
         ia,ja,ka,ia1,ia2,ia3,ja1,ja2,ja3, & 
         ispec,inode,jnode,knode,   & 
         ispec11,ilayer
    real(dp) :: zeta2,zeta,zeta2m2,zetac,rr,rrho,rmu, & 
         rkappa,tmp1,tmp2,rgrav,ifpi,drho
    integer(i4b), dimension(3) :: ispec_lay

    
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zeta = sqrt(zeta2)
    zetac = (4.0_dp/3.0_dp)*zeta2-2.0_dp
    ifpi = 1.0_dp/(4.0_dp*pibigg)
    
    kl = 3*(ngll-1)+2
    ku = 3*(ngll-1)+2
    ldab = 2*kl+ku+1
 
    ! initialize the array
    aa = 0.0_dp

    ! if mesh doesn't start at center add inner boundary
    ! term to Poisson's equation
    if(ispec1 /= 1) then         
       ia = 1
       ja = 1
       ka = kl+ku+1+ia-ja
       rr = r_node(1,ispec1)
       aa(ka,ja) = aa(ka,ja) + l*ifpi*rr          
    end if

    do ilayer = 1,nsect
       ispec_lay = sect_ind(:,ilayer)
       if (ispec1 > ispec_lay(2)) cycle

       ispec11 = max(ispec1,ispec_lay(1))

       ! begin loop over the spectral elements
       do ispec = ispec11,ispec_lay(2)

          ! begin loop over the internal nodes       
          do inode = 1,ngll
          
             !----------------------------------!
             !  deal with non-derivative terms  !
             !----------------------------------!
             
             ! get some parameters
             rr = r_node(inode,ispec)
             rrho = rho_node(inode,ispec)
             rmu = mu_node(inode,ispec)
             rkappa = kappa_node(inode,ispec)
             rgrav = grav_node(inode,ispec)

             ia1 = gvn_ssg(inode,ispec,ispec1,1)
             ia2 = ia1 + 1
             ia3 = ia2 + 1


             !----------------------------------------!
             !    Terms common to solid and fluid     !
             !----------------------------------------!

             ! get indices for the phi-phi' term
             ia = ia1
             ja = ia1
             ka = kl+ku+1+ia-ja

             ! add contribution to the matrix from phi-phi'
             aa(ka,ja) = aa(ka,ja) & 
                         +ifpi*zeta2*wgll(inode) & 
                         *jac_element(ispec)

             !---------------------------------------!
             !           Fluid only term             !
             !---------------------------------------!
             
             if (ispec_lay(3) == 3) then
                drho = drho_node(inode,ispec)
                if (rgrav /= 0.0_dp) then
                   aa(ka,ja) = aa(ka,ja) &
                               + rr*rr*drho*wgll(inode) &
                                 *jac_element(ispec)/rgrav
                end if           

             !---------------------------------------!
             !           Solid only terms            !
             !---------------------------------------!

             else
             
                ! get indices for the phi-v'
                ia = ia1
                ja = ia3
                ka = kl+ku+1+ia-ja

                ! add contribution to the matrix from phi-v'
                aa(ka,ja) = aa(ka,ja) & 
                        +zeta2*rrho*rr*wgll(inode) &
                            *jac_element(ispec)


                ! get indices for the u-u' term
                ia = ia2
                ja = ia2
                ka = kl+ku+1+ia-ja
          
                ! add contribution to the matrix from u-u'
                aa(ka,ja) = aa(ka,ja) + & 
                            zeta2*rmu*wgll(inode)*jac_element(ispec)
                aa(ka,ja) = aa(ka,ja) + 4.0_dp*rrho*(pibigg*rrho*rr-rgrav) & 
                                      * rr*wgll(inode)*jac_element(ispec)

             

                ! get indices for the u-v' term
                ia = ia2
                ja = ia3
                ka = kl+ku+1+ia-ja

                
                ! add contribution to the matrix from u-v'
                aa(ka,ja) = aa(ka,ja) + zeta2*rrho*rgrav*rr & 
                                        * wgll(inode)*jac_element(ispec)
                         
             
                ! get indices for the v-phi'
                ia = ia3
                ja = ia1
                ka = kl+ku+1+ia-ja

                ! add contribution to the matrix from v-phi'
                aa(ka,ja) = aa(ka,ja) & 
                            +zeta2*rrho*rr*wgll(inode) &
                             *jac_element(ispec)
             
                ! get indices for the v-u' term
                ia = ia3
                ja = ia2
                ka = kl+ku+1+ia-ja
          
                ! add contribution to the matrix from v-u'
                aa(ka,ja) = aa(ka,ja) + zeta2*rrho*rgrav*rr & 
                                        *wgll(inode)*jac_element(ispec)

                         
                ! get indices for the v-v' term
                ia = ia3
                ja = ia3
                ka = kl+ku+1+ia-ja
             
                ! add contribution to the matrix from v-v'
                aa(ka,ja) = aa(ka,ja)                         & 
                            +(zeta2*zeta2*rkappa              &
                            +zeta2*zetac*rmu)                 & 
                            *wgll(inode)*jac_element(ispec)             
             end if

             !--------------------------------!
             !   deal with derivative terms   !
             !--------------------------------!
             
             do jnode = 1,ngll
              
                ja1 = gvn_ssg(jnode,ispec,ispec1,1)
                ja2 = ja1 + 1
                ja3 = ja2 + 1

                !----------------------------------------!
                !    Terms common to solid and fluid     !
                !----------------------------------------!

                ! get indices for the phi-phi' term
                ia = ia1
                ja = ja1
                ka = kl+ku+1+ia-ja

                ! add contribution to matrix from phi-phi'
                do knode = 1,ngll
                   rr = r_node(knode,ispec)
                   aa(ka,ja) = aa(ka,ja)  & 
                               +ifpi*hprime(knode,inode) & 
                               *hprime(knode,jnode)      & 
                               *rr*rr*wgll(knode)        & 
                               *ijac_element(ispec)                   

                   
                end do

                !---------------------------------------!
                !           Solid only terms            !
                !---------------------------------------!

                if (ispec_lay(3) /= 3) then
                
                   ! get indices for the phi-u' term
                   ia = ia1
                   ja = ja2
                   ka = kl+ku+1+ia-ja

                   ! add contribution to the matrix from phi-u'
                   rr = r_node(jnode,ispec)
                   rrho = rho_node(jnode,ispec)
                   aa(ka,ja) = aa(ka,ja) & 
                               +rrho*rr*rr & 
                               *hprime(jnode,inode) &
                               *wgll(jnode)
                
                   ! get indices for the u-phi' term
                   ia = ia2
                   ja = ja1
                   ka = kl+ku+1+ia-ja

                   ! add contribution to the matrix from u-phi'
                   rr = r_node(inode,ispec)
                   rrho = rho_node(inode,ispec)
                   aa(ka,ja) = aa(ka,ja) & 
                               +rrho*rr*rr & 
                               *hprime(inode,jnode) &
                               *wgll(inode)

                   ! get indices for the u-u' terms
                   ia = ia2
                   ja = ja2
                   ka = kl+ku+1+ia-ja
                
                   ! add contribution to the matrix from u-u'
                   do knode = 1,ngll
                      rr = r_node(knode,ispec)
                      rkappa = kappa_node(knode,ispec)
                      rmu = mu_node(knode,ispec)
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(inode == knode) tmp1 = tmp1+2.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(jnode == knode) tmp2 = tmp2+2.0_dp
                      aa(ka,ja) = aa(ka,ja) + rkappa*tmp1*tmp2 & 
                                  *wgll(knode)*jac_element(ispec)               
                   
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(inode == knode) tmp1 = tmp1-1.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(jnode == knode) tmp2 = tmp2-1.0_dp
                      aa(ka,ja) = aa(ka,ja) + (4.0_dp/3.0_dp) & 
                                  *rmu*tmp1*tmp2*wgll(knode)*jac_element(ispec)  
                   
                   end do
                
                
                   ! get indices for the u-v' term
                   ia = ia2
                   ja = ja3
                   ka = kl+ku+1+ia-ja
             
                   ! add contribution to the matrix from u-v'
                   rkappa = kappa_node(jnode,ispec)
                   rr = r_node(jnode,ispec)
                   rmu = mu_node(jnode,ispec)
                   tmp1 = rr*hprime(jnode,inode)*ijac_element(ispec)
                   tmp2 = tmp1
                   if(inode == jnode) tmp1 = tmp1+2.0_dp
                   aa(ka,ja) = aa(ka,ja)-zeta2*rkappa*tmp1 & 
                               *wgll(jnode)*jac_element(ispec)
                
                   if(inode == jnode) tmp2 = tmp2-1.0_dp
                   aa(ka,ja) = aa(ka,ja) + (2.0_dp/3.0_dp) & 
                               *zeta2*rmu*tmp2*wgll(jnode)*jac_element(ispec)
                
                   rmu = mu_node(inode,ispec)
                   rr = r_node(inode,ispec)
                   tmp1 = rr*hprime(inode,jnode)*ijac_element(ispec)
                   if(inode == jnode) tmp1 = tmp1-1.0_dp
                   aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1 & 
                                           *wgll(inode)*jac_element(ispec)
                
             
             
             
                   ! get indices for the v-u' term
                   ia = ia3
                   ja = ja2
                   ka = kl+ku+1+ia-ja
             
                   ! add contribution to the matrix from v-u'
                   rkappa = kappa_node(inode,ispec)
                   rr = r_node(inode,ispec) 
                   rmu = mu_node(inode,ispec)
                   tmp2 = rr*hprime(inode,jnode)*ijac_element(ispec)
                   if(inode == jnode) tmp2 = tmp2+2.0_dp
                   aa(ka,ja) = aa(ka,ja) - zeta2*rkappa*tmp2 & 
                                           *wgll(inode)*jac_element(ispec)
                
                   if(inode == jnode) tmp2 = tmp2-3.0_dp
                   aa(ka,ja) = aa(ka,ja) + (2.0_dp/3.0_dp)*zeta2*rmu & 
                                           *tmp2*wgll(inode)*jac_element(ispec)
                
                   rmu = mu_node(jnode,ispec)
                   rr = r_node(jnode,ispec)
                   tmp2 = rr*hprime(jnode,inode)*ijac_element(ispec)
                   if(inode == jnode) tmp2 = tmp2-1.0_dp
                   aa(ka,ja) = aa(ka,ja) +zeta2*rmu*tmp2 & 
                               *wgll(jnode)*jac_element(ispec)
             
             
                   ! get indices for the v-v' term
                   ia = ia3
                   ja = ja3
                   ka = kl+ku+1+ia-ja
             
               
                   ! add contribution to the matrix from v-v'
                   do knode = 1,ngll
                      rr  = r_node(knode,ispec)
                      rmu = mu_node(knode,ispec)                
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(knode == inode) tmp1 = tmp1-1.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(knode == jnode) tmp2 = tmp2-1.0_dp          
                      aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1*tmp2 &
                                              *wgll(knode)*jac_element(ispec)
                   
                   end do

                end if


             end do


          end do
          ! end loop over internal nodes
          
          
       end do
       ! end loop over the spectral elements
       
       !=========================================!
       !           Add boundary terms            !
       !=========================================!

       ! Add term from top of layer currently in
       if (ilayer /= nsect) then

          if ((ispec_lay(3) == 3) .and. (sect_ind(3,ilayer+1) /= 3)) then
             ! Fluid-solid boundary

             rr = r_node(ngll,ispec_lay(2))
             ! want density on fluid side
             rrho = rho_node(ngll,ispec_lay(2)) 
             rgrav = grav_node(ngll,ispec_lay(2))
       
             ! contribution from phi-u'
             ia = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,1)
             ja = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + rr*rr*rrho
             
       
             ! contribution from u-phi'
             ia = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ja = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,1)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + rr*rr*rrho
       
       
             ! contribution from u-u'
             ia = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ja = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + rr*rr*rgrav*rrho  

          else if ((sect_ind(3,ilayer) /= 3) .and. (sect_ind(3,ilayer+1) == 3)) then
             ! Solid-fluid boundary

             rr = r_node(ngll,ispec_lay(2))
             ! want density on fluid side
             rrho = rho_node(1,sect_ind(1,ilayer+1))  
             rgrav = grav_node(ngll,ispec_lay(2))
             
             ! contribution from phi-u'
             ia = gvn_ssg(ngll,ispec_lay(2),ispec1,1)
             ja = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) - rr*rr*rrho  
       
       
             ! contribution from u-phi'
             ia = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ja = gvn_ssg(ngll,ispec_lay(2),ispec1,1)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) - rr*rr*rrho  
       
       
             ! contribution from u-u'
             ia = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ja = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ka = kl+ku+1+ia-ja
             aa(ka,ja) = aa(ka,ja) - rr*rr*rgrav*rrho 

          end if
       end if  

    end do
    

    ! contribution to Poisson equation from exterior
    ! ispec_nsl???
    ia = gvn_ssg(ngll,nspec,ispec1,1)
    ja = gvn_ssg(ngll,nspec,ispec1,1)
    ka = kl+ku+1+ia-ja
    rr = r_node(ngll,nspec)
    aa(ka,ja) = aa(ka,ja)+(l+1)*ifpi*rr


    
    return
  end subroutine amat_ssg


 subroutine amat_ssg_sym(ispec1,l,aa)
    use nrtype
    use module_model
    implicit none
    integer(i4b), intent(in) :: ispec1
    integer(i4b), intent(in) :: l
    complex(dpc), dimension(:,:), intent(out) :: aa


    integer(i4b) :: kd,ldab, & 
         ia,ja,ka,ia1,ia2,ia3,ja1,ja2,ja3, & 
         ispec,inode,jnode,knode,   & 
         ispec11,ilayer
    real(dp) :: zeta2,zeta,zeta2m2,zetac,rr,rrho,rmu, & 
         rkappa,tmp1,tmp2,rgrav,ifpi,drho
    integer(i4b), dimension(3) :: ispec_lay

    ! Store upper diagonal components only
    ! j >= i

    
    zeta2 = l*(l+1)
    zeta2m2 = zeta2-2.0_dp
    zeta = sqrt(zeta2)
    zetac = (4.0_dp/3.0_dp)*zeta2-2.0_dp
    ifpi = 1.0_dp/(4.0_dp*pibigg)
    
    kd = 3*(ngll-1)+2
    ldab = kd+1
 
    ! initialize the array
    aa = 0.0_dp

    ! if mesh doesn't start at center add inner boundary
    ! term to Poisson's equation
    if(ispec1 /= 1) then         
       ia = 1
       ja = 1
       ka = kd+1+ia-ja
       rr = r_node(1,ispec1)
       aa(ka,ja) = aa(ka,ja) + l*ifpi*rr          
    end if

    do ilayer = 1,nsect
       ispec_lay = sect_ind(:,ilayer)
       if (ispec1 > ispec_lay(2)) cycle

       ispec11 = max(ispec1,ispec_lay(1))

       ! begin loop over the spectral elements
       do ispec = ispec11,ispec_lay(2)

          ! begin loop over the internal nodes       
          do inode = 1,ngll
          
             !----------------------------------!
             !  deal with non-derivative terms  !
             !----------------------------------!
             
             ! get some parameters
             rr = r_node(inode,ispec)
             rrho = rho_node(inode,ispec)
             rmu = mu_node(inode,ispec)
             rkappa = kappa_node(inode,ispec)
             rgrav = grav_node(inode,ispec)

             ia1 = gvn_ssg(inode,ispec,ispec1,1)
             ia2 = ia1 + 1
             ia3 = ia2 + 1


             !----------------------------------------!
             !    Terms common to solid and fluid     !
             !----------------------------------------!

             ! get indices for the phi-phi' term
             ia = ia1
             ja = ia1
             ka = kd+1+ia-ja

             ! add contribution to the matrix from phi-phi'
             aa(ka,ja) = aa(ka,ja) & 
                         +ifpi*zeta2*wgll(inode) & 
                         *jac_element(ispec)

             !---------------------------------------!
             !           Fluid only term             !
             !---------------------------------------!
             
             if (ispec_lay(3) == 3) then
                drho = drho_node(inode,ispec)
                if (rgrav /= 0.0_dp) then
                   aa(ka,ja) = aa(ka,ja) &
                               + rr*rr*drho*wgll(inode) &
                                 *jac_element(ispec)/rgrav
                end if           

             !---------------------------------------!
             !           Solid only terms            !
             !---------------------------------------!

             else
             
                ! get indices for the phi-v'
                ia = ia1
                ja = ia3
                ka = kd+1+ia-ja

                ! add contribution to the matrix from phi-v'
                aa(ka,ja) = aa(ka,ja) & 
                        +zeta2*rrho*rr*wgll(inode) &
                            *jac_element(ispec)


                ! get indices for the u-u' term
                ia = ia2
                ja = ia2
                ka = kd+1+ia-ja
          
                ! add contribution to the matrix from u-u'
                aa(ka,ja) = aa(ka,ja) + & 
                            zeta2*rmu*wgll(inode)*jac_element(ispec)
                aa(ka,ja) = aa(ka,ja) + 4.0_dp*rrho*(pibigg*rrho*rr-rgrav) & 
                                      * rr*wgll(inode)*jac_element(ispec)

             

                ! get indices for the u-v' term
                ia = ia2
                ja = ia3
                ka = kd+1+ia-ja

                
                ! add contribution to the matrix from u-v'
                aa(ka,ja) = aa(ka,ja) + zeta2*rrho*rgrav*rr & 
                                        * wgll(inode)*jac_element(ispec)
                         
                ! get indices for the v-v' term
                ia = ia3
                ja = ia3
                ka = kd+1+ia-ja
             
                ! add contribution to the matrix from v-v'
                aa(ka,ja) = aa(ka,ja)                         & 
                            +(zeta2*zeta2*rkappa              &
                            +zeta2*zetac*rmu)                 & 
                            *wgll(inode)*jac_element(ispec)             
             end if

             !--------------------------------!
             !   deal with derivative terms   !
             !--------------------------------!
             
             do jnode = inode,ngll
              
                ja1 = gvn_ssg(jnode,ispec,ispec1,1)
                ja2 = ja1 + 1
                ja3 = ja2 + 1

                !----------------------------------------!
                !    Terms common to solid and fluid     !
                !----------------------------------------!

                ! get indices for the phi-phi' term
                ia = ia1
                ja = ja1
                ka = kd+1+ia-ja

                ! add contribution to matrix from phi-phi'
                do knode = 1,ngll
                   rr = r_node(knode,ispec)
                   aa(ka,ja) = aa(ka,ja)  & 
                               +ifpi*hprime(knode,inode) & 
                               *hprime(knode,jnode)      & 
                               *rr*rr*wgll(knode)        & 
                               *ijac_element(ispec)                   

                   
                end do

                !---------------------------------------!
                !           Solid only terms            !
                !---------------------------------------!

                if (ispec_lay(3) /= 3) then
                
                   ! get indices for the phi-u' term
                   ia = ia1
                   ja = ja2
                   ka = kd+1+ia-ja

                   ! add contribution to the matrix from phi-u'
                   rr = r_node(jnode,ispec)
                   rrho = rho_node(jnode,ispec)
                   aa(ka,ja) = aa(ka,ja) & 
                               +rrho*rr*rr & 
                               *hprime(jnode,inode) &
                               *wgll(jnode)
               

                   ! get indices for the u-u' terms
                   ia = ia2
                   ja = ja2
                   ka = kd+1+ia-ja
                
                   ! add contribution to the matrix from u-u'
                   do knode = 1,ngll
                      rr = r_node(knode,ispec)
                      rkappa = kappa_node(knode,ispec)
                      rmu = mu_node(knode,ispec)
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(inode == knode) tmp1 = tmp1+2.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(jnode == knode) tmp2 = tmp2+2.0_dp
                      aa(ka,ja) = aa(ka,ja) + rkappa*tmp1*tmp2 & 
                                  *wgll(knode)*jac_element(ispec)               
                   
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(inode == knode) tmp1 = tmp1-1.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(jnode == knode) tmp2 = tmp2-1.0_dp
                      aa(ka,ja) = aa(ka,ja) + (4.0_dp/3.0_dp) & 
                                  *rmu*tmp1*tmp2*wgll(knode)*jac_element(ispec)  
                   
                   end do
                
                
                   ! get indices for the u-v' term
                   ia = ia2
                   ja = ja3
                   ka = kd+1+ia-ja
             
                   ! add contribution to the matrix from u-v'
                   rkappa = kappa_node(jnode,ispec)
                   rr = r_node(jnode,ispec)
                   rmu = mu_node(jnode,ispec)
                   tmp1 = rr*hprime(jnode,inode)*ijac_element(ispec)
                   tmp2 = tmp1
                   if(inode == jnode) tmp1 = tmp1+2.0_dp
                   aa(ka,ja) = aa(ka,ja)-zeta2*rkappa*tmp1 & 
                               *wgll(jnode)*jac_element(ispec)
                
                   if(inode == jnode) tmp2 = tmp2-1.0_dp
                   aa(ka,ja) = aa(ka,ja) + (2.0_dp/3.0_dp) & 
                               *zeta2*rmu*tmp2*wgll(jnode)*jac_element(ispec)
                
                   rmu = mu_node(inode,ispec)
                   rr = r_node(inode,ispec)
                   tmp1 = rr*hprime(inode,jnode)*ijac_element(ispec)
                   if(inode == jnode) tmp1 = tmp1-1.0_dp
                   aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1 & 
                                           *wgll(inode)*jac_element(ispec)
             
             
                   ! get indices for the v-v' term
                   ia = ia3
                   ja = ja3
                   ka = kd+1+ia-ja
             
               
                   ! add contribution to the matrix from v-v'
                   do knode = 1,ngll
                      rr  = r_node(knode,ispec)
                      rmu = mu_node(knode,ispec)                
                      tmp1 = rr*hprime(knode,inode)*ijac_element(ispec)
                      if(knode == inode) tmp1 = tmp1-1.0_dp
                      tmp2 = rr*hprime(knode,jnode)*ijac_element(ispec)
                      if(knode == jnode) tmp2 = tmp2-1.0_dp          
                      aa(ka,ja) = aa(ka,ja) + zeta2*rmu*tmp1*tmp2 &
                                              *wgll(knode)*jac_element(ispec)
                   
                   end do

                   if (jnode == inode) cycle

                   ! get indices for the u-phi' term
                   ia = ia2
                   ja = ja1
                   ka = kd+1+ia-ja

                   ! add contribution to the matrix from u-phi'
                   rr = r_node(inode,ispec)
                   rrho = rho_node(inode,ispec)
                   aa(ka,ja) = aa(ka,ja) & 
                               +rrho*rr*rr & 
                               *hprime(inode,jnode) &
                               *wgll(inode)

             
                   ! get indices for the v-u' term
                   ia = ia3
                   ja = ja2
                   ka = kd+1+ia-ja
             
                   ! add contribution to the matrix from v-u'
                   rkappa = kappa_node(inode,ispec)
                   rr = r_node(inode,ispec) 
                   rmu = mu_node(inode,ispec)
                   tmp2 = rr*hprime(inode,jnode)*ijac_element(ispec)
                   if(inode == jnode) tmp2 = tmp2+2.0_dp
                   aa(ka,ja) = aa(ka,ja) - zeta2*rkappa*tmp2 & 
                                           *wgll(inode)*jac_element(ispec)
                
                   if(inode == jnode) tmp2 = tmp2-3.0_dp
                   aa(ka,ja) = aa(ka,ja) + (2.0_dp/3.0_dp)*zeta2*rmu & 
                                           *tmp2*wgll(inode)*jac_element(ispec)
                
                   rmu = mu_node(jnode,ispec)
                   rr = r_node(jnode,ispec)
                   tmp2 = rr*hprime(jnode,inode)*ijac_element(ispec)
                   if(inode == jnode) tmp2 = tmp2-1.0_dp
                   aa(ka,ja) = aa(ka,ja) +zeta2*rmu*tmp2 & 
                               *wgll(jnode)*jac_element(ispec)


                end if


             end do


          end do
          ! end loop over internal nodes
          
          
       end do
       ! end loop over the spectral elements
       
       !=========================================!
       !           Add boundary terms            !
       !=========================================!

       ! Add term from top of layer currently in
       if (ilayer /= nsect) then

          if ((ispec_lay(3) == 3) .and. (sect_ind(3,ilayer+1) /= 3)) then
             ! Fluid-solid boundary

             rr = r_node(ngll,ispec_lay(2))
             ! want density on fluid side
             rrho = rho_node(ngll,ispec_lay(2)) 
             rgrav = grav_node(ngll,ispec_lay(2))
       
             ! contribution from phi-u'
             ia = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,1)
             ja = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ka = kd+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + rr*rr*rrho
                    
       
             ! contribution from u-u'
             ia = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ja = gvn_ssg(1,sect_ind(1,ilayer+1),ispec1,2)
             ka = kd+1+ia-ja
             aa(ka,ja) = aa(ka,ja) + rr*rr*rgrav*rrho  

          else if ((sect_ind(3,ilayer) /= 3) .and. (sect_ind(3,ilayer+1) == 3)) then
             ! Solid-fluid boundary

             rr = r_node(ngll,ispec_lay(2))
             ! want density on fluid side
             rrho = rho_node(1,sect_ind(1,ilayer+1))  
             rgrav = grav_node(ngll,ispec_lay(2))
             
             ! contribution from phi-u'
             ia = gvn_ssg(ngll,ispec_lay(2),ispec1,1)
             ja = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ka = kd+1+ia-ja
             aa(ka,ja) = aa(ka,ja) - rr*rr*rrho  
             
       
             ! contribution from u-u'
             ia = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ja = gvn_ssg(ngll,ispec_lay(2),ispec1,2)
             ka = kd+1+ia-ja
             aa(ka,ja) = aa(ka,ja) - rr*rr*rgrav*rrho 

          end if
       end if  

    end do
    

    ! contribution to Poisson equation from exterior
    ! ispec_nsl???
    ia = gvn_ssg(ngll,nspec,ispec1,1)
    ja = gvn_ssg(ngll,nspec,ispec1,1)
    ka = kd+1+ia-ja
    rr = r_node(ngll,nspec)
    aa(ka,ja) = aa(ka,ja)+(l+1)*ifpi*rr


    
    return
  end subroutine amat_ssg_sym


end module module_mat
