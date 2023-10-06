module module_fourier

  use nrtype

  integer(i4b), save :: ngl,nphi,lmax
  real(dp), save :: dphi
  real(dp), dimension(:), allocatable, save :: xgl,wgl,tgl,aphi
  
  real(dp), dimension(:,:), allocatable, save :: dnm_0
  real(dp), dimension(:,:,:), allocatable, save :: dnm_all_cf,dnm_all_fc

  contains

    subroutine fourrow_dp(data,isign)
      use nrtype; use nrutil, only : assert,swap
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: isign
      integer(i4b) :: n,i,istep,j,m,mmax,n2
      real(dp) :: theta
      complex(dpc), dimension(size(data,1)) :: temp
      complex(dpc) :: w,wp
      complex(dpc) :: ws
      n=size(data,2)
      call assert(iand(n,n-1) == 0, 'n must be a power of 2 in fourrow_dp')
      n2=n/2
      j=n2
      do i=1,n-2
         if(j > i) call swap(data(:,j+1),data(:,i+1))
         m=n2
         do
            if(m < 2 .or. j < m) exit
            j=j-m
            m=m/2
         end do
         j=j+m
      end do
      mmax=1
      do 
         if(n <= mmax) exit
         istep=2*mmax
         theta=pi_d/(isign*mmax)
         wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
         w=cmplx(1.0_dp,0.0_dp,kind=dpc)
         do m=1,mmax
            ws=w
            do i=m,n,istep
               j=i+mmax
               temp=ws*data(:,j)
               data(:,j)=data(:,i)-temp
               data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
          end do
         mmax=istep
      end do
    end subroutine fourrow_dp


    subroutine fourrow_vec_dp(data,isign)
      use nrtype; use nrutil, only : assert,swap
      implicit none
      complex(dpc), dimension(:,:,:), intent(inout) :: data
      integer(i4b), intent(in) :: isign
      integer(i4b) :: n,i,istep,j,m,mmax,n2
      real(dp) :: theta
      complex(dpc), dimension(size(data,1),size(data,2)) :: temp
      complex(dpc) :: w,wp
      complex(dpc) :: ws
      n=size(data,3)
      call assert(iand(n,n-1) == 0, 'n must be a power of 2 in fourrow_dp')
      n2=n/2
      j=n2
      do i=1,n-2
         if(j > i) call swap(data(:,:,j+1),data(:,:,i+1))
         m=n2
         do
            if(m < 2 .or. j < m) exit
            j=j-m
            m=m/2
         end do
         j=j+m
      end do
      mmax=1
      do 
         if(n <= mmax) exit
         istep=2*mmax
         theta=pi_d/(isign*mmax)
         wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
         w=cmplx(1.0_dp,0.0_dp,kind=dpc)
         do m=1,mmax
            ws=w
            do i=m,n,istep
               j=i+mmax
               temp=ws*data(:,:,j)
               data(:,:,j)=data(:,:,i)-temp
               data(:,:,i)=data(:,:,i)+temp
            end do
            w=w*wp+w
          end do
         mmax=istep
      end do
    end subroutine fourrow_vec_dp


    subroutine set_sh_grid

      use module_function

      implicit none

      integer(i4b) :: iphi
      
      !!!!!!!!!!!!!!!
      ngl = lmax + 1
      !ngl = 256
      allocate(xgl(ngl),wgl(ngl),tgl(ngl))
      call gauleg(-1.0_dp,1.0_dp,xgl,wgl)
      tgl = acos(xgl) 

      ! set number of longitudinal grid points to be lowest power of 2
      ! greater than 2*lmax
      nphi = 2.0_dp**(floor(log(2.0_dp*lmax)/log(2.0_dp))+1)
      !nphi = 512

      ! longitudinal grid points
      allocate(aphi(nphi))
      dphi = twopi_d/nphi
      do iphi = 1,nphi
         aphi(iphi) = (iphi-1)*dphi
      end do


    end subroutine set_sh_grid


    subroutine calc_sph_harm

      implicit none

      integer(i4b) :: l,igl,n,lmtot
      real(dp), dimension(-2:2,-lmax:lmax) :: dnm_tmp
      

      !allocate(dnm_all(-lmax:lmax,ngl,0:lmax,-2:2))
      !dnm_all = 0.0_dp

      allocate(dnm_0(ngl,(lmax+1)**2))
      dnm_0 = 0.0_dp

      allocate(dnm_all_cf(-2:2,ngl,(lmax+1)**2))
      dnm_all_cf = 0.0_dp

      allocate(dnm_all_fc(-2:2,(lmax+1)**2,ngl))
      dnm_all_fc = 0.0_dp

      lmtot = 1
      do l = 0,lmax
         do igl = 1,ngl
            if (l == 0) then
               call rotmx2(0,l,tgl(igl),dnm_tmp(0,-l:l),1,2*l+1)
            else if (l == 1) then
               call rotmx2(1,l,tgl(igl),dnm_tmp(-1:1,-l:l),3,2*l+1)
            else
               call rotmx2(2,l,tgl(igl),dnm_tmp(-2:2,-l:l),5,2*l+1)
            end if

            dnm_0(igl,lmtot:lmtot+2*l) = dnm_tmp(0,-l:l)

            do n = -2,2            
               dnm_all_fc(n,lmtot:lmtot+2*l,igl) = dnm_tmp(n,-l:l)
               dnm_all_cf(n,igl,lmtot:lmtot+2*l) = dnm_tmp(n,-l:l)
            end do

         end do
         lmtot = lmtot + 2*l + 1
      end do


    end subroutine calc_sph_harm



    subroutine coefs_from_fun(data,coefs)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension(ngl,nphi), intent(in) :: data
      complex(dpc), dimension((lmax+1)**2), intent(out) :: coefs

      integer(i4b) :: l,m,ilong,igl,im,lmtot,lmcur
      real(dp) :: con
      complex(dpc), dimension(ngl,nphi) :: dataloc

      dataloc = data

      ! compute phi FFT at each theta
      call fourrow_dp(dataloc,-1)
      dataloc = dataloc/nphi

      coefs = 0.0_dp

      ! compute theta integral for each l
      lmtot = 1
      do l = 0,lmax
         con = twopi_d*sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
         do m = -l,l
            !lmcur = lmtot + m + l + 1
            
            ! find position of m in data array
            if (m < 0) then
               im = nphi + m + 1
            else
               im = m + 1
            end if

            do igl = 1,ngl               
               coefs(lmtot+l+m) = coefs(lmtot+l+m) + dataloc(igl,im) &
                                             *con*dnm_0(igl,lmtot+l+m)*wgl(igl)
            end do
         end do
         lmtot = lmtot + 2*l + 1
      end do

    end subroutine coefs_from_fun


    subroutine coefs_from_fun_ten_tl(tens,coefs)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension(5,ngl,nphi), intent(in) :: tens
      complex(dpc), dimension(5,(lmax+1)**2), intent(out) :: coefs

      integer(i4b) :: l,m,ilong,igl,im,lmtot,lmcur,i
      real(dp) :: con
      complex(dpc), dimension(5,ngl,nphi) :: tensloc

      tensloc = tens

      ! compute phi FFT at each theta

      call fourrow_vec_dp(tensloc(:,:,:),-1)
      tensloc = tensloc/nphi

      coefs = 0.0_dp

      ! compute theta integral for each l
      ! !$OMP PARALLEL DO
      do l = 0,lmax
         lmtot = l**2 + 1
         con = twopi_d*sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
         do m = -l,l               
            ! find position of m in data array
            if (m < 0) then
               im = nphi + m + 1
            else
               im = m + 1
            end if
            do igl = 1,ngl    
               coefs(1,lmtot+l+m) = coefs(1,lmtot+l+m) + tensloc(1,igl,im)*wgl(igl) &
                                                         *con*dnm_all_cf(-2,igl,lmtot+l+m)
               coefs(2,lmtot+l+m) = coefs(2,lmtot+l+m) + tensloc(2,igl,im)*wgl(igl) &
                                                         *con*dnm_all_cf(-1,igl,lmtot+l+m)
               coefs(3,lmtot+l+m) = coefs(3,lmtot+l+m) + tensloc(3,igl,im)*wgl(igl) &
                                                         *con*dnm_all_cf(0,igl,lmtot+l+m)
               coefs(4,lmtot+l+m) = coefs(4,lmtot+l+m) + tensloc(4,igl,im)*wgl(igl) &
                                                         *con*dnm_all_cf(1,igl,lmtot+l+m)
               coefs(5,lmtot+l+m) = coefs(5,lmtot+l+m) + tensloc(5,igl,im)*wgl(igl) &
                                                         *con*dnm_all_cf(2,igl,lmtot+l+m)
            end do
         end do
      end do
      ! !$OMP END PARALLEL DO
      
    end subroutine coefs_from_fun_ten_tl



    subroutine coefs_from_fun_parr(data,coefs)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension(ngl,nphi), intent(in) :: data
      complex(dpc), dimension((lmax+1)**2), intent(out) :: coefs

      integer(i4b) :: l,m,ilong,igl,im,lmtot,lmcur
      real(dp) :: con
      complex(dpc), dimension(ngl,nphi) :: dataloc

      dataloc = data

      ! compute phi FFT at each theta
      call fourrow_dp(dataloc,-1)
      dataloc = dataloc/nphi

      coefs = 0.0_dp

      ! compute theta integral for each l
      ! !$OMP PARALLEL DO
      do l = 0,lmax
         lmtot = l**2 + 1
         con = twopi_d*sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
            do m = -l,l
               ! find position of m in data array
               if (m < 0) then
                  im = nphi + m + 1
               else
                  im = m + 1
               end if
               do igl = 1,ngl               
               coefs(lmtot+l+m) = coefs(lmtot+l+m) + dataloc(igl,im) &
                                                     *con*dnm_0(igl,lmtot+l+m)*wgl(igl)
            end do
         end do
      end do
      ! !$OMP END PARALLEL DO

    end subroutine coefs_from_fun_parr


    subroutine zero_coef_from_fun(data,coef)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension(ngl,nphi), intent(in) :: data
      complex(dpc), intent(out) :: coef

      integer(i4b) :: l,m,ilong,igl,im,lmtot,lmcur
      real(dp) :: con
      complex(dpc), dimension(ngl,nphi) :: dataloc

      dataloc = data

      ! compute phi FFT at each theta
      call fourrow_dp(dataloc,-1)
      dataloc = dataloc/nphi

      coef = 0.0_dp

      ! compute theta integral
      con = twopi_d/sqrt(fourpi_d)
      do igl = 1,ngl
         coef = coef + dataloc(igl,1)*con*dnm_0(igl,1)*wgl(igl)
      end do


    end subroutine zero_coef_from_fun


    subroutine fun_from_coefs(coefs,displ)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension((lmax+1)**2), intent(in) :: coefs
      complex(dpc), dimension(ngl,nphi), intent(out) :: displ

      integer(i4b) :: l,m,igl,im,lmtot,lmcur
      real(dp) :: con

      displ = 0.0_dp

      lmtot = 1
      do l = 0,lmax
         con = sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
         do m = -l,l
            !lmcur = lmtot + m + l + 1
               
            ! find position of m in data array
            if (m < 0) then
               im = nphi + m + 1
            else
               im = m + 1
            end if
            do igl = 1,ngl
               displ(igl,im) = displ(igl,im) + coefs(lmtot+l+m)*dnm_0(igl,lmtot+l+m)*con
            end do
         end do
         lmtot = lmtot + 2*l + 1
      end do

      call fourrow_dp(displ,1)      

    end subroutine fun_from_coefs


    subroutine fun_from_coefs_parr(coefs,displ)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension((lmax+1)**2), intent(in) :: coefs
      complex(dpc), dimension(ngl,nphi), intent(out) :: displ

      integer(i4b) :: l,m,igl,im,lmtot,lmcur
      real(dp) :: con

      displ = 0.0_dp

      lmtot = 1
      ! !$OMP PARALLEL DO
      do igl = 1,ngl
         do l = 0,lmax
            lmtot = l**2 + 1
            con = sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
            do m = -l,l
               ! find position of m in data array
               if (m < 0) then
                  im = nphi + m + 1
               else
                  im = m + 1
               end if

               displ(igl,im) = displ(igl,im) + coefs(lmtot+l+m)*dnm_0(igl,lmtot+l+m)*con
            end do

         end do
      end do
      ! !$OMP END PARALLEL DO

      call fourrow_dp(displ,1)      

    end subroutine fun_from_coefs_parr


    subroutine fun_from_coefs_ten(coefs,tens)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension((lmax+1)**2,6), intent(in) :: coefs
      complex(dpc), dimension(ngl,nphi,6), intent(out) :: tens

      integer(i4b) :: l,m,igl,im,lmtot,lmcur,i,n
      real(dp) :: con

      ! t(1) = t(-,-)
      ! t(2) = t(-,0)
      ! t(3) = t(-,+)
      ! t(4) = t(0,0)
      ! t(5) = t(0,+)
      ! t(6) = t(+,+)

      tens = 0.0_dp

      n = -2
      do i = 1,6
         lmtot = 1
         do l = 0,lmax
            con = sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
            do igl = 1,ngl
               do m = -l,l
                  !lmcur = lmtot + m + l + 1
               
                  ! find position of m in data array
                  if (m < 0) then
                     im = nphi + m + 1
                  else
                     im = m + 1
                  end if


                  tens(igl,im,i) = tens(igl,im,i) + coefs(lmtot+l+m,i) &
                                                    *dnm_all_cf(igl,lmtot+l+m,n)*con

               end do
            end do
            lmtot = lmtot + 2*l + 1
         end do
         call fourrow_dp(tens(:,:,i),1)      

         if (i /= 3) n = n+1

      end do

    end subroutine fun_from_coefs_ten


    subroutine fun_from_coefs_ten_tl(coefs,tens)

      use nrtype
      use module_util

      implicit none

      complex(dpc), dimension(5,(lmax+1)**2), intent(in) :: coefs
      complex(dpc), dimension(5,ngl,nphi), intent(out) :: tens

      integer(i4b) :: l,m,igl,im,lmtot,lmcur,i,n
      real(dp) :: con
      complex(dpc), dimension(5,nphi,ngl) :: tens_tmp

      ! t(1) = t(-,-)
      ! t(2) = t(-,0)
      ! t(3) = t(-,+) = 0.5*t(0,0)
      ! t(4) = t(0,+)
      ! t(5) = t(+,+)

      tens = 0.0_dp
      tens_tmp = 0.0_dp
      
      ! !$OMP PARALLEL DO
      do igl = 1,ngl
         do l = 0,lmax
            lmtot = l**2 + 1
            con = sqrt((2.0_dp*l+1.0_dp)/fourpi_d)
            do m = -l,l
            
               ! find position of m in data array
               if (m < 0) then
                  im = nphi + m + 1
               else
                  im = m + 1
               end if
               tens_tmp(1,im,igl) = tens_tmp(1,im,igl) + coefs(1,lmtot+l+m) &
                                                         *dnm_all_fc(-2,lmtot+l+m,igl)*con
               tens_tmp(2,im,igl) = tens_tmp(2,im,igl) + coefs(2,lmtot+l+m) &
                                                         *dnm_all_fc(-1,lmtot+l+m,igl)*con
               tens_tmp(3,im,igl) = tens_tmp(3,im,igl) + coefs(3,lmtot+l+m) &
                                                         *dnm_all_fc(0,lmtot+l+m,igl)*con
               tens_tmp(4,im,igl) = tens_tmp(4,im,igl) + coefs(4,lmtot+l+m) &
                                                         *dnm_all_fc(1,lmtot+l+m,igl)*con
               tens_tmp(5,im,igl) = tens_tmp(5,im,igl) + coefs(5,lmtot+l+m) &
                                                         *dnm_all_fc(2,lmtot+l+m,igl)*con

            end do
         end do
      end do
      ! !$OMP END PARALLEL DO      

      ! !$OMP PARALLEL DO
      do i = 1,5
         tens(i,:,:) = transpose(tens_tmp(i,:,:))
      end do
      ! !$OMP END PARALLEL DO

      call fourrow_vec_dp(tens(:,:,:),1)      

      
    end subroutine fun_from_coefs_ten_tl


    subroutine high_pass(data,i1,i2)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: i1,i2
      integer(i4b) :: i
      real(dp) :: r1

      if(i1 >= i2) stop 'bad input to low_pass'

      do i=1,min(i2,size(data,2)/2+1)
         if(i < i1) then
            data(:,i)=0.0_dp
         else
            r1=pi_d*real(i-i1)/real(i2-i1)
            r1=0.5_dp*(1.0_dp-cos(r1))
            data(:,i)=r1*data(:,i)
         end if
      end do

      return
    end subroutine high_pass



    subroutine low_pass(data,i1,i2)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      integer(i4b), intent(in) :: i1,i2
      integer(i4b) :: i
      real(dp) :: r1

      if(i1 > i2) stop 'bad input to low_pass'

      do i=i1,size(data,2)/2+1
         if(i > i2) then
            data(:,i)= 0.0_dp
         else
            r1=pi_d*real(i2-i)/real(i2-i1)
            r1=0.5_dp*(1.0_dp-cos(r1))
            data(:,i)=r1*data(:,i)
         end if
      end do


         return
    end subroutine low_pass

    subroutine time_chop(data,dt,ts,t1,t2,t3,t4)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt,ts,t1,t2,t3,t4
      integer(i4b) :: i,n,i1,i2
      real(dp) :: arg,t

      n=size(data,2)

      do i=1,n
         t=ts+(i-1)*dt
         if(t < t1) then
            data(:,i)=0.0_dp
         else if(t >= t1 .and. t < t2) then
            arg=pi_d*(t-t1)/(t2-t1)
            arg=0.5_dp*(1.0_dp-cos(arg))
            data(:,i)=arg*data(:,i)
         else if(t > t3 .and. t <= t4) then
            arg=pi_d*(t4-t)/(t4-t3)
            arg=0.5_dp*(1.0_dp-cos(arg))
            data(:,i)=arg*data(:,i)
         else if(t > t4) then
            data(:,i)=0.0_dp
         end if
         
      end do
      
      return
    end subroutine time_chop


    subroutine convolve(data,dt,tau)
      use nrtype; use module_util
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), intent(in) :: tau
      real(dp) :: w,dw,sinc

      n=size(data,2)
      dw=1.0_dp/(dt*n)

      do i=2,n/2+1
         w=twopi_d*(i-1)*dw
         sinc=sin(w*tau)/(w*tau)
         data(:,i)=data(:,i)*sinc
      end do
         
      return
    end subroutine convolve


  subroutine rick_convolve(data,dt,tau)
      use nrtype; use module_util
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), intent(in) :: tau
      real(dp) :: w,dw,sinc,t
      complex(dpc), dimension(:,:), allocatable :: ra

      n=size(data,2)
      dw=1.0_dp/(dt*n)
      allocate(ra(1,n))
      do i=1,n
         t=(i-1)*dt
         ra(1,i)=rick(t,tau)
      end do
      call fourrow_dp(ra,1)

      do i=2,n/2+1
         data(:,i)=data(:,i)*ra(1,i)
      end do
         
      return
    end subroutine rick_convolve


    subroutine accel(data,dt)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt
      integer(i4b) :: i,n
      real(dp), parameter :: tau=10.0_dp
      real(dp) :: w,dw,sinc

      n=size(data,2)
      dw=1.0_dp/(dt*n)

      do i=2,n/2+1
         w=twopi_d*(i-1)*dw
         data(:,i)=-w**2*data(:,i)
      end do
         
      return
    end subroutine accel


    subroutine filter(data,dt,ts,w1,w2,w3,w4,acc,box,ric)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(inout) :: data
      real(dp), intent(in) :: dt,ts,w1,w2,w3,w4
      logical(lgt), intent(in), optional :: acc
      real(dp), intent(in), optional :: box,ric
      integer(i4b) :: n,i1,i2,i3,i4,i
      real(dp) :: dw,w,tau,t1,t2,t3,t4
      real(dp), parameter :: frac=0.8
      complex(dpc) :: filt


      n=size(data,2)
      dw=1.0_dp/(n*dt)
      i1=floor(w1/dw)+1
      i2=floor(w2/dw)+1
      i3=floor(w3/dw)+1
      i4=floor(w4/dw)+1
      t4=ts+dt*(size(data,2)-1)
      t3=frac*t4
      t1=ts
      t2=(1.0_dp-frac)*t4
      call time_chop(data,dt,ts,t1,t2,t3,t4)
      call fourrow_dp(data,1)
      if(i1 /= i2) then
         call high_pass(data,i1,i2)
      end if
      if(i3 /= i4) then
         call low_pass(data,i3,i4)
      end if
      if(acc) then
         call accel(data,dt)
      end if
      if(present(box)) then
         if(box > 0.0_dp) then
            tau=box
            call convolve(data,dt,tau)
         end if
      end if
      if(present(ric)) then
         if(ric > 0.0_dp) then
            tau=ric
            call rick_convolve(data,dt,tau)
         end if
      end if
      data(:,n/2+2:n)=conjg(data(:,n/2:2:-1))
      call fourrow_dp(data,-1)
      data=data/n

      return
    end subroutine filter



    subroutine get_power_spectra(data,dt,power,freq)
      use nrtype
      implicit none
      complex(dpc), dimension(:,:), intent(in) :: data
      real(dp), intent(in) :: dt
      real(dp), dimension(:), allocatable, intent(out) :: freq
      real(dp), dimension(:,:), allocatable, intent(out) :: power
      integer(i4b) :: i,n
      real(dp) :: nfreq
      complex(dpc), dimension(:,:), allocatable :: data_tmp
  
      n=size(data,2)
      allocate(data_tmp(size(data,1),n))
      data_tmp=data
      call fourrow_dp(data_tmp,1)
      allocate(power(size(data,1),n/2+1),freq(n/2+1))
      power(:,1:n/2+1)=data_tmp(:,1:n/2+1)*conjg(data_tmp(:,1:n/2+1))
      do i=1,n/2+1
         freq(i)=(i-1)/(n*dt)
      end do
      deallocate(data_tmp)
      return
    end subroutine get_power_spectra

    
    function  rick(t,tau)
      use nrtype
      implicit none
      real(dp) :: rick
      real(dp), intent(in)  :: t,tau
      real(dp) :: alpha = 0.75_dp
      rick=-(t/tau-alpha)*exp(-2.0*pi_d**2*(t/tau-alpha)**2)
      rick = twopi_d*rick*exp(0.5_dp)
      return
    end function rick
    
    
    function hann(t,t11,t12,t21,t22)
      use nrtype
      implicit none
      logical(lgt) :: ltmp
      real(dp) :: hann
      real(dp), intent(in) :: t,t11,t12,t21,t22
      
      ltmp = (t11 == 0.0_dp .and. t12 == 0.0_dp & 
           .and. t21 == 0.0_dp .and. t22 == 0.0_dp)
      if(.not.ltmp) then
         if(t < t11) then
            hann = 0.0_dp
         else if(t >= t11 .and. t < t12) then
            hann = pi_d*(t-t11)/(t12-t11)
            hann = 0.5_dp*(1.0_dp-cos(hann))
         else if(t >= t12 .and. t < t21) then
            hann = 1.0_dp
         else if(t >= t21 .and. t < t22) then
            hann = pi_d*(t22-t)/(t22-t21)
            hann = 0.5_dp*(1.0_dp-cos(hann))
         else if(t >= t22) then
            hann = 0.0_dp
         end if
      end if
      return
    end function hann
    


end module module_fourier
