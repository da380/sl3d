c--------------------------------------------------------------
c
c     subroutine prot(l,theta,nmax,mmax,dd,w1,w2)
c
c     Uses an l-stepping recursion to calculate elements of
c     the rotation matrices d^{(l)}_{Nm}(\theta) for N between
c     -nmax and nmax and m between -mmax and mmax
c     the result is in dd. w1 and w2 are work arrays that are
c     used to store the results for the value of l and l-1 just
c     calculated. When prot is called it checks to see if l is
c     greater than the previous l, and if so it does only the
c     necessary number of iterations to increment l from the
c     previous value to the new value. If l is not greater than
c     the previous value it reinitializes and iterates from l=0.
c     Therefore, if successive calls to prot are for successive
c     values of l, only one l-step is performed. It is
c     assumed that the contents of dd,w1,w2 and the values
c     nmax,mmax remain unchanged between calls. If theta changes
c     it reinitializes and iterates from l=0. The routine can
c     be forced to reinitialize by calling it with l=0.
c
      subroutine prot(l,theta,nmax,mmax,dd,w1,w2)
      implicit double precision (a-h,o-z)
      save llast,thlast
      data llast/-1/
      dimension dd(-nmax:nmax,-mmax:mmax),
     1          w1(-nmax:nmax,-mmax:mmax),
     1          w2(-nmax:nmax,-mmax:mmax)

      cost=dcos(theta)
      cht=dcos(theta/2)
      sht=dsin(theta/2)
      if(l.le.llast.or.theta.ne.thlast) llast=-1 ! reinitialize

      do ll=llast+1,l
        xll=ll
        if(ll.eq.0) then ! initialize
          do i=-nmax,nmax
          do j=-mmax,mmax
            dd(i,j)=0.d0
            w1(i,j)=0.d0
            w2(i,j)=0.d0
          enddo
          enddo
          dd(0,0)=1.d0
        else
          do j=0,mmax
          do i=-j,j
            xi=i
            xj=j
            if(j.lt.ll) then
              if(ll.eq.1) then ! l-stepping formula doesnt work
                dd(i,j)=cost
              else
                dd(i,j)=
     1          (-(xll*dsqrt((-1- xi+xll)*(-1+xi+xll)*(-1-xj 
     1                             + xll)*(-1+xj+xll))*w1(i,j))
     1          + (-1+2*xll)*(-(xi*xj) + (-1+xll)*xll*cost)*w2(i,j))/
     1          ((-1+xll)*dsqrt((-xi+xll)*(xi+xll)*(-xj+xll)*(xj+xll)))

              endif
            else if(j.eq.ll.and.(i.ge.-j.and.i.le.j)) then
              dd(i,j)=(-1)**(ll-i)
     1            *dsqrt(fac(2*ll)/(fac(ll-i)*fac(ll+i)))
     1            *cht**(ll+i)*sht**(ll-i) 
            endif
          enddo
          enddo
        endif
        do j=0,mmax
        do i=-j,j
          w1(i,j)=w2(i,j)
          w2(i,j)=dd(i,j)
        enddo
        enddo
      enddo
c     fill the rest of the array
      do j=1,mmax
        isg=1
        do i=-j,j
          dd(j,i)=isg*dd(i,j)
          dd(-i,-j)=isg*dd(i,j)
          dd(-j,-i)=dd(i,j)
          isg=-isg
        enddo
      enddo


      llast=l
      thlast=theta
      return
      end

      double precision function fac(n)
c
c     factorial function
c
      if(n.lt.0) then
        fac=0
      else
        fac=1
        do i=1,n
          fac=fac*i
        enddo
      endif

      return
      end


    
