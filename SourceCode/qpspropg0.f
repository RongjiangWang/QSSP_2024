      subroutine qpspropg0(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of spheroidal response (with gravity) for degree 0
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,istp,ly,lyswap,ily,nly,ldeg,nstep,key
      real*8 f,rr1,rr2,dlnr,h
      complex*16 beta
      complex*16 y0(3),yup(3),ylw(3),coef(3,3),b(3,2)
      logical*2 sety0up,sety0lw
      external qpdifmat0
c
      if(lylwa.gt.lylw)return
c
      ldeg=0
      f=dreal(comi)/PI2
c
c===============================================================================
c
      sety0up=.false.
c
c     propagation from surface to source
c
      if(freesurf.and.lyup.eq.1)then
        yup(1)=(1.d0,0.d0)
        yup(2)=(0.d0,0.d0)
        yup(3)=(0.d0,0.d0)
      else if(lyup.lt.lyos)then
        beta=cdsqrt((2.25d0,0.d0)/crrup(lyup)**2
     &      -(comi2+(4.d0,0.d0)*cgrup(lyup)/crrup(lyup))
     &                           /cvpup(lyup)**2)
        yup(1)=(1.d0,0.d0)
        yup(2)=croup(lyup)*cvpup(lyup)**2
     &        *((1.5d0,0.d0)-beta*crrup(lyup))
        yup(3)=(0.d0,0.d0)
      else
        call qpstart0g(lyup,2,yup)
      endif
c
      do ly=lyup,lys-1
        if(lyr.eq.ly)then
          call cmemcpy(yup,y0,3)
          sety0up=.true.
        endif
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrlw(ly)/rrup(ly))/dble(nly)
        rr2=rrup(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrup(ly)*dexp(dble(ily)*dlnr)
          call ruku(yup,3,1,ly,ldeg,qpdifmat0,rr1,rr2,nruku(ldeg,ly))
        enddo
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c     yup(3)=yup(3)
c
c===============================================================================
c
      sety0lw=.false.
c
c     propagation from bottom to source
c
      if(lylw.eq.lylwb)then
        call qpstart0g(lylw,1,ylw)
      else
        call qpstart0g(lylw,0,ylw)
      endif
      if(lylw.eq.lyr.and.lylw.ge.lys)then
        call cmemcpy(ylw,y0,3)
        sety0lw=.true.
      endif
c
      do ly=lylw-1,lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call ruku(ylw,3,1,ly,ldeg,qpdifmat0,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr)then
          call cmemcpy(ylw,y0,3)
          sety0lw=.true.
        endif
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      if(sety0up.and..not.sety0lw.or..not.sety0up.and.sety0lw)then
        y0(1)=y0(1)/crrup(lyr)
        y0(2)=y0(2)/crrup(lyr)**2
c       y0(3)=y0(3)
      else
        stop 'Error #1 in qpspropg0!'
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(3,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      b(3,2)=(0.d0,0.d0)
      do i=1,3
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
        coef(i,3)=(0.d0,0.d0)
      enddo
c
c     a constant will be added to potential for region below the source
c
      coef(3,3)=-(1.d0,0.d0)
c
      key=0
      call cdsvd500(coef,b,3,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpspropg0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.lt.lys)then
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
          ypsv(3,istp)=ypsv(3,istp)+b(3,istp)
        enddo
      endif
c
      if(lylwa.le.0)then
        do istp=1,2
          ypsv(5,istp)=ypsv(3,istp)
          ypsv(6,istp)=ypsv(3,istp)/crrup(lyr)
          ypsv(3,istp)=(0.d0,0.d0)
        enddo
        return
      endif
c
c
c===============================================================================
c
      sety0lw=.false.
c
c     propagation from bottom to source
c
      lyswap=lylwa
      call qpstart0g(lyswap,1,ylw)
c
      if(lylwa.eq.lyr.and.lylwa.ge.lys)then
        call cmemcpy(ylw,y0,3)
        sety0lw=.true.
      endif
c
      do ly=lylwa-1,lys,-1
        h=rrup(ly)-rrlw(ly)
        nly=1+idint(h*(f/vpup(ly)+0.5d0/rrlw(ly)))
        dlnr=dlog(rrup(ly)/rrlw(ly))/dble(nly)
        rr2=rrlw(ly)
        do ily=1,nly
          rr1=rr2
          rr2=rrlw(ly)*dexp(dble(ily)*dlnr)
          call ruku(ylw,3,1,ly,ldeg,qpdifmat0,rr1,rr2,nruku(ldeg,ly))
        enddo
        if(ly.eq.lyr.and.ly.ge.lys)then
          call cmemcpy(ylw,y0,3)
          sety0lw=.true.
        endif
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c     ylw(3)=ylw(3)
c
      if(sety0lw)then
        y0(1)=y0(1)/crrup(lyr)
        y0(2)=y0(2)/crrup(lyr)**2
c       y0(3)=y0(3)
      else if(.not.sety0up)then
        stop 'Error #2 in qpspropg0!'
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(3,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      b(3,2)=(0.d0,0.d0)
      do i=1,3
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
        coef(i,3)=(0.d0,0.d0)
      enddo
      coef(3,3)=-(1.d0,0.d0)
      key=0
      call cdsvd500(coef,b,3,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpspropg0: anormal exit from cdsvd500!'
        return
      endif
      if(lyr.lt.lys)then
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=ypsv(i,istp)-b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,3
            ypsv(i,istp)=ypsv(i,istp)-b(2,istp)*y0(i)
          enddo
          ypsv(3,istp)=ypsv(3,istp)-b(3,istp)
        enddo
      endif
c
      do istp=1,2
        ypsv(5,istp)=ypsv(3,istp)
        ypsv(6,istp)=ypsv(3,istp)/crrup(lyr)
        ypsv(3,istp)=(0.d0,0.d0)
      enddo
c
      return
      end