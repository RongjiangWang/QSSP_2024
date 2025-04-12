      subroutine qpsprop(ypsv,ldeg,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of spheroidal response
c     ypsv(6,4): solution vector (complex)
c
      integer*4 ldeg,lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,j,j0,istp,ly,lyswap,key
      real*8 y4max
      complex*16 cdet,cyswap,uc,cros
      complex*16 c0(6,3),c1(6,3),cc0(4,2),cc1(4,2)
      complex*16 y0(6,3),y1(6,3)
      complex*16 yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      complex*16 wave(6),orth(3,3),orthc(2,2)
      complex*16 coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
      logical*2 sety0up,sety0lw
c
      complex*16 c2,c3,c6
      data c2,c3,c6/(2.d0,0.d0),(3.d0,0.d0),(6.d0,0.d0)/
c
c===============================================================================
c
      sety0up=.false.
c
c     propagation from surface to atmosphere/ocean bottom
c
      if(lylwa.gt.lylw)return
c
      if(lyob.gt.lyup)then
        do j=1,2
          do i=1,4
            yupc(i,j)=(0.d0,0.d0)
          enddo
        enddo
c
        if(freesurf.and.lyup.eq.1)then
          yupc(1,1)=(1.d0,0.d0)
          yupc(3,2)=(1.d0,0.d0)
        else
          call qpstart4a(ldeg,lyup,2,yupc)
        endif
      endif
c
      do ly=lyup,min0(lyob,lys)-1
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=yupc(1,j)
            y0(2,j)=yupc(2,j)
            y0(3,j)=-yupc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=yupc(3,j)
            y0(6,j)=yupc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0up=.true.
        endif
c
        do j=1,3,2
          wave(j)=cdexp(cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),yupc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(2,1)*cc0(4,2)-cc0(4,1)*cc0(2,2)
        orthc(1,1)=cc0(4,2)/cdet
        orthc(1,2)=-cc0(2,2)/cdet
        orthc(2,1)=-cc0(4,1)/cdet
        orthc(2,2)=cc0(2,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0up)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
        cc1(1,1)=cc1(1,1)*wave(1)*wave(2)
        cc1(2,1)=(1.d0,0.d0)
        cc1(3,1)=cc1(3,1)*wave(3)*wave(2)
        cc1(4,1)=(0.d0,0.d0)
c
        cc1(1,2)=cc1(1,2)*wave(1)*wave(4)
        cc1(2,2)=(0.d0,0.d0)
        cc1(3,2)=cc1(3,2)*wave(3)*wave(4)
        cc1(4,2)=(1.d0,0.d0)
c
        call caxcb(mas4x4lw(1,1,ly),cc1,4,4,2,yupc)
      enddo
c
      if(lys.ge.lyob)then
        if(lyup.lt.lyob)then
          do j=1,2
            yup(1,j)=yupc(1,j)
            yup(2,j)=yupc(2,j)
            yup(3,j)=(0.d0,0.d0)
            yup(4,j)=(0.d0,0.d0)
            yup(5,j)=yupc(3,j)
            yup(6,j)=yupc(4,j)
          enddo
          yup(1,3)=(0.d0,0.d0)
          yup(2,3)=(0.d0,0.d0)
          yup(3,3)=(1.d0,0.d0)
          yup(4,3)=(0.d0,0.d0)
          yup(5,3)=(0.d0,0.d0)
          yup(6,3)=(0.d0,0.d0)
        else if(lyup.lt.lycm)then
          if(freesurf.and.lyup.eq.1)then
            do j=1,3
              do i=1,6
                yup(i,j)=(0.d0,0.d0)
              enddo
            enddo
            yup(1,1)=(1.d0,0.d0)
            yup(3,2)=(1.d0,0.d0)
            yup(5,3)=(1.d0,0.d0)
          else
            call qpstart6a(ldeg,lyup,2,yup)
          endif
        endif
      endif
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source or core-mantle boundary
c
      do ly=max0(lyup,lyob),min0(lycm,lys)-1
        if(ly.eq.lyr)then
          call cmemcpy(yup,y0,18)
          sety0up=.true.
        endif
c
        do j=1,5,2
          wave(j)=cdexp(cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),yup,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0up)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
        c1(1,1)=c1(1,1)*wave(1)*wave(2)
        c1(2,1)=(1.d0,0.d0)
        c1(3,1)=c1(3,1)*wave(3)*wave(2)
        c1(4,1)=(0.d0,0.d0)
        c1(5,1)=c1(5,1)*wave(5)*wave(2)
        c1(6,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*wave(1)*wave(4)
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*wave(3)*wave(4)
        c1(4,2)=(1.d0,0.d0)
        c1(5,2)=c1(5,2)*wave(5)*wave(4)
        c1(6,2)=(0.d0,0.d0)
c
        c1(1,3)=c1(1,3)*wave(1)*wave(6)
        c1(2,3)=(0.d0,0.d0)
        c1(3,3)=c1(3,3)*wave(3)*wave(6)
        c1(4,3)=(0.d0,0.d0)
        c1(5,3)=c1(5,3)*wave(5)*wave(6)
        c1(6,3)=(1.d0,0.d0)
c
        call caxcb(mas6x6lw(1,1,ly),c1,6,6,3,yup)
      enddo
c===============================================================================
c
c     propagation from core-mantle boundary to source or solid core surface
c
      if(lys.ge.lycm)then
        if(lyup.lt.lycm)then
          y4max=cdabs(yup(4,3))
          j0=3
          do j=1,2
            if(y4max.lt.cdabs(yup(4,j)))then
              y4max=cdabs(yup(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=yup(i,j0)
            yup(i,j0)=yup(i,3)
            yup(i,3)=cyswap
          enddo
          do j=1,2
            yupc(1,j)=yup(1,j)-yup(4,j)*yup(1,3)/yup(4,3)
            yupc(2,j)=yup(2,j)-yup(4,j)*yup(2,3)/yup(4,3)
            yupc(3,j)=yup(5,j)-yup(4,j)*yup(5,3)/yup(4,3)
            yupc(4,j)=yup(6,j)-yup(4,j)*yup(6,3)/yup(4,3)
          enddo
          if(sety0up)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)-yup(4,j)*y0(i,3)/yup(4,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        else if(lyup.lt.lycc)then
          call qpstart4a(ldeg,lyup,2,yupc)
        endif
      endif
c
      do ly=lycm,min0(lycc,lys)-1
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=yupc(1,j)
            y0(2,j)=yupc(2,j)
            y0(3,j)=-yupc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=yupc(3,j)
            y0(6,j)=yupc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0up=.true.
        endif
c
        do j=1,3,2
          wave(j)=cdexp(cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),yupc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(2,1)*cc0(4,2)-cc0(4,1)*cc0(2,2)
        orthc(1,1)=cc0(4,2)/cdet
        orthc(1,2)=-cc0(2,2)/cdet
        orthc(2,1)=-cc0(4,1)/cdet
        orthc(2,2)=cc0(2,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0up)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
        cc1(1,1)=cc1(1,1)*wave(1)*wave(2)
        cc1(2,1)=(1.d0,0.d0)
        cc1(3,1)=cc1(3,1)*wave(3)*wave(2)
        cc1(4,1)=(0.d0,0.d0)
c
        cc1(1,2)=cc1(1,2)*wave(1)*wave(4)
        cc1(2,2)=(0.d0,0.d0)
        cc1(3,2)=cc1(3,2)*wave(3)*wave(4)
        cc1(4,2)=(1.d0,0.d0)
c
        call caxcb(mas4x4lw(1,1,ly),cc1,4,4,2,yupc)
      enddo
c===============================================================================
c
c     propagation from solid core surface to source
c
      if(lys.ge.lycc)then
        if(lyup.lt.lycc)then
          do j=1,2
            yup(1,j)=yupc(1,j)
            yup(2,j)=yupc(2,j)
            yup(3,j)=(0.d0,0.d0)
            yup(4,j)=(0.d0,0.d0)
            yup(5,j)=yupc(3,j)
            yup(6,j)=yupc(4,j)
          enddo
          if(lys.ge.lycc)then
            yup(1,3)=(0.d0,0.d0)
            yup(2,3)=(0.d0,0.d0)
            yup(3,3)=(1.d0,0.d0)
            yup(4,3)=(0.d0,0.d0)
            yup(5,3)=(0.d0,0.d0)
            yup(6,3)=(0.d0,0.d0)
          endif
        else if(lyup.lt.ly0)then
          call qpstart6a(ldeg,lyup,2,yup)
        endif
      endif
c
      do ly=lycc,lys-1
        if(ly.eq.lyr)then
          call cmemcpy(yup,y0,18)
          sety0up=.true.
        endif
c
        do j=1,5,2
          wave(j)=cdexp(cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),yup,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0up)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
        c1(1,1)=c1(1,1)*wave(1)*wave(2)
        c1(2,1)=(1.d0,0.d0)
        c1(3,1)=c1(3,1)*wave(3)*wave(2)
        c1(4,1)=(0.d0,0.d0)
        c1(5,1)=c1(5,1)*wave(5)*wave(2)
        c1(6,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*wave(1)*wave(4)
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*wave(3)*wave(4)
        c1(4,2)=(1.d0,0.d0)
        c1(5,2)=c1(5,2)*wave(5)*wave(4)
        c1(6,2)=(0.d0,0.d0)
c
        c1(1,3)=c1(1,3)*wave(1)*wave(6)
        c1(2,3)=(0.d0,0.d0)
        c1(3,3)=c1(3,3)*wave(3)*wave(6)
        c1(4,3)=(0.d0,0.d0)
        c1(5,3)=c1(5,3)*wave(5)*wave(6)
        c1(6,3)=(1.d0,0.d0)
c
        call caxcb(mas6x6lw(1,1,ly),c1,6,6,3,yup)
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        if(lys.gt.lyup.and.lys.ne.lycm)then
          cros=cro(lys-1)
        else
          cros=cro(lys)
        endif
        do j=1,2
          yup(1,j)=yupc(1,j)
          yup(2,j)=yupc(2,j)
          yup(3,j)=-yupc(2,j)/(cros*comi2*crrup(lys)**2)
          yup(4,j)=(0.d0,0.d0)
          yup(5,j)=yupc(3,j)
          yup(6,j)=yupc(4,j)
        enddo
        do i=1,6
          yup(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        yup(1,j)=yup(1,j)/crrup(lys)
        yup(2,j)=yup(2,j)/crrup(lys)**2
        yup(3,j)=yup(3,j)/crrup(lys)
        yup(4,j)=yup(4,j)/crrup(lys)**2
c       yup(5,j)=yup(5,j)
        yup(6,j)=yup(6,j)/crrup(lys)
      enddo
c
c===============================================================================
c
      sety0lw=.false.
c
c     propagation within inner core
c
      if(lylw.ge.lycc)then
c
c       lowest layer is within inner core
c
        if(lylw.eq.lylwb)then
          call qpstart6a(ldeg,lylw,1,ylw)
        else
          do j=1,3
            do i=1,6
              ylw(i,j)=mas6x6up(i,2*j-1,lylw)
            enddo
          enddo
        endif
        if(lyr.eq.lylw.and.lylw.ge.lys)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      endif
c
      do ly=lylw-1,max0(lycc,lys),-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(ly.eq.lyr)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lys.lt.lycc)then
        if(lylw.ge.lycc)then
c
c         interface conditions: solid to liquid
c
          y4max=cdabs(ylw(4,3))
          j0=3
          do j=1,2
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
            ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
            ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
            ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
          enddo
          if(sety0lw)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        else if(lylw.ge.lycm)then
c
c         lowest layer is within the liquid core
c
          if(lylw.eq.lylwb)then
            call qpstart4a(ldeg,lylw,1,ylwc)
          else
            do j=1,2
              do i=1,4
                ylwc(i,j)=mas4x4up(i,2*j-1,lylw)
              enddo
            enddo
          endif
c
          if(lylw.eq.lyr.and.lylw.ge.lys)then
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=ylwc(2,j)
              y0(3,j)=-ylwc(2,j)/(cro(lylw)*comi2*crrup(lylw)**2)
              y0(4,j)=(0.d0,0.d0)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylw,lycc)-1,max0(lycm,lys),-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            y0(3,j)=-ylwc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lys.lt.lycm)then
        if(lylw.ge.lycm)then
c
c         interface conditions: liquid to solid
c
          do j=1,2
            ylw(1,j)=ylwc(1,j)
            ylw(2,j)=ylwc(2,j)
            ylw(3,j)=(0.d0,0.d0)
            ylw(4,j)=(0.d0,0.d0)
            ylw(5,j)=ylwc(3,j)
            ylw(6,j)=ylwc(4,j)
          enddo
          ylw(1,3)=(0.d0,0.d0)
          ylw(2,3)=(0.d0,0.d0)
          ylw(3,3)=(1.d0,0.d0)
          ylw(4,3)=(0.d0,0.d0)
          ylw(5,3)=(0.d0,0.d0)
          ylw(6,3)=(0.d0,0.d0)
c
          if(lyr.eq.lycm.and.lycm.ge.lys)then
            call cmemcpy(ylw,y0,18)
            sety0lw=.true.
          endif
        else if(lylw.ge.lyob)then
c
c         lowest layer is within the mantle
c
          if(lylw.eq.lylwb)then
            call qpstart6a(ldeg,lylw,1,ylw)
          else
            do j=1,3
              do i=1,6
                ylw(i,j)=mas6x6up(i,2*j-1,lylw)
              enddo
            enddo
          endif
          if(lyr.eq.lylw.and.lylw.ge.lys)then
            call cmemcpy(ylw,y0,18)
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylw,lycm)-1,max0(lyob,lys),-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(lyr.eq.ly)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lys.lt.lyob)then
        if(lylw.ge.lyob)then
c
c         interface conditions: solid to liquid
c
          y4max=cdabs(ylw(4,3))
          j0=3
          do j=1,2
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
            ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
            ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
            ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
          enddo
          if(sety0lw)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        else
c
c         lowest layer is within the atmosphere
c
          if(lylw.eq.lylwb)then
            call qpstart4a(ldeg,lylw,1,ylwc)
          else
            do j=1,2
              do i=1,4
                ylwc(i,j)=mas4x4up(i,2*j-1,lylw)
              enddo
            enddo
          endif
c
          if(lylw.eq.lyr.and.lylw.ge.lys)then
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=ylwc(2,j)
              y0(3,j)=-ylwc(2,j)/(cro(lylw)*comi2*crrup(lylw)**2)
              y0(4,j)=(0.d0,0.d0)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylw,lyob)-1,lys,-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            y0(3,j)=-ylwc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0lw=.true.
        endif
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=ylwc(2,j)
          ylw(3,j)=-ylwc(2,j)/(cro(lys)*comi2*crrup(lys)**2)
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        do i=1,6
          ylw(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
c
      if(sety0up.and..not.sety0lw.or..not.sety0up.and.sety0lw)then
        do j=1,3
          y0(1,j)=y0(1,j)/crrup(lyr)
          y0(2,j)=y0(2,j)/crrup(lyr)**2
          y0(3,j)=y0(3,j)/crrup(lyr)
          y0(4,j)=y0(4,j)/crrup(lyr)**2
c          y0(5,j)=y0(5,j)
          y0(6,j)=y0(6,j)/crrup(lyr)
        enddo
      else
        stop 'Error #1 in qpsprop!'
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      if(vsup(lys).le.0.d0)then
        do istp=1,2
          do i=1,4
            b4(i,istp)=(0.d0,0.d0)
          enddo
          b4(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,2
          do i=1,2
            coef4(i,j)=yup(i,j)
            coef4(i,j+2)=-ylw(i,j)
          enddo
          do i=3,4
            coef4(i,j)=yup(i+2,j)
            coef4(i,j+2)=-ylw(i+2,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef4,b4,4,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.lt.lys)then
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
        do istp=3,4
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            b6(i,istp)=(0.d0,0.d0)
          enddo
          b6(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup(i,j)
            coef6(i,j+3)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.lt.lys)then
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j+3,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
      if(lylwa.le.0)return
c
c===============================================================================
c
      sety0lw=.false.
c
c     propagation within inner core
c
      if(lylwa.ge.lycc)then
c
c       lowest layer is within inner core
c
        lyswap=lylwa
        call qpstart6a(ldeg,lyswap,1,ylw)
        if(lylwa.eq.lyr.and.lylwa.ge.lys)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      endif
c
      do ly=lylwa-1,max0(lycc,lys),-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(ly.eq.lyr)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation within outer core
c
      if(lys.lt.lycc)then
        if(lylwa.ge.lycc)then
c
c         interface conditions: solid to liquid
c
          y4max=cdabs(ylw(4,3))
          j0=3
          do j=1,2
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
            ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
            ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
            ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
          enddo
          if(sety0lw)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        else if(lylwa.ge.lycm)then
c
c         lowest layer is within the liquid core
c
          lyswap=lylwa
          call qpstart4a(ldeg,lyswap,1,ylwc)
c
          if(lylwa.eq.lyr.and.lylwa.ge.lys)then
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=ylwc(2,j)
              y0(3,j)=-ylwc(2,j)/(cro(lylwa)*comi2*crrup(lylwa)**2)
              y0(4,j)=(0.d0,0.d0)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylwa,lycc)-1,max0(lycm,lys),-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            y0(3,j)=-ylwc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lys.lt.lycm)then
        if(lylwa.ge.lycm)then
c
c         interface conditions: liquid to solid
c
          do j=1,2
            ylw(1,j)=ylwc(1,j)
            ylw(2,j)=ylwc(2,j)
            ylw(3,j)=(0.d0,0.d0)
            ylw(4,j)=(0.d0,0.d0)
            ylw(5,j)=ylwc(3,j)
            ylw(6,j)=ylwc(4,j)
          enddo
          ylw(1,3)=(0.d0,0.d0)
          ylw(2,3)=(0.d0,0.d0)
          ylw(3,3)=(1.d0,0.d0)
          ylw(4,3)=(0.d0,0.d0)
          ylw(5,3)=(0.d0,0.d0)
          ylw(6,3)=(0.d0,0.d0)
c
          if(lycm.eq.lyr.and.lycm.ge.lys)then
            call cmemcpy(ylw,y0,18)
            sety0lw=.true.
          endif
        else if(lylwa.ge.lyob)then
c
c         lowest layer is within the mantle
c
          lyswap=lylwa
          call qpstart6a(ldeg,lyswap,1,ylw)
          if(lylwa.eq.lyr.and.lylwa.ge.lys)then
            call cmemcpy(ylw,y0,18)
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylwa,lycm)-1,max0(lyob,lys),-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(ly.eq.lyr)then
          call cmemcpy(ylw,y0,18)
          sety0lw=.true.
        endif
      enddo
c
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lys.lt.lyob)then
        if(lylwa.ge.lyob)then
c
c         interface conditions: solid to liquid
c
          y4max=cdabs(ylw(4,3))
          j0=3
          do j=1,2
            if(y4max.lt.cdabs(ylw(4,j)))then
              y4max=cdabs(ylw(4,j))
              j0=j
            endif
          enddo
          do i=1,6
            cyswap=ylw(i,j0)
            ylw(i,j0)=ylw(i,3)
            ylw(i,3)=cyswap
          enddo
          do j=1,2
            ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
            ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
            ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
            ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
          enddo
          if(sety0lw)then
            do i=1,6
              cyswap=y0(i,j0)
              y0(i,j0)=y0(i,3)
              y0(i,3)=cyswap
            enddo
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
              enddo
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        else
c
c         lowest layer is within the atmosphere
c
          lyswap=lylwa
          call qpstart4a(ldeg,lyswap,1,ylwc)
c
          if(lylwa.eq.lyr.and.lylwa.ge.lys)then
            do j=1,2
              y0(1,j)=ylwc(1,j)
              y0(2,j)=ylwc(2,j)
              y0(3,j)=-ylwc(2,j)/(cro(lylwa)*comi2*crrup(lylwa)**2)
              y0(4,j)=(0.d0,0.d0)
              y0(5,j)=ylwc(3,j)
              y0(6,j)=ylwc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
            sety0lw=.true.
          endif
        endif
      endif
c
      do ly=min0(lylwa,lyob)-1,lys,-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(sety0lw)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            y0(3,j)=-ylwc(2,j)/(cro(ly)*comi2*crrup(ly)**2)
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
          sety0lw=.true.
        endif
      enddo
c
      if(lys.lt.lyob.or.
     &   lys.ge.lycm.and.lys.lt.lycc)then
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=ylwc(2,j)
          ylw(3,j)=-ylwc(2,j)/(cro(lys)*comi2*crrup(lys)**2)
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        do i=1,6
          ylw(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
c
      if(sety0lw)then
        do j=1,3
          y0(1,j)=y0(1,j)/crrup(lyr)
          y0(2,j)=y0(2,j)/crrup(lyr)**2
          y0(3,j)=y0(3,j)/crrup(lyr)
          y0(4,j)=y0(4,j)/crrup(lyr)**2
c          y0(5,j)=y0(5,j)
          y0(6,j)=y0(6,j)/crrup(lyr)
        enddo
      else if(.not.sety0up)then
        stop 'Error #2 in qpsprop!'
      endif
c
c===============================================================================
c     source function
c===============================================================================
c
      if(vsup(lys).le.0.d0)then
        do istp=1,2
          do i=1,4
            b4(i,istp)=(0.d0,0.d0)
          enddo
          b4(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,2
          do i=1,2
            coef4(i,j)=yup(i,j)
            coef4(i,j+2)=-ylw(i,j)
          enddo
          do i=3,4
            coef4(i,j)=yup(i+2,j)
            coef4(i,j+2)=-ylw(i+2,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef4,b4,4,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.lt.lys)then
          do istp=1,2
            do i=1,6
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)-b4(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,2
            do i=1,6
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)-b4(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      else
        do istp=1,4
          do i=1,6
            b6(i,istp)=(0.d0,0.d0)
          enddo
          b6(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup(i,j)
            coef6(i,j+3)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdsvd500(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop: anormal exit from cdsvd500!'
          return
        endif
        if(lyr.lt.lys)then
          do istp=1,4
            do i=1,6
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)-b6(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,6
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)-b6(j+3,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
      return
      end