      subroutine qpqmodel(f)
      use qpalloc
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      real*8 f
c
      integer*4 ly
      real*8 fac,pup,plw,sup,slw,qswap,zero
      complex*16 cqpup,cqsup,cqplw,cqslw
c
      if(.not.dispersion.or.f.ge.FSBREF)then
        fac=0.d0
      else if(f.le.FSBLW)then
        fac=dlog(FSBLW/FSBREF)/PI
      else
        fac=dlog(f/FSBREF)/PI
      endif
c
      if(f.gt.0.d0)then
        zero=1.d0
      else
        zero=0.d0
      endif
c
      if(1.d0+fac/qsmin.le.0.d0)then
        stop 'Error in qpqmodel: too small Qs value!'
      endif
      do ly=1,ly0
        cqpup=dcmplx(1.d0,0.5d0*zero/qpup(ly))
        cqplw=dcmplx(1.d0,0.5d0*zero/qplw(ly))
c
        if(ly.ge.lyob.and.ly.lt.lycm.or.ly.ge.lycc)then
          pup=1.d0+fac/qpup(ly)
          plw=1.d0+fac/qplw(ly)
        else
          pup=1.d0
          plw=1.d0
        endif
c
        cvpup(ly)=dcmplx(vpup(ly)*pup/cdabs(cqpup),0.d0)*cqpup
        cvplw(ly)=dcmplx(vplw(ly)*plw/cdabs(cqplw),0.d0)*cqplw
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
c
        if(ly.ge.lyob.and.ly.lt.lycm.or.ly.ge.lycc)then
          sup=1.d0+fac/qsup(ly)
          slw=1.d0+fac/qslw(ly)
          cqsup=dcmplx(1.d0,0.5d0*zero/qsup(ly))
          cqslw=dcmplx(1.d0,0.5d0*zero/qslw(ly))
          cvsup(ly)=dcmplx(vsup(ly)*sup/cdabs(cqsup),0.d0)*cqsup
          cvslw(ly)=dcmplx(vslw(ly)*slw/cdabs(cqslw),0.d0)*cqslw
          cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        else
          cvsup(ly)=(0.d0,0.d0)
          cvslw(ly)=(0.d0,0.d0)
          cvs(ly)=(0.d0,0.d0)
        endif
c
        cmuup(ly)=croup(ly)*cvsup(ly)**2
        claup(ly)=croup(ly)*cvpup(ly)**2-(2.d0,0.d0)*cmuup(ly)
        cmulw(ly)=crolw(ly)*cvslw(ly)**2
        clalw(ly)=crolw(ly)*cvplw(ly)**2-(2.d0,0.d0)*cmulw(ly)
        cmu(ly)=(0.5d0,0.d0)*(cmuup(ly)+cmulw(ly))
        cla(ly)=(0.5d0,0.d0)*(claup(ly)+clalw(ly))
      enddo
      return
      end
