      subroutine qpdifmatl(ly,ldeg,rr,mat)
      use qpalloc
      implicit none
c
      integer*4 ly,ldeg
      real*8 rr
      complex*16 mat(6,6)
c
c     4x4 coefficient matrix for spheroidal mode l > 0 in liquid layer
c
      real*8 f,flm,fn2,mass,rorr,beta,dro,ro1,n2up,n2lw
      real*8 vp2rr,vp2it,vp2ab
      complex*16 cldeg,clp1,cll1
      complex*16 cup,clw,crr,crorr,cgrrr,cvprr,cvp2rr
      complex*16 cgarr,cn2rr
c
      complex*16 c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      f=dreal(comi)/PI2
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
c
      crr=dcmplx(rr,0.d0)
      cup=(crr-crrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      cvprr=cup*cvpup(ly)+clw*cvplw(ly)
      cvp2rr=cvprr**2
      beta=dlog(roup(ly)/rolw(ly))/(rrup(ly)-rrlw(ly))
c
      if(rr.le.rrlw(ly))then
        mass=0.d0
        crorr=crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
      else if(ly.ge.lyos)then
        crorr=cup*croup(ly)+clw*crolw(ly)
        rorr=dreal(crorr)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        dro=(rorr-rolw(ly))/(rr-rrlw(ly))
        ro1=rorr-dro*rrlw(ly) 
        mass=PI*(rr-rrlw(ly))*((4.d0/3.d0)*ro1
     &        *(rr**2+rr*rrlw(ly)+rrlw(ly)**2)
     &        +dro*(rr**3+rr**2*rrlw(ly)
     &        +rr*rrlw(ly)**2+rrlw(ly)**3))
      else
        rorr=rolw(ly)*dexp(dlog(roup(ly)/rolw(ly))
     &                  *(rr-rrlw(ly))/(rrup(ly)-rrlw(ly)))
        crorr=dcmplx(rorr,0.d0)
        cgarr=dcmplx(2.d0*PI2*BIGG*rorr,0.d0)
c
        mass=2.d0*PI2/beta**3
     &      *(rorr*(beta*rr*(beta*rr-2.d0)+2.d0)
     &       -rolw(ly)*(beta*rrlw(ly)*(beta*rrlw(ly)-2.d0)+2.d0))
      endif
      cgrrr=cgrlw(ly)*(crrlw(ly)/crr)**2
     &     +dcmplx(BIGG*mass/rr**2,0.d0)
c
c     Brunt-Väisälä frequency
c
      if(isotherm.or.ly.ge.lyos)then
        cn2rr=(0.d0,0.d0)
      else
        vp2ab=cdabs(cvp2rr)
        vp2it=-cdabs(cgrrr)/beta
c
c       ab: adiabatic, it: isothermal
c
        n2up=-cdabs(cgrup(ly))*(beta+cdabs(cgrup(ly))/vpup(ly)**2)
        n2lw=-cdabs(cgrlw(ly))*(beta+cdabs(cgrlw(ly))/vplw(ly)**2)
        if(n2up.lt.0.d0.or.n2lw.lt.0.d0)then
          cvp2rr=dcmplx(vp2it/cdabs(cvp2rr),0.d0)*cvp2rr
          cn2rr=(0.d0,0.d0)
        else
          flm=0.25d0*fbvatm+vlamb*dsqrt(cdabs(cll1))/(PI2*rearth)
          fn2=fbvatm+0.25d0*vlamb*dsqrt(cdabs(cll1))/(PI2*rearth)
          if(f.ge.flm)then
            vp2rr=vp2ab
          else
            vp2rr=vp2it+(vp2ab-vp2it)*(dsin(0.5d0*PI*f/flm))**2
          endif
          cvp2rr=dcmplx(vp2rr/cdabs(cvp2rr),0.d0)*cvp2rr
          cn2rr=-cgrrr*(dcmplx(beta,0.d0)+cgrrr/cvp2rr)
          if(f.le.fn2)then
            cn2rr=cn2rr*dcmplx((dsin(0.5d0*PI*f/fn2))**3,0.d0)
          endif
        endif
      endif
c
      mat(1,1)=cgrrr/cvp2rr-c1/crr
      mat(1,2)=cll1/crr-comi2*crr/cvp2rr
      mat(1,3)=-crr/cvp2rr
      mat(1,4)=(0.d0,0.d0)
c
      mat(2,1)=(c1-cn2rr/comi2)/crr
      mat(2,2)=cn2rr/cgrrr
      mat(2,3)=cn2rr/(comi2*cgrrr)
      mat(2,4)=(0.d0,0.d0)
c
      mat(3,1)=cgarr/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=-clp1/crr
      mat(3,4)=c1/crr
c
      mat(4,1)=cgarr*clp1/crr
      mat(4,2)=-cgarr*cll1/crr
      mat(4,3)=(0.d0,0.d0)
      mat(4,4)=cldeg/crr
c
      return
      end