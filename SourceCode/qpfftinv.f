      subroutine qpfftinv(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
      integer*4 lf,mf,ir
      real*8 f,t
      complex*16 s
c
      real*8,allocatable:: y0(:),y(:),dswap(:)
      complex*16,allocatable:: cy(:,:),lpf(:),hpf(:),bpf(:),cswap(:)
c
      allocate(y0(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: y0 not allocated!'
      allocate(y(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: y not allocated!'
      allocate(cy(nf,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: cy not allocated!'
      allocate(bpf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: bpf not allocated!'
      allocate(lpf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: lpf not allocated!'
      allocate(hpf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: hpf not allocated!'
      allocate(dswap(4*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: dswap not allocated!'
      allocate(cswap(2*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpfftinv: cswap not allocated!'
c
      do lf=1,nf
        f=dble(lf-1)*df
        s=dcmplx(-PI2*fi,PI2*f)
        do ir=1,nr
          gm(lf,ir)=(-gz(lf,ir)+dcmplx(freeairgrd,0.d0)*uz(lf,ir))/s
     &             +uz(lf,ir)*s
c
c         note: s-domain divided by s => t-domain differeated over time
c
        enddo
      enddo
c
c     low-pass filter
c
      if(nlpf.gt.0.or.nhpf.gt.0)then
        if(nlpf.gt.0.and.f2corner.gt.0.d0)then
          call butterworth(nlpf,f2corner,df,nf,cswap)
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
          call four1w(cswap,dswap,2*nf,+1)
          do lf=1,2*nf
            t=dble(lf-1)*dt
            cswap(lf)=dcmplx(dswap(2*lf-1),0.d0)
     &               *dcmplx(dexp(PI2*fi*t)*df,0.d0)
          enddo
          call four1w(cswap,dswap,2*nf,-1)
          do lf=1,nf
            lpf(lf)=cswap(lf)*dcmplx(dt,0.d0)
          enddo
          if(ilph.eq.2)then
            do lf=1,nf
              lpf(lf)=lpf(lf)*dconjg(lpf(lf))
            enddo
          endif
        else
          do lf=1,nf
            lpf(lf)=(1.d0,0.d0)
          enddo
        endif
        if(nhpf.gt.0.and.f1corner.gt.0.d0)then
          call bandpass(nhpf,f1corner,f2corner,df,nf,cswap)
          mf=1
          do lf=2*nf,nf+2,-1
            mf=mf+1
            cswap(lf)=dconjg(cswap(mf))
          enddo
          cswap(nf+1)=(0.d0,0.d0)
          call four1w(cswap,dswap,2*nf,+1)
          do lf=1,2*nf
            t=dble(lf-1)*dt
            cswap(lf)=dcmplx(dswap(2*lf-1),0.d0)
     &               *dcmplx(dexp(PI2*fi*t)*df,0.d0)
          enddo
          call four1w(cswap,dswap,2*nf,-1)
          do lf=1,nf
            hpf(lf)=cswap(lf)*dcmplx(dt,0.d0)
          enddo
          if(ihph.eq.2)then
            do lf=1,nf
              hpf(lf)=hpf(lf)*dconjg(hpf(lf))
            enddo
          endif
        else
          do lf=1,nf
            hpf(lf)=(1.d0,0.d0)
          enddo
        endif
c
        do lf=1,nf
          bpf(lf)=lpf(lf)*hpf(lf)
        enddo
c
        do lf=1,nf
          do ir=1,nr
            ue(lf,ir)=ue(lf,ir)*bpf(lf)
            un(lf,ir)=un(lf,ir)*bpf(lf)
            uz(lf,ir)=uz(lf,ir)*bpf(lf)
c
            ge(lf,ir)=ge(lf,ir)*bpf(lf)
            gn(lf,ir)=gn(lf,ir)*bpf(lf)
            gz(lf,ir)=gz(lf,ir)*bpf(lf)
c
            roe(lf,ir)=roe(lf,ir)*bpf(lf)
            ron(lf,ir)=ron(lf,ir)*bpf(lf)
            roz(lf,ir)=roz(lf,ir)*bpf(lf)
c
            uee(lf,ir)=uee(lf,ir)*bpf(lf)
            uen(lf,ir)=uen(lf,ir)*bpf(lf)
            uez(lf,ir)=uez(lf,ir)*bpf(lf)
            unn(lf,ir)=unn(lf,ir)*bpf(lf)
            unz(lf,ir)=unz(lf,ir)*bpf(lf)
            uzz(lf,ir)=uzz(lf,ir)*bpf(lf)
c
            see(lf,ir)=see(lf,ir)*bpf(lf)
            sen(lf,ir)=sen(lf,ir)*bpf(lf)
            sez(lf,ir)=sez(lf,ir)*bpf(lf)
            snn(lf,ir)=snn(lf,ir)*bpf(lf)
            snz(lf,ir)=snz(lf,ir)*bpf(lf)
            szz(lf,ir)=szz(lf,ir)*bpf(lf)
c
            gm(lf,ir)=gm(lf,ir)*bpf(lf)
          enddo
        enddo
      endif
c
      if(icmp(1).eq.1)then
c
c       output seismograms: displacement vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ue,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,dispout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,un,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,dispout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,dispout(3))
      endif
      if(icmp(2).eq.1)then
c
c       output seismograms: velocity vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ue,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,veloout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,un,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,veloout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,veloout(3))
      endif
      if(icmp(3).eq.1)then
c
c       output seismograms: acceleration vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ue,cy,cswap,
     &                 dswap,2,ntcutout,rname,y0,y,acceout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,un,cy,cswap,
     &                 dswap,2,ntcutout,rname,y0,y,acceout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uz,cy,cswap,
     &                 dswap,2,ntcutout,rname,y0,y,acceout(3))
      endif
c
      if(icmp(4).eq.1)then
c
c       output seismograms: strain tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uee,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uen,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uez,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(3))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,unn,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(4))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,unz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(5))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uzz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,strainout(6))
      endif
c
      if(icmp(5).eq.1)then
c
c       output seismograms: strain rate tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uee,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uen,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uez,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(3))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,unn,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(4))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,unz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(5))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,uzz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,strainrateout(6))
      endif
      if(icmp(6).eq.1)then
c
c       output seismograms: stress tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,see,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,sen,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,sez,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(3))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,snn,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(4))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,snz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(5))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,szz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,stressout(6))
      endif
c
      if(icmp(7).eq.1)then
c
c       output seismograms: stress rate tensor
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,see,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,sen,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,sez,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(3))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,snn,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(4))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,snz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(5))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,szz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,stressrateout(6))
      endif
c
      if(icmp(8).eq.1)then
c
c       output seismograms: rotation vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,roe,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,rotaout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ron,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,rotaout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,roz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,rotaout(3))
      endif
c
      if(icmp(9).eq.1)then
c
c       output seismograms: rotation rate vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,roe,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,rotarateout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ron,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,rotarateout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,roz,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,rotarateout(3))
      endif
      if(icmp(10).eq.1)then
c
c       output seismograms: local gravitational acceleration vector
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,ge,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,gravout(1))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,gn,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,gravout(2))
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,gz,cy,cswap,
     &                 dswap,0,ntcutout,rname,y0,y,gravout(3))
      endif
c
      if(icmp(11).eq.1)then
c
c       output seismograms: gravimeter response
c
        call transfs2t(nf,nr,df,togsmin,dt,dtout,fi,tred,gm,cy,cswap,
     &                 dswap,1,ntcutout,rname,y0,y,grmout)
          
      endif
      return
      end
