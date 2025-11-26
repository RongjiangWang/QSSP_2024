      subroutine qpgetinp(unit)
      use qpalloc
      implicit none
      integer*4 unit
c
c     work space
c
      integer*4 i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo,sdfsel,ierr
      integer*4 ipath,isurf
      real*8 twindow,twinout,suppress,munit,sfe,sfn,sfz,omi
      real*8 strike,dip,rake,depdif,dswap(14)
      character*80 grndir,outfile,fswap
      logical*2 lswap
c
c     uniform receiver depth
c     ======================
c
      lyadd=1
c
      call skipdoc(unit)
      read(unit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skipdoc(unit)
      read(unit,*)twindow,dtout
      ntcut=1+idnint(twindow/dtout)
      nt=2
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      dt=twindow/dble(nt-1)
      nf=nt/2
      df=1.d0/(dble(nt)*dt)
c
      call skipdoc(unit)
      read(unit,*)fcut
      nfcut=min0(nf,1+idnint(fcut/df))
      fcut=dble(nfcut-1)*df
      call skipdoc(unit)
      read(unit,*)slwmax
      if(slwmax.le.0.d0)then
        stop ' Error in qpgetinp: bad selection of max. slowness!'
      else
        slwmax=slwmax/KM2M
      endif
c
      call skipdoc(unit)
      read(unit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        suppress=dexp(-1.d0)
      endif
      fi=dlog(suppress)*df/PI2
      call skipdoc(unit)
      read(unit,*)ipath,minpath,maxpath
      minpath=minpath*KM2M
      maxpath=maxpath*KM2M
      call skipdoc(unit)
      read(unit,*)rearth,isurf
      rearth=rearth*KM2M
      freesurf=isurf.eq.1
      if(ipath.ne.1)then
        ipatha=0
        ipathb=0
        minpath=0.d0
        maxpath=rearth
      else
        if(minpath.lt.dpr.or.minpath.ge.dmax1(rearth,maxpath))then
          ipatha=0
          minpath=0.d0
        else
          ipatha=1
          lyadd=lyadd+1
        endif
        if(maxpath.ge.rearth.or.maxpath.le.minpath)then
          ipathb=0
          maxpath=rearth
        else
          ipathb=1
          lyadd=lyadd+1
        endif
      endif
c
c     cutoffs of spectra
c     ==================
c
      call skipdoc(unit)
      read(unit,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop ' Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      call skipdoc(unit)
      read(unit,*)i,j,ldegmin,ldegcut
      selpsv=i.eq.1
      selsh=j.eq.1
      if(.not.(selpsv.or.selsh))then
        stop ' Error in qpgetinp: none of PSV and SH is selected!'
      endif
      if(ldegcut.lt.ldegmin)ldegcut=ldegmin
      ldegmin=min0(max0(1+ndmax,ldegmin),ldegcut)
      ldegmax=ldegcut+1+ndmax
c
      omi=PI*fcut
      if(idint(rearth*omi*slwmax).gt.ldegmax)then
        write(*,'(a)')' Warning from qpgrnspec:'
        write(*,'(a)')' Maximum harmonic degree not large enough!'
        write(*,'(a,i6)')' Suggestion: increase the max. degree to >= ',
     &                   1+idint(rearth*omi*slwmax)
        write(*,'(a,f10.5,a)')' ... or decrease frequency cutoff to <=',
     &                   dble(ldegcut)/(rearth*PI2*slwmax),' Hz'
        write(*,'(a,f10.5,a)')' ... or decrease slowness cutoff to <=',
     &                   dble(ldegcut)*KM2M/(rearth*omi),' s/km'
        stop
      endif
c
      allocate(plm(0:ldegmax,0:2),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: plm not allocated!'
c
      allocate(ul0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ul0 not allocated!'
      allocate(vl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: vl0 not allocated!'
      allocate(wl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: wl0 not allocated!'
      allocate(el0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: el0 not allocated!'
      allocate(fl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: fl0 not allocated!'
      allocate(gl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gl0 not allocated!'
c
      allocate(pl0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: pl0 not allocated!'
      allocate(ql0(0:ldegmax,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ql0 not allocated!'
c
      allocate(urlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: urlm not allocated!'
      allocate(utlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: utlm not allocated!'
      allocate(uplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: uplm not allocated!'
c
      allocate(grlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: grlm not allocated!'
      allocate(gtlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gtlm not allocated!'
      allocate(gplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: gplm not allocated!'
c
      allocate(errlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: errlm not allocated!'
      allocate(ertlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ertlm not allocated!'
      allocate(erplm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: erplm not allocated!'
      allocate(etrlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etrlm not allocated!'
      allocate(ett0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ett0lm not allocated!'
      allocate(ettalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ettalm not allocated!'
      allocate(ettblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ettblm not allocated!'
      allocate(etp0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etp0lm not allocated!'
      allocate(etpalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etpalm not allocated!'
      allocate(etpblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: etpblm not allocated!'
      allocate(eprlm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eprlm not allocated!'
      allocate(ept0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: ept0lm not allocated!'
      allocate(eptalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eptalm not allocated!'
      allocate(eptblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eptblm not allocated!'
      allocate(epp0lm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: epp0lm not allocated!'
      allocate(eppalm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eppalm not allocated!'
      allocate(eppblm(0:ldegmax,6,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgetinp: eppblm not allocated!'
c
      allocate(lyupp(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyupp not allocated!'
      allocate(lyups(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyups not allocated!'
      allocate(lyupt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lyupt not allocated!'
      allocate(lylwp(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylwp not allocated!'
      allocate(lylws(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylws not allocated!'
      allocate(lylwt(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lylwt not allocated!'
c
c     Green's function files
c     ======================
c
      call skipdoc(unit)
      read(unit,*)ngrn,grndir
      if(ngrn.le.0)then
        stop ' bad number of source depths!'
      endif
      lyadd=lyadd+ngrn
c
      allocate(grndep(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: grndep not allocated!'
      allocate(grnsel(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: grnsel not allocated!'
      allocate(lygrn(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lygrn not allocated!'
c
      allocate(specfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: specfile not allocated!'
      allocate(uspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: uspecfile not allocated!'
      allocate(vspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vspecfile not allocated!'
      allocate(wspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: wspecfile not allocated!'
      allocate(especfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: especfile not allocated!'
      allocate(fspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: fspecfile not allocated!'
      allocate(gspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: gspecfile not allocated!'
      allocate(pspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: pspecfile not allocated!'
      allocate(qspecfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qspecfile not allocated!'
c
      allocate(isg1(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: isg1 not allocated!'
      allocate(isg2(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: isg2 not allocated!'
      allocate(nsg(ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: nsg not allocated!'
      allocate(grnrr0(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: grnrr0 not allocated!'
      allocate(disk(0:ldegmax),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: disk not allocated!'
c
      do ig=1,ngrn
        call skipdoc(unit)
        read(unit,*)grndep(ig),grnrr0(ig),specfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' bad Green function selection!'
        endif
        grndep(ig)=grndep(ig)*KM2M
        grnrr0(ig)=grnrr0(ig)*KM2M
      enddo
c
c     sort green function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            dswap(1)=grndep(i)
            dswap(2)=grnrr0(i)
            fswap=specfile(i)
            iswap=grnsel(i)
c
            grndep(i)=grndep(j)
            grnrr0(i)=grnrr0(j)
            specfile(i)=specfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            grnrr0(j)=dswap(2)
            specfile(j)=fswap
            grnsel(j)=iswap
          endif
        enddo
      enddo
c
      do flen=80,1,-1
        if(grndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do ig=1,ngrn
        uspecfile(ig)=grndir(1:flen)//'U_'//specfile(ig)
        vspecfile(ig)=grndir(1:flen)//'V_'//specfile(ig)
        wspecfile(ig)=grndir(1:flen)//'W_'//specfile(ig)
        especfile(ig)=grndir(1:flen)//'E_'//specfile(ig)
        fspecfile(ig)=grndir(1:flen)//'F_'//specfile(ig)
        gspecfile(ig)=grndir(1:flen)//'G_'//specfile(ig)
        pspecfile(ig)=grndir(1:flen)//'P_'//specfile(ig)
        qspecfile(ig)=grndir(1:flen)//'Q_'//specfile(ig)
      enddo
c
c     multi-event source parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)neq,sdfsel
c
      if(neq.le.0)goto 201
c
      allocate(mrreq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrreq not allocated!'
      allocate(mtteq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtteq not allocated!'
      allocate(mppeq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mppeq not allocated!'
      allocate(mrteq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrteq not allocated!'
      allocate(mpreq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mpreq not allocated!'
      allocate(mtpeq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtpeq not allocated!'
      allocate(latseq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: latseq not allocated!'
      allocate(lonseq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lonseq not allocated!'
      allocate(depseq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: depseq not allocated!'
      allocate(togseq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: togseq not allocated!'
      allocate(trsseq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: trsseq not allocated!'
      allocate(istfeq(neq),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: istfeq not allocated!'
c
      if(sdfsel.eq.1)then
        do is=1,neq
          call skipdoc(unit)
c
c         the six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
c
          read(unit,*)munit,mrreq(is),mtteq(is),mppeq(is),
     &                      mrteq(is),mpreq(is),mtpeq(is),
     &                      latseq(is),lonseq(is),depseq(is),
     &                      togseq(is),trsseq(is),istfeq(is)
          mtteq(is)=mtteq(is)*munit
          mppeq(is)=mppeq(is)*munit
          mrreq(is)=mrreq(is)*munit
          mtpeq(is)=mtpeq(is)*munit
          mpreq(is)=mpreq(is)*munit
          mrteq(is)=mrteq(is)*munit
          depseq(is)=depseq(is)*KM2M
          if(istfeq(is).lt.0)then
            stop ' bad selection for earthquake stf wavelet!'
          endif
        enddo
      else if(sdfsel.eq.2)then
        do is=1,neq
          call skipdoc(unit)
          read(unit,*)munit,strike,dip,rake,
     &                latseq(is),lonseq(is),depseq(is),
     &                togseq(is),trsseq(is),istfeq(is)
          call moments(munit,strike,dip,rake,
     &                 mtteq(is),mppeq(is),mrreq(is),
     &                 mtpeq(is),mpreq(is),mrteq(is))
          depseq(is)=depseq(is)*KM2M
          if(istfeq(is).lt.0)then
            stop ' bad selection for earthquake stf wavelet!'
          endif
        enddo
      else
        stop ' bad selection for dislocation source data format!'
      endif
c
201   continue
c
      call skipdoc(unit)
      read(unit,*)nsf
c
      allocate(sfrsf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfrsf not allocated!'
      allocate(sftsf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sftsf not allocated!'
      allocate(sfpsf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfpsf not allocated!'
      allocate(latssf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: latssf not allocated!'
      allocate(lonssf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lonssf not allocated!'
      allocate(depssf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: depssf not allocated!'
      allocate(togssf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: togssf not allocated!'
      allocate(trsssf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: trsssf not allocated!'
      allocate(istfsf(nsf),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: istfsf not allocated!'
c
      do is=1,nsf
        call skipdoc(unit)
        read(unit,*)munit,sfe,sfn,sfz,
     &              latssf(is),lonssf(is),depssf(is),
     &              togssf(is),trsssf(is),istfsf(is)
        sfpsf(is)=sfe*munit
        sftsf(is)=-sfn*munit
        sfrsf(is)=sfz*munit
        depssf(is)=depssf(is)*KM2M
        if(istfsf(is).lt.0)then
          stop ' bad selection for single force stf wavelet!'
        endif
      enddo
c
c     summarize all sources together
c
      ns=nsf+neq
c
      allocate(sfr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfr not allocated!'
      allocate(sft(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sft not allocated!'
      allocate(sfp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: sfp not allocated!'
c
      allocate(mrr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrr not allocated!'
      allocate(mtt(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtt not allocated!'
      allocate(mpp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mpp not allocated!'
      allocate(mrt(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mrt not allocated!'
      allocate(mpr(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mpr not allocated!'
      allocate(mtp(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: mtp not allocated!'
      allocate(lats(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lats not allocated!'
      allocate(lons(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lons not allocated!'
      allocate(deps(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: deps not allocated!'
      allocate(togs(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: togs not allocated!'
      allocate(trss(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: trss not allocated!'
      allocate(istf(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: istf not allocated!'
      allocate(eqs(ns),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: eqs not allocated!'
c
      do is=1,neq
        mtt(is)=mtteq(is)
        mpp(is)=mppeq(is)
        mrr(is)=mrreq(is)
        mtp(is)=mtpeq(is)
        mpr(is)=mpreq(is)
        mrt(is)=mrteq(is)
        sfr(is)=0.d0
        sft(is)=0.d0
        sfp(is)=0.d0
        lats(is)=latseq(is)
        lons(is)=lonseq(is)
        deps(is)=depseq(is)
        togs(is)=togseq(is)
        trss(is)=trsseq(is)
        istf(is)=istfeq(is)
        eqs(is)=.true.
      enddo
      do is=neq+1,ns
        mtt(is)=0.d0
        mpp(is)=0.d0
        mrr(is)=0.d0
        mtp(is)=0.d0
        mpr(is)=0.d0
        mrt(is)=0.d0
        sfr(is)=sfrsf(is-neq)
        sft(is)=sftsf(is-neq)
        sfp(is)=sfpsf(is-neq)
        lats(is)=latssf(is-neq)
        lons(is)=lonssf(is-neq)
        deps(is)=depssf(is-neq)
        togs(is)=togssf(is-neq)
        trss(is)=trsssf(is-neq)
        istf(is)=istfsf(is-neq)
        eqs(is)=.false.
      enddo
c
      togsmin=togs(1)
      do is=1,ns
        togsmin=dmin1(togsmin,togs(is))
      enddo
      do is=1,ns
        togs(is)=togs(is)-togsmin
      enddo
c
c     sort sub-events by depth
c
      do i=1,ns
        do j=i+1,ns
          if(deps(j).lt.deps(i))then
            dswap(1)=lats(i)
            dswap(2)=lons(i)
            dswap(3)=deps(i)
            dswap(4)=mtt(i)
            dswap(5)=mpp(i)
            dswap(6)=mrr(i)
            dswap(7)=mtp(i)
            dswap(8)=mpr(i)
            dswap(9)=mrt(i)
            dswap(10)=sft(i)
            dswap(11)=sfp(i)
            dswap(12)=sfr(i)
            dswap(13)=togs(i)
            dswap(14)=trss(i)
            iswap=istf(i)
            lswap=eqs(i)
c
            lats(i)=lats(j)
            lons(i)=lons(j)
            deps(i)=deps(j)
            mtt(i)=mtt(j)
            mpp(i)=mpp(j)
            mrr(i)=mrr(j)
            mtp(i)=mtp(j)
            mpr(i)=mpr(j)
            mrt(i)=mrt(j)
            sft(i)=sft(j)
            sfp(i)=sfp(j)
            sfr(i)=sfr(j)
            togs(i)=togs(j)
            trss(i)=trss(j)
            istf(i)=istf(j)
            eqs(i)=eqs(j)
c
            lats(j)=dswap(1)
            lons(j)=dswap(2)
            deps(j)=dswap(3)
            mtt(j)=dswap(4)
            mpp(j)=dswap(5)
            mrr(j)=dswap(6)
            mtp(j)=dswap(7)
            mpr(j)=dswap(8)
            mrt(j)=dswap(9)
            sft(j)=dswap(10)
            sfp(j)=dswap(11)
            sfr(j)=dswap(12)
            togs(j)=dswap(13)
            trss(j)=dswap(14)
            istf(j)=iswap
            eqs(j)=lswap
          endif
        enddo
      enddo
c
      isg1(1)=1
      is1=1
      do ig=1,ngrn-1
        depdif=0.5d0*(grndep(ig)+grndep(ig+1))
        isg2(ig)=isg1(ig)-1
        do is=is1,ns
          if(deps(is).lt.depdif)then
            isg2(ig)=is
          endif
        enddo
        isg1(ig+1)=isg2(ig)+1
        is1=isg1(ig+1)
      enddo
      isg2(ngrn)=ns
c
      do ig=1,ngrn
        nsg(ig)=max0(0,1+isg2(ig)-isg1(ig))
      enddo
c
c     receiver parameters
c     ===================
c
      call skipdoc(unit)
      read(unit,*)(icmp(i),i=1,11)
      call skipdoc(unit)
      read(unit,*)outfile
      call skipdoc(unit)
      read(unit,*)twinout
      ntcutout=1+idnint(twinout/dtout)
      call skipdoc(unit)
      read(unit,*)nhpf,f1corner,ilph
      if(nhpf.gt.0.and.(f1corner.le.0.d0.or.ilph.lt.1.or.ilph.gt.2))then
        stop ' bad highpass filter parameters!'
      endif
      call skipdoc(unit)
      read(unit,*)nlpf,f2corner,ihph
      if(nlpf.gt.0.and.(f2corner.le.0.d0.or.ilph.lt.1.or.ilph.gt.2))then
        stop ' bad lowpass filter parameters!'
      endif
      if(nhpf.gt.0.and.nlpf.gt.0.and.f1corner.ge.f2corner)then
        stop ' bad filter parameters!'
      endif
      call skipdoc(unit)
      read(unit,*)slwlwcut,slwupcut
      slwlwcut=slwlwcut/KM2M
      slwupcut=slwupcut/KM2M
      if(slwupcut.gt.slwmax.or.slwlwcut.ge.slwupcut)then
        slwlwcut=0.d0
        slwupcut=slwmax
      endif
      call skipdoc(unit)
      read(unit,*)nr
c
      allocate(latr(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: latr not allocated!'
      allocate(lonr(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: lonr not allocated!'
      allocate(tred(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: tred not allocated!'
      allocate(rname(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: rname not allocated!'
c
      do ir=1,nr
        call skipdoc(unit)
        read(unit,*)latr(ir),lonr(ir),rname(ir),tred(ir)
      enddo
c
      do flen=80,1,-1
        if(outfile(flen:flen).ne.' ')goto 300
      enddo
300   continue
      dispout(1)=outfile(1:flen)//'_disp_e.dat'
      dispout(2)=outfile(1:flen)//'_disp_n.dat'
      dispout(3)=outfile(1:flen)//'_disp_z.dat'
      veloout(1)=outfile(1:flen)//'_velo_e.dat'
      veloout(2)=outfile(1:flen)//'_velo_n.dat'
      veloout(3)=outfile(1:flen)//'_velo_z.dat'
      acceout(1)=outfile(1:flen)//'_acce_e.dat'
      acceout(2)=outfile(1:flen)//'_acce_n.dat'
      acceout(3)=outfile(1:flen)//'_acce_z.dat'
c
      rotaout(1)=outfile(1:flen)//'_rota_e.dat'
      rotaout(2)=outfile(1:flen)//'_rota_n.dat'
      rotaout(3)=outfile(1:flen)//'_rota_z.dat'
      rotarateout(1)=outfile(1:flen)//'_rota_rate_e.dat'
      rotarateout(2)=outfile(1:flen)//'_rota_rate_n.dat'
      rotarateout(3)=outfile(1:flen)//'_rota_rate_z.dat'
c
      strainout(1)=outfile(1:flen)//'_strain_ee.dat'
      strainout(2)=outfile(1:flen)//'_strain_en.dat'
      strainout(3)=outfile(1:flen)//'_strain_ez.dat'
      strainout(4)=outfile(1:flen)//'_strain_nn.dat'
      strainout(5)=outfile(1:flen)//'_strain_nz.dat'
      strainout(6)=outfile(1:flen)//'_strain_zz.dat'
      strainrateout(1)=outfile(1:flen)//'_strain_rate_ee.dat'
      strainrateout(2)=outfile(1:flen)//'_strain_rate_en.dat'
      strainrateout(3)=outfile(1:flen)//'_strain_rate_ez.dat'
      strainrateout(4)=outfile(1:flen)//'_strain_rate_nn.dat'
      strainrateout(5)=outfile(1:flen)//'_strain_rate_nz.dat'
      strainrateout(6)=outfile(1:flen)//'_strain_rate_zz.dat'
c
      stressout(1)=outfile(1:flen)//'_stress_ee.dat'
      stressout(2)=outfile(1:flen)//'_stress_en.dat'
      stressout(3)=outfile(1:flen)//'_stress_ez.dat'
      stressout(4)=outfile(1:flen)//'_stress_nn.dat'
      stressout(5)=outfile(1:flen)//'_stress_nz.dat'
      stressout(6)=outfile(1:flen)//'_stress_zz.dat'
      stressrateout(1)=outfile(1:flen)//'_stress_rate_ee.dat'
      stressrateout(2)=outfile(1:flen)//'_stress_rate_en.dat'
      stressrateout(3)=outfile(1:flen)//'_stress_rate_ez.dat'
      stressrateout(4)=outfile(1:flen)//'_stress_rate_nn.dat'
      stressrateout(5)=outfile(1:flen)//'_stress_rate_nz.dat'
      stressrateout(6)=outfile(1:flen)//'_stress_rate_zz.dat'
c
      gravout(1)=outfile(1:flen)//'_gravitation_e.dat'
      gravout(2)=outfile(1:flen)//'_gravitation_n.dat'
      gravout(3)=outfile(1:flen)//'_gravitation_z.dat'
c
      grmout=outfile(1:flen)//'_gravimeter.dat'
c
c     multilayered model parameters
c     =============================
c
      call skipdoc(unit)
      read(unit,*)l,i,j
      dispersion=i.eq.1
      isotherm=j.eq.1
c
      l0=l+1
      allocate(dp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0 not allocated!'
      allocate(vp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0 not allocated!'
      allocate(vs0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0 not allocated!'
      allocate(ro0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0 not allocated!'
      allocate(qp0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0 not allocated!'
      allocate(qs0(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0 not allocated!'
c
      allocate(dp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0up not allocated!'
      allocate(vp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0up not allocated!'
      allocate(vs0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0up not allocated!'
      allocate(ro0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0up not allocated!'
      allocate(qp0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0up not allocated!'
      allocate(qs0up(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0up not allocated!'
c
      allocate(dp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: dp0lw not allocated!'
      allocate(vp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vp0lw not allocated!'
      allocate(vs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: vs0lw not allocated!'
      allocate(ro0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: ro0lw not allocated!'
      allocate(qp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qp0lw not allocated!'
      allocate(qs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpgettinp: qs0lw not allocated!'
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=rearth
      do i=1,l
        call skipdoc(unit)
        read(unit,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(vp0(i).le.0.d0.or.vs0(i).lt.0.d0.or.ro0(i).le.0.d0.or.
     &     qp0(i).le.0.d0.or.qs0(i).lt.0.d0)then
          stop ' Error in qpgetinp: bad seismic parameter!'
        endif
c
        if(i.eq.1)then
          if(dp0(1).gt.0.d0)then
            stop ' Error in qpgetinp: bad start depth!'
          endif
          depatmos=-KM2M*dp0(1)
          rratmos=rearth+depatmos
        endif
c
        dp0(i)=KM2M*dp0(i)+depatmos
        vp0(i)=KM2M*vp0(i)
        vs0(i)=KM2M*vs0(i)
        ro0(i)=KM2M*ro0(i)
        if(vs0(i).gt.0.d0)qsmin=dmin1(qsmin,qs0(i))
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop ' Error in qpgetinp: bad layering of earth model!'
          endif
        endif
      enddo
c      open(39,file='model.out',status='unknown')
c      write(39,'(a)')'  no   depth[km]    vp[km/s]    vs[km/s]'
c     &                 //'  ro[g/cm^3]          qp          qs'
c      do i=1,l
c        write(39,'(i4,4f12.3,2f12.1)')i,dp0(i)/KM2M,vp0(i)/KM2M,
c     &          vs0(i)/KM2M,ro0(i)/KM2M,qp0(i),qs0(i)
c      enddo
c      close(39)
c
      if(dp0(l).gt.rratmos)then
        stop ' Error in qpgetinp: earth radius larger than pre-defined!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        dp0(l)=rratmos
        vp0(l)=vp0(l-1)
        vs0(l)=vs0(l-1)
        ro0(l)=ro0(l-1)
        qp0(l)=qp0(l-1)
        qs0(l)=qs0(l-1)
      endif
c
      dpr=dpr+depatmos
      if(dpr.lt.0.d0.or.dpr.gt.dp0(l))then
        stop ' Error in qpgetinp: receiver too shallow or too deep!'
      endif
      do i=1,ngrn
        grndep(i)=grndep(i)+depatmos
        if(grndep(i).lt.0.d0.or.grndep(i).ge.dp0(l))then
          stop ' Error in qpgetinp: source too shallow or too deep!'
        endif
      enddo
      do is=1,ns
        deps(is)=deps(is)+depatmos
      enddo
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then

          l0=l0+1
          dp0up(l0)=dp0(i-1)
          vp0up(l0)=vp0(i-1)
          vs0up(l0)=vs0(i-1)
          ro0up(l0)=ro0(i-1)
          qp0up(l0)=qp0(i-1)
          qs0up(l0)=qs0(i-1)
c
          dp0lw(l0)=dp0(i)
          vp0lw(l0)=vp0(i)
          vs0lw(l0)=vs0(i)
          ro0lw(l0)=ro0(i)
          qp0lw(l0)=qp0(i)
          qs0lw(l0)=qs0(i)
        endif
      enddo
c
c     end of inputs
c     =============
c
      return
      end
