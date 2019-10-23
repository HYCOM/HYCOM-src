      subroutine geopar
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
!
! --- set up model parameters related to geography
!
! --- hycom version 2.1
!
      implicit none
!
      real*4,  allocatable, dimension(:,:) :: g_sc
!
      real      dp0kf,dpm,dpms,ds0kf,dsm,dsms
      real      hmina,hminb,hmaxa,hmaxb
      real*8    sum_ip,sum_is,sum_isa
      integer   i,ios,j,k,ktr,l,nishlf
      character preambl(5)*79,cline*80
!
      real       aspmax
      parameter (aspmax=2.0)  ! maximum grid aspect ratio for diffusion
!     parameter (aspmax=1.0)  ! ignore  grid aspect ratio in  diffusion
!
! --- read grid location,spacing,coriolis arrays
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading grid file from ', &
                               trim(flnmgrd),'.[ab]'
        open (unit=uoff+9,file=trim(flnmgrd)//'.b', &
              status='old')
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read(cline,*) i
!
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read (cline,*) j
!
      if     (i.ne.itdm .or. j.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
          'error - wrong array size in grid file'
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      if     (mnproc.eq.1) then
      write (lp,'(a)') trim(cline)
      endif
      read (cline,*) mapflg
!
      call zaiopf(trim(flnmgrd)//'.a','old', 9)
!
      do k= 1,15
        call zagetc(cline,ios, uoff+9)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
          endif !1st tile
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) hminb,hmaxb
        if     (mnproc.eq.1) then
        write (lp,'(a)') trim(cline)
        endif
        call xcsync(flush_lp)
!
        if     (k.eq.1) then
          call zaiord(plon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.2) then
          call zaiord(plat, ip,.false., hmina,hmaxa, 9)
          do i= 1,2  !skip qlon,qlat
            call zagetc(cline,ios, uoff+9)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)') &
                  'geopar: I/O error from zagetc, iunit,ios = ', &
                  uoff+9,ios
              endif !1st tile
              call xcstop('(geopar)')
                     stop '(geopar)'
            endif
#if defined (USE_NUOPC_CESMBETA)
!           qlon, qlat
            j = index(cline,'=')
            read (cline(j+1:),*) hminb,hmaxb
            if     (mnproc.eq.1) then
              write (lp,'(a)') trim(cline)
            endif
            call xcsync(flush_lp)
            call zaiord(util2, ip,.false., hmina,hmaxa, 9)
            if (i.eq.1) then
              qlon(:,:)=util2(:,:)
            else
              qlat(:,:)=util2(:,:)
            endif
#else
!           skip qlon, qlat
            call zaiosk(9)
#endif
          enddo
        elseif (k.eq.3) then
          call zaiord(ulon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.4) then
          call zaiord(ulat, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.5) then
          call zaiord(vlon, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.6) then
          call zaiord(vlat, ip,.false., hmina,hmaxa, 9)
          call zagetc(cline,ios, uoff+9)
          if     (ios.ne.0) then
            if     (mnproc.eq.1) then
              write(lp,'(/ a,i4,i9 /)') &
                'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
            endif !1st tile
            call xcstop('(geopar)')
                   stop '(geopar)'
          endif
#if defined(ESPC_COUPLE) || defined (USE_NUOPC_CESMBETA)
!         pang
          i = index(cline,'=')
          read (cline(i+1:),*) hminb,hmaxb
          if     (mnproc.eq.1) then
          write (lp,'(a)') trim(cline)
          endif
          call xcsync(flush_lp)
          call zaiord(pang, ip,.false., hmina,hmaxa, 9)
#else
!         skip pang
          call zaiosk(9)
#endif
        elseif (k.eq.7) then
          call zaiord(scpx, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.8) then
          call zaiord(scpy, ip,.false., hmina,hmaxa, 9)
        elseif (k.eq.9) then
          call zaiord(scqx, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.10) then
          call zaiord(scqy, iq,.false., hmina,hmaxa, 9)
        elseif (k.eq.11) then
          call zaiord(scux, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.12) then
          call zaiord(scuy, iu,.false., hmina,hmaxa, 9)
        elseif (k.eq.13) then
          call zaiord(scvx, iv,.false., hmina,hmaxa, 9)
        elseif (k.eq.14) then
          call zaiord(scvy, iv,.false., hmina,hmaxa, 9)
        else
          call zaiord(corio,iq,.false., hmina,hmaxa, 9)
        endif
!
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. &
                abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            '.a,.b min = ',hmina,hminb,hmina-hminb, &
            '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          endif
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
      enddo
!
      call zaiocl(9)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=uoff+9)
      endif
!
      if (itest.gt.0 .and. jtest.gt.0) then
        i=itest
        j=jtest
        write (lp,'(/ a,2i5,a,f8.3,a,f12.9,2f10.2/)') &
         ' i,j=',i+i0,j+j0, &
         ' plat=',plat(i,j), &
         ' corio,scux,vy=',corio(i,j),scux(i,j),scvy(i,j)
      endif
      call xcsync(flush_lp)
!
! --- read basin depth array
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        write (lp,'(3a)') ' reading bathymetry file from ', &
                               trim(flnmdep),'.[ab]'
        open (unit=uoff+9,file=trim(flnmdep)//'.b', &
              status='old')
        read (     uoff+9,'(a79)')  preambl
      endif
      call xcsync(flush_lp)
      call zagetc(cline,ios, uoff+9)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
        endif !1st tile
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close(unit=uoff+9)
        write (lp,'(/(1x,a))') preambl,cline
      endif
!
      call zaiopf(trim(flnmdep)//'.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)
      call zaiocl(9)
!
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. &
              abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
          'error - .a and .b files not consistent:', &
          '.a,.b min = ',hmina,hminb,hmina-hminb, &
          '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
        endif
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j= 1,jj
        do i= 1,ii
          if     (depths(i,j).gt.0.5*hugel) then
            depths(i,j) = 0.0
          endif
        enddo
      enddo
!
! --- determine do-loop limits for u,v,p,q points, and update halo for depths
      call bigrid(depths, mapflg, util1,util2,util3)
!cc      call prtmsk(ip,depths,util1,idm,ii,jj,0.0,1.0,
!cc     &     'bottom depth (m)')
!
!     now safe to apply halo to arrays.
!
      vland = 1.0
      call xctilr(plon,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(plat,  1,1, nbdy,nbdy, halo_ps)
#if defined(USE_NUOPC_CESMBETA)
      call xctilr(qlon,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(qlat,  1,1, nbdy,nbdy, halo_ps)
#endif
#if defined(ESPC_COUPLE) || defined (USE_NUOPC_CESMBETA)
      call xctilr(pang,  1,1, nbdy,nbdy, halo_ps)
#endif
      if     (momtyp.eq.2) then
      call xctilr(scpx,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scpy,  1,1, nbdy,nbdy, halo_ps)
      call xctilr(scqx,  1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqy,  1,1, nbdy,nbdy, halo_qs)
      endif !momtyp=2
      call xctilr(corio, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(ulon,  1,1, nbdy,nbdy, halo_us)
      call xctilr(ulat,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scux,  1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy,  1,1, nbdy,nbdy, halo_us)
      call xctilr(vlon,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(vlat,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvx,  1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy,  1,1, nbdy,nbdy, halo_vs)
      vland = 0.0
!
! --- momtum4 needs sc[pq][xy] defined everywhere
!
      if     (momtyp.ne.2) then
        allocate( g_sc(itdm,jtdm) )
        call zaiopf(trim(flnmgrd)//'.a','old', 9)
        do k= 1,9
          call zaiosk(9)
        enddo
        call zaiordg(g_sc, 9)
        call geopar_halo(scpx, g_sc, halo_ps)
        call zaiordg(g_sc, 9)
        call geopar_halo(scpy, g_sc, halo_ps)
        call zaiordg(g_sc, 9)
        call geopar_halo(scqx, g_sc, halo_qs)
        call zaiordg(g_sc, 9)
        call geopar_halo(scqy, g_sc, halo_qs)
        call zaiocl(9)
        deallocate(g_sc)
      endif !momtyp/=2
!
! --- area of grid cells (length x width) at u,v,p,q points resp.
!
!*****!$OMP PARALLEL DO PRIVATE(j,i)
!*****!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          scu2(i,j)=scux(i,j)*scuy(i,j)
          scv2(i,j)=scvx(i,j)*scvy(i,j)
          scp2(i,j)=scpx(i,j)*scpy(i,j)
          scq2(i,j)=scqx(i,j)*scqy(i,j)
!
          scuxi(i,j)=1.0/max(scux(i,j),epsil)
          scvyi(i,j)=1.0/max(scvy(i,j),epsil)
          scp2i(i,j)=1.0/max(scp2(i,j),epsil)
          scq2i(i,j)=1.0/max(scq2(i,j),epsil)
!
! ---     largest grid spacing (within limits) used in all diffusion
! ---     coefficients: min(max(sc?x,sc?y),sc?x*aspmax,sc?y*aspmax)
          aspux(i,j)=min(max(scux(i,j),scuy(i,j)), &
                         min(scux(i,j),scuy(i,j))*aspmax) &
                     /max(scux(i,j),epsil)
          aspuy(i,j)=min(max(scux(i,j),scuy(i,j)), &
                         min(scux(i,j),scuy(i,j))*aspmax) &
                     /max(scuy(i,j),epsil)
          aspvx(i,j)=min(max(scvx(i,j),scvy(i,j)), &
                         min(scvx(i,j),scvy(i,j))*aspmax) &
                     /max(scvx(i,j),epsil)
          aspvy(i,j)=min(max(scvx(i,j),scvy(i,j)), &
                         min(scvx(i,j),scvy(i,j))*aspmax) &
                     /max(scvy(i,j),epsil)
!
          util1(i,j)=depths(i,j)*scp2(i,j)
        enddo
      enddo
!
! --- read ice shelf depth array
!
      if     (ishelf.eq.0) then
        ishlf(:,:) = ip(:,:)  !no ice shelf
      else
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          write (lp,'(3a)') ' reading ice shelf file from ', &
                                 trim(flnmshlf),'.[ab]'
          open (unit=uoff+9,file=trim(flnmshlf)//'.b', &
                status='old')
          read (     uoff+9,'(a79)')  preambl
        endif
        call xcsync(flush_lp)
        call zagetc(cline,ios, uoff+9)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
          endif !1st tile
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*)   hminb,hmaxb
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          close(unit=uoff+9)
          write (lp,'(/(1x,a))') preambl,cline
        endif
!
        call zaiopf(trim(flnmshlf)//'.a','old', 9)
        call zaiord(util3,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
!
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. &
                abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            '.a,.b min = ',hmina,hminb,hmina-hminb, &
            '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          endif
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.0) then
              util3(i,j) = 0.0  !land
            elseif (util3(i,j).gt.0.5*hugel) then
              util3(i,j) = 0.0  !ice shelf over ocean
            elseif (util3(i,j).le.0.0) then
              util3(i,j) = 0.0  !ice shelf over ocean
            else
              util3(i,j) = 1.0  !open ocean
            endif
          enddo
        enddo
        call xctilr(util3,1,1, nbdy,nbdy, halo_ps)
        ishlf(:,:) = 0  !for jj:jdm and ii:idm
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            ishlf(i,j) = util3(i,j)
            util2(i,j) = ip(i,j)
          enddo
        enddo
!
        call xcsum(sum_is,  util3,ip)
        call xcsum(sum_ip,  util2,ip)
        call xcsum(sum_isa, scp2, ishlf)
        call xcsum(area,    scp2, ip)
        nishlf = nint(sum_ip) - nint(sum_is)
        if     (mnproc.eq.1) then
        write (lp,'(/a,i9,f10.2)') &
               ' number of ice shelf points and area (10^6 km^2):', &
               nishlf,(area-sum_isa)*1.d-12
        endif
        call xcsync(flush_lp)
      endif !ishelf
!
! --- In arctic (tripole) domain, top row of mass points is redundent,
! ---  so always use ipa, based on ishlf, for mass sums
#if defined(ARCTIC)
      ipa(:,:) = ishlf(:,:)
      if     (jj+j0.eq.jtdm) then
! ---   mask top row of mass points
        ipa(:,jj:jj+nbdy) = 0
      endif
#else
! --- Not a tripole domain, so ipa=ishlf
      ipa(:,:) = ishlf(:,:)
#endif
!
      call xcsum(avgbot, util1,ipa)
      call xcsum(area,   scp2, ipa)
      avgbot=avgbot/area
      if     (mnproc.eq.1) then
      write (lp,'(/a,f9.1,f10.2)') &
             ' mean basin depth (m) and area (10^6 km^2):', &
             avgbot,area*1.e-12
      endif
      call xcsync(flush_lp)
!
! --- calculate dp0k and ds0k?
      if     (dp00.lt.0.0) then
! ---   dp0k and ds0k already input
        dp00 =onem*dp0k(1)
        dp00x=onem*dp0k(kk-1)
        dp00i=onem*dp00i
        dpms = 0.0
        do k=1,kk
          dpm     = dp0k(k)
          dpms    = dpms + dpm
          dp0k(k) = dp0k(k)*onem
          if     (mnproc.eq.1) then
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
          call xcsync(flush_lp)
        enddo !k
        dsms = 0.0
        do k=1,nsigma
          dsm     = ds0k(k)
          dsms    = dsms + dsm
          ds0k(k) = ds0k(k)*onem
          if     (mnproc.eq.1) then
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
          call xcsync(flush_lp)
        enddo !k
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
      else
! ---   calculate dp0k and ds0k
!
! ---   logorithmic k-dependence of dp0 (deep z's)
        dp00 =onem*dp00
        dp00x=onem*dp00x
        dp00i=onem*dp00i
        if     (isopyc) then
          dp0k(1)=thkmin*onem
        else
          dp0k(1)=dp00
        endif
        dpm  = dp0k(1)*qonem
        dpms = dpm
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,135) 1,dp0k(1)*qonem,dpm,dpms
        endif
 135    format('dp0k(',i2,') =',f7.2,' m', &
                  '    thkns =',f7.2,' m', &
                  '    depth =',f8.2,' m')
        call xcsync(flush_lp)
!
        dp0kf=1.0
        do k=2,kk
          dp0kf=dp0kf*dp00f
          if     (k.le.nhybrd) then
            if     (dp00f.ge.1.0) then
              dp0k(k)=min(dp00*dp0kf,dp00x)
            else
              dp0k(k)=max(dp00*dp0kf,dp00x)
            endif
          else
            dp0k(k)=0.0
          endif
          dpm  = dp0k(k)*qonem
          dpms = dpms + dpm
          if     (mnproc.eq.1) then
          write(lp,135) k,dp0k(k)*qonem,dpm,dpms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: dp0kf  = ',dp0kf,    mnproc
            write(6,*) 'geopar: dp0k   = ',dp0k(k),k,mnproc
          endif
          call xcsync(flush_lp)
        enddo !k
!
! ---   logorithmic k-dependence of ds0 (shallow z-s)
        ds00 =onem*ds00
        ds00x=onem*ds00x
        if     (isopyc) then
          ds0k(1)=thkmin*onem
        else
          ds0k(1)=ds00
        endif
        dsm  = ds0k(1)*qonem
        dsms = dsm
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,130) 1,ds0k(1)*qonem,dsm,dsms
        endif
 130    format('ds0k(',i2,') =',f7.2,' m', &
                  '    thkns =',f7.2,' m', &
                  '    depth =',f8.2,' m')
        call xcsync(flush_lp)
!
        ds0kf=1.0
        do k=2,nsigma
          ds0kf=ds0kf*ds00f
          if     (ds00f.ge.1.0) then
            ds0k(k)=min(ds00*ds0kf,ds00x)
          else
            ds0k(k)=max(ds00*ds0kf,ds00x)
          endif
          dsm  = ds0k(k)*qonem
          dsms = dsms + dsm
          if     (mnproc.eq.1) then
          write(lp,130) k,ds0k(k)*qonem,dsm,dsms
          endif
          if     (mnproc.eq.-99) then  ! bugfix that prevents optimization
            write(6,*) 'geopar: ds0kf  = ',ds0kf,    mnproc
            write(6,*) 'geopar: ds0k   = ',ds0k(k),k,mnproc
          endif
          call xcsync(flush_lp)
        enddo !k
        if     (mnproc.eq.1) then
        write(lp,*)
        endif
      endif !input:calculate dp0k,ds0k
!
! --- start and stop depths for terrain following coordinate
      if     (nsigma.eq.0) then
        dpns    = dp0k(1)
        dsns    = 0.0
        ds0k(1) = dp0k(1)
        do k= 2,kk
          ds0k(k)=0.0
        enddo !k
      else
        dpns = 0.0
        dsns = 0.0
        do k=1,nsigma
          dpns = dpns + dp0k(k)
          dsns = dsns + ds0k(k)
        enddo !k
        do k= nsigma+1,kk
          ds0k(k)=0.0
        enddo !k
      endif !nsigma
      dpns = dpns*qonem  !depths is in m
      dsns = dsns*qonem  !depths is in m
!
      if     (mnproc.eq.1) then
      write(lp,131) nsigma,dpns,dsns
      endif
 131  format('nsigma = ',i2, &
             '    deep    =',f8.2,' m', &
             '    shallow =',f8.2,' m' )
      call flush(lp)
!
! --- initialize thermobaric reference state arrays.
!
      if     (kapref.eq.-1) then
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          write (lp,'(3a)') ' reading thermobaric reference file from ', &
                                 trim(flnmforw), 'tbaric.[ab]'
          open (unit=uoff+9,file=trim(flnmforw)//'tbaric.b', &
                status='old')
          read (     uoff+9,'(a79)')  preambl
        endif
        call xcsync(flush_lp)
        call zagetc(cline,ios, uoff+9)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'geopar: I/O error from zagetc, iunit,ios = ',uoff+9,ios
          endif !1st tile
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*)   hminb,hmaxb
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          close(unit=uoff+9)
          write (lp,'(/(1x,a))') preambl,cline
        endif
!
! ---   input field is between 1.0 and 3.0 and indicates the
! ---   relative strength of the two nearest reference states,
! ---     e.g. 1.7 is 70% ref2 and 30% ref1
! ---     and  2.3 is 70% ref2 and 30% ref3.
!
        call zaiopf(trim(flnmforw)//'tbaric.a','old', 9)
        call zaiord(util1,ip,.false., hmina,hmaxa, 9)
        call zaiocl(9)
!
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or. &
                abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            '.a,.b min = ',hmina,hminb,hmina-hminb, &
            '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
          endif
          call xcstop('(geopar)')
                 stop '(geopar)'
        endif
!
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.0) then
              util1(i,j) = 1.0 !land
            endif
          enddo
        enddo
!
        vland = 1.0
        call xctilr(util1,  1,1, nbdy,nbdy, halo_ps)
        vland = 0.0
!
!       kapi is the 2nd reference state (1st is always 2)
!       skap is the scale factor (0.0-1.0) for the 1st reference state
!
!       assumes that reference states 1 and 3 are never next to each other.
!
        do j= 1,jj
          do i= 1,ii
            if     (max(util1(i,  j), &
                        util1(i-1,j), &
                        util1(i+1,j), &
                        util1(i,  j-1), &
                        util1(i,  j+1) ).gt.2.0) then
              util2(i,j) = 3.0              !kapi
               skap(i,j) = 3.0 - util1(i,j)
            else
              util2(i,j) = 1.0              !kapi
               skap(i,j) = util1(i,j) - 1.0
            endif
          enddo
        enddo
        vland = 1.0
        call xctilr(util2, 1,1, nbdy,nbdy, halo_ps)
        call xctilr(skap,  1,1, nbdy,nbdy, halo_ps)
        vland = 0.0
!
        kapi(:,:) = util2(:,:)
      else
        skap(:,:) = 1.0     !for diagnostics only
        kapi(:,:) = kapref  !for diagnostics only
      endif !kapref.eq.-1:else
!
! --- initialize some arrays
! --- set depthu,dpu,utotn,pgfx,depthv,dpv,vtotn,pgfy to zero everywhere,
! --- so that they can be used at "lateral neighbors" of u and v points.
! --- similarly for pbot,dp at neighbors of q points.
!
      disp_count=0
!
!$OMP PARALLEL DO PRIVATE(j,i,k,ktr) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          p(     i,j,1)=0.0
          pu(    i,j,1)=0.0
          pv(    i,j,1)=0.0
          utotn( i,j)=0.0
          vtotn( i,j)=0.0
          pgfx(  i,j)=0.0
          pgfy(  i,j)=0.0
          gradx( i,j)=0.0
          grady( i,j)=0.0
          depthu(i,j)=0.0
          depthv(i,j)=0.0
          pbot(  i,j)=0.0
!
          displd_mn(i,j)=0.0
          dispqd_mn(i,j)=0.0
          tidepg_mn(i,j)=0.0
!
          psikk( i,j,1)=0.0
          psikk( i,j,2)=0.0
          thkk(  i,j,1)=0.0
          thkk(  i,j,2)=0.0
!
          ubavg( i,j,1)=hugel
          ubavg( i,j,2)=hugel
          ubavg( i,j,3)=hugel
          vbavg( i,j,1)=hugel
          vbavg( i,j,2)=hugel
          vbavg( i,j,3)=hugel
!
           umix( i,j)=hugel
           vmix( i,j)=hugel
          utotm( i,j)=hugel
          vtotm( i,j)=hugel
          uflux( i,j)=hugel
          vflux( i,j)=hugel
          uflux2(i,j)=hugel
          vflux2(i,j)=hugel
          uflux3(i,j)=hugel
          vflux3(i,j)=hugel
          uja(   i,j)=hugel
          ujb(   i,j)=hugel
          via(   i,j)=hugel
          vib(   i,j)=hugel
          do k=1,kk
            dp( i,j,k,1)=0.0
            dp( i,j,k,2)=0.0
            dpo(i,j,k,1)=0.0
            dpo(i,j,k,2)=0.0
            dpu(i,j,k,1)=0.0
            dpu(i,j,k,2)=0.0
            dpv(i,j,k,1)=0.0
            dpv(i,j,k,2)=0.0
!
            u(  i,j,k,1)=hugel
            u(  i,j,k,2)=hugel
            v(  i,j,k,1)=hugel
            v(  i,j,k,2)=hugel
!
            uflx(  i,j,k)=hugel
            vflx(  i,j,k)=hugel
!
            dpav(  i,j,k)=0.0
            uflxav(i,j,k)=0.0
            vflxav(i,j,k)=0.0
            diaflx(i,j,k)=0.0
!
            do ktr= 1,ntracr
              tracer(i,j,k,1,ktr)=0.0
              tracer(i,j,k,2,ktr)=0.0
            enddo
          enddo
        enddo
      enddo
!$OMP END PARALLEL DO
!
!$OMP PARALLEL DO PRIVATE(j,l,i,k) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l=1,isp(j) !ok
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
            ubavg(i,j,1)=0.0
            ubavg(i,j,2)=0.0
            ubavg(i,j,3)=0.0
            umix  (i,j)=0.0
            utotm (i,j)=0.0
            uflux (i,j)=0.0
            uflux2(i,j)=0.0
            uflux3(i,j)=0.0
            uja(i,j)=0.0
            ujb(i,j)=0.0
!
            do k=1,kk
              uflx(i,j,k)=0.0
              u(i,j,k,1)=0.0
              u(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
!
      call xctilr(ubavg,    1,   3, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(umix,     1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(utotm,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux,    1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux2,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uflux3,   1,   1, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(uja,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(ujb,      1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflx,     1,  kk, nbdy,nbdy, halo_us)  ! note scalar
      call xctilr(u,        1,2*kk, nbdy,nbdy, halo_us)  ! note scalar
!
!$OMP PARALLEL DO PRIVATE(i,l,j,k) &
!$OMP          SCHEDULE(STATIC)
      do i=1,ii
        do l=1,jsp(i) !ok
          do j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
            vbavg(i,j,1)=0.0
            vbavg(i,j,2)=0.0
            vbavg(i,j,3)=0.0
            vmix  (i,j)=0.0
            vtotm (i,j)=0.0
            vflux (i,j)=0.0
            vflux2(i,j)=0.0
            vflux3(i,j)=0.0
            via(i,j)=0.0
            vib(i,j)=0.0
!
            do k=1,kk
              vflx(i,j,k)=0.0
              v(i,j,k,1)=0.0
              v(i,j,k,2)=0.0
            enddo
          enddo
        enddo
      enddo
!
      call xctilr(vbavg,    1,   3, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vmix,     1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vtotm,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux,    1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux2,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflux3,   1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(via,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vib,      1,   1, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(vflx,     1,  kk, nbdy,nbdy, halo_vs)  ! note scalar
      call xctilr(v,        1,2*kk, nbdy,nbdy, halo_vs)  ! note scalar
!
      return
      end

      subroutine geopar_halo(sc, g_sc, halo_type)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer     halo_type
      real        sc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      real*4      g_sc(itdm,jtdm)
!
! --- define the halo of sc, from g_sc
! --- needed because xctilr may use vland in the halo
!
      integer i,ia,ig,j,ja,jg
!
#if defined(ARCTIC)
! --- periodic in x, closed in South, dipole patch in North
!
      do j=1,jj
        jg = j0+j
        do i=1-nbdy,0
          ig = mod(i0+i-1+itdm,itdm)+1
          sc(i,j) = g_sc(ig,jg)
        enddo !i
        do i=ii+1,ii+nbdy
          ig = mod(i0+i-1+itdm,itdm)+1
          sc(i,j) = g_sc(ig,jg)
        enddo !i
      enddo !j
!
      do i=1-nbdy,ii+nbdy
        ig = mod(i0+i-1+itdm,itdm)+1
        do j=1-nbdy,0
          jg = max(1,j0+j)
          sc(i,j) = g_sc(ig,jg)
        enddo !j
        do j=jj+1,jj+nbdy
          jg = j0+j
          if     (jg.le.jtdm) then
            sc(i,j) = g_sc(ig,jg)
          else
            if     (halo_type.eq.halo_ps) then
              ja = jtdm-1-(jg-jtdm)
              ia = itdm-mod(ig-1,itdm)
              sc(i,j) = g_sc(ia,ja)
            else   !halo_type.eq.halo_qs
              ja = jtdm-(jg-jtdm)
              ia = mod(itdm-(ig-1),itdm)+1
              sc(i,j) = g_sc(ia,ja)
            endif
          endif
        enddo !j
      enddo !i
#else
      do j=1,jj
        jg = j0+j
        do i=1-nbdy,0
          if     (nreg.le.2) then
            ig = max(   1,i0+i)           !closed
          else
            ig = mod(i0+i-1+itdm,itdm)+1  !periodic
          endif
          sc(i,j) = g_sc(ig,jg)
        enddo !i
        do i=ii+1,ii+nbdy
          if     (nreg.le.2) then
            ig = min(itdm,i0+i)           !closed
          else
            ig = mod(i0+i-1+itdm,itdm)+1  !periodic
          endif
          sc(i,j) = g_sc(ig,jg)
        enddo !i
      enddo !j
!
      do i=1-nbdy,ii+nbdy
        if     (nreg.le.2) then
          ig = max(1, min( itdm, i0+i ) ) !closed
        else
          ig = mod(i0+i-1+itdm,itdm)+1    !periodic
        endif
        do j=1-nbdy,0
          if     (nreg.eq.0 .or. nreg.eq.4) then
            jg = max(   1,j0+j)           !closed
          else
            jg = mod(j0+j-1+jtdm,jtdm)+1  !periodic
          endif
          sc(i,j) = g_sc(ig,jg)
        enddo !j
        do j=jj+1,jj+nbdy
          if     (nreg.eq.0 .or. nreg.eq.4) then
            jg = min(jtdm,j0+j)           !closed
          else
            jg = mod(j0+j-1+jtdm,jtdm)+1  !periodic
          endif
          sc(i,j) = g_sc(ig,jg)
        enddo !j
      enddo !i
#endif
!diag 
!diag if (itest.gt.0 .and. jtest.gt.0) then
!diag   do j= 1,jtdm
!diag     write(lp,'(a,2i5,e16.8)') &
!diag     'g_sc =',ittest,j,g_sc(ittest,j)
!diag   enddo
!diag   do i= 1,itdm
!diag     write(lp,'(a,2i5,e16.8)') &
!diag     'g_sc =',i,jttest,g_sc(i,jttest)
!diag   enddo
!diag   i = itest
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'sc   (km) =', &
!diag             i+i0,1+j0, &
!diag             (sc(i,j)*1.e-3,j=1,-3,-1)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'g_sc (km) =', &
!diag             i+i0,1+j0, &
!diag             (g_sc(i+i0,j+j0)*1.e-3,j=1,-3,-1)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'sc   (km) =', &
!diag             i+i0,jj+j0, &
!diag             (sc(i,j)*1.e-3,j=jj,jj+4)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'g_sc (km) =', &
!diag             i+i0,jj+j0, &
!diag             (g_sc(i+i0,j+j0)*1.e-3,j=jj,jj+4)
!diag   j = jtest
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'sc   (km) =', &
!diag             1+i0,j+j0, &
!diag             (sc(i,j)*1.e-3,i=1,-3,-1)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'g_sc (km) =', &
!diag             1+i0,j+j0, &
!diag             (g_sc(i+i0,j+j0)*1.e-3,i=1,-3,-1)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'sc   (km) =', &
!diag             ii+i0,j+j0, &
!diag             (sc(i,j)*1.e-3,i=ii,ii+4)
!diag   write(lp,'(a,2i5,5f10.5)') &
!diag             'g_sc (km) =', &
!diag             ii+i0,j+j0, &
!diag             (g_sc(i+i0,j+j0)*1.e-3,i=ii,ii+4)
!diag endif !test
      return
      end
!
!
!> Revision history:
!>
!> May  1997 - extended list of variables set to 'hugel' on land
!> Oct. 1999 - added code that defines the vertical distribution of dp0
!>             used in hybgen
!> Jan. 2000 - added mapflg logic for different projections
!> Feb. 2000 - added dp00f for logorithmic z-level spacing
!> Mar. 2000 - added dp00s for sigma-spacing in shallow water
!> May  2000 - conversion to SI units (still wrong corio)
!> Feb. 2001 - removed rotated grid option
!> Jan. 2002 - more flexible Z-sigma-Z vertical configuration
!> Jan. 2002 - all grids now via array input
!> Sep. 2004 - define kapi and skap for thermobaricity
!> Oct. 2008 - dp0k and ds0k can now be input, see blkdat.F
!> Mar. 2012 - replaced dssk with dpns and dsns
!> Apr. 2014 - added ishlf
!> Apr. 2014 - added ipa
!> Feb. 2015 - added pang for coupled cases
!> July 2017 - for momtum4, calculate accurate halos for sc[pq][xy]
!> Aug. 2018 - initialize umix and vmix
!> Dec. 2018 - add /* USE_NUOPC_CESMBETA */ macro for pang (for coupled simulation)
