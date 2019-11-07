#if defined (USE_NUOPC_CESMBETA)
      subroutine blkdat(hycom_start_dtg)
#else
      subroutine blkdat
#endif
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_incupd     ! HYCOM incremental update (for data assimilation)
      use mod_floats     ! HYCOM synthetic floats, drifters and moorings
      use mod_tides      ! HYCOM tides
#if defined(STOKES)
      use mod_stokes     ! HYCOM Stokes Drift effects from Wave fields
#endif
      implicit none
!
      real      day1,hybrlx,cplifq,frzifq
      integer   k,kdmblk,mlflag,thflag,trcflg1
      integer   lngblk
      character sigfmt*26
#if defined (USE_NUOPC_CESMBETA)
      real,         intent(in):: hycom_start_dtg
#endif
!
# include "stmt_fns.h"
!
! --- initialize common variables.
!
#if defined(OCEANS2)
      if     (nocean.eq.2) then
! ---   slave HYCOM works from ./OCEAN2
        flnminp = './OCEAN2/'
      else
! ---   master HYCOM
        flnminp = './'
      endif
#else
      flnminp = './'
#endif
      open(unit=uoff+99,file=trim(flnminp)//'blkdat.input')  !on all nodes
!
! --- 'lp' = logical unit number for printer output
!     lp = 6  !now defined in xcspmd
!
! --- four lines (80-characters) describing the simulation
      read( uoff+99,'(a80/a80/a80/a80)') ctitle
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,'(a80/a80/a80/a80)') ctitle
      call flush(lp)
      endif !1st tile
!
! --- 'iversn' = hycom version number x10
! --- 'iexpt'  = experiment number x10
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iversn,'iversn')
      call blkini(iexpt, 'iexpt ')
!
      if (iversn.lt.23 .or. iversn.gt.23) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3,a,i3 /)')  &
          'error - iversn must be between',23,' and',23
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
!
! --- 'idm   ' = longitudinal array size
! --- 'jdm   ' = latitudinal  array size
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(itest ,'idm   ')
      call blkini(jtest ,'jdm   ')
!
      if     (itdm.eq.-1) then !serial relo case only
        itdm = itest
        idm  = itest
        ii   = itest
        jtdm = jtest
        jdm  = jtest
        jj   = jtest
      endif
!
      if     (itest.ne.itdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)')  &
          'error - expected idm =',itdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (jtest.ne.jtdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)')  &
          'error - expected jdm =',jtdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
!
! --- 'itest,jtest' = grid point where detailed diagnostics are desired
! ---                 itest=jtest=0 turns off all detailed diagnostics
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(ittest,'itest ')
      call blkini(jttest,'jtest ')
!
      if (ittest.gt.itdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)')  &
          'error - maximum itest is',itdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if (jttest.gt.jtdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i5 /)')  &
          'error - maximum jtest is',jtdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
!
! --- map global ittest,jttest to local itest,jtest
      if     (ittest.gt.i0 .and. ittest.le.i0+ii .and. &
              jttest.gt.j0 .and. jttest.le.j0+jj      ) then
        itest = ittest - i0
        jtest = jttest - j0
      else
        itest = -99
        jtest = -99
      endif
!
!     if (mnproc.eq.1) then
!     write(lp,*)
!     endif !1st tile
!     call xcsync(flush_lp)
!     do k= 1,ijpr
!       if     (mnproc.eq.k) then
!         write(lp,'(a,3i5)') 'mnproc,[ij]test =',mnproc,itest,jtest
!       endif
!       call xcsync(flush_lp)
!     enddo !k
!
! --- 'kdm   ' = number of layers
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(kdmblk,'kdm   ')
!
#if defined(RELO)
      kdm = kdmblk
      kk  = kdmblk
!
      allocate( &
          dp0k(kdm), &
          ds0k(kdm), &
         sigma(kdm), &
        salmin(kdm) )
      call mem_stat_add( 4*kdm )
#else
      if     (kdmblk.ne.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)')  &
          'error - expected kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
#endif
!
! --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
! --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
! --- 'dp00'   = deep    z-level spacing minimum thickness (m)
! --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
! --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
! --- 'ds00'   = shallow z-level spacing minimum thickness (m)
! --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
! --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
! --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
! --- 'isotop' = shallowest depth for isopycnal layers     (m), <0 from file
!
! --- the above specifies a vertical coord. that is isopycnal or:
! ---  near surface z in    deep water, based on dp00,dp00x,dp00f
! ---  near surface z in shallow water, based on ds00,ds00x,ds00f and nsigma
! ---               sigma between them, based on ds00,ds00x,ds00f and nsigma
!
! --- d[ps].k/d[ps].1 = d[ps]00f**(k-1) unless limited by d[ps]00x.
! --- if d[ps]00f>1, d[ps]00x (> d[ps]00) is the maximum thickness.
! --- if d[ps]00f<1, d[ps]00x (< d[ps]00) is the minimum thickness.
!
! --- near the surface (i.e. shallower than isotop), layers are always fixed
! --- depth (z or sigma).  layer 1 is always fixed, so isotop=0.0 is not
! --- allowed.  if isotop<0.0 then isotop(1:idm,1:jdm) is a spacially
! --- varying array, input from the file iso.top.[ab].
!
! --- away from the surface, the minimum layer thickness is dp00i.
! --- to recover original hybrid behaviour,   set dp00i=dp00x
! --- for z-only, sigma-only or sigma-z only, set dp00i=dp00x
!
! --- for z-only set nsigma=0 (and ds00,ds00x,ds00f=dp00,dp00x,dp00f)
! --- for sigma-z (shallow-deep) use a very small ds00
! ---  (pure sigma-z also has ds00f=dp00f and ds00x=dp00x*ds00/dp00)
! --- for z-sigma (shallow-deep) use a very large dp00 (not recommended)
! --- for sigma-only set nsigma=kdm, dp00 large, and ds00 small
!
! --- for an entirely fixed vertical coordinate (no isopycnal layers), set
! --- isotop large or make all target densities (sigma(k), below) very small.
!
! --- or, in place of 'dp00','dp00x','dp00f','ds00','ds00x','ds00f' specify:
! --- 'dp0k  ' = layer k deep    z-level spacing minimum thickness (m)
! ---              k=1,kdm; dp0k must be zero for k>nhybrd
! --- 'ds0k  ' = layer k shallow z-level spacing minimum thickness (m)
! ---              k=1,nsigma
!
! --- uses version 2.2.58 definition of deep and shallow z-levels.
! --- terrain following starts at depth dpns=sum(dp0k(k),k=1,nsigma) and
! --- ends at depth dsns=sum(ds0k(k),k=1,nsigma), and the depth of the
! --- k-th layer interface varies linearly with total depth between
! --- these two reference depths.
!
! --- previous to 2.2.58, it was layer thickness (not layer interface
! --- depth) that varied linearly with total depth.  These two approachs
! --- are identical for "pure sigma-z", but differ if ds00f/=dp00f.
!
      call blkini(nhybrd,'nhybrd')
      call blkini(nsigma,'nsigma')
      call blkinr2(dp00,k, &
                         'dp00  ','(a6," =",f10.4," m")', &
                         'dp0k  ','(a6," =",f10.4," m")' )
      if     (k.eq.1) then !dp00
        call blkinr(dp00x, 'dp00x ','(a6," =",f10.4," m")')
        call blkinr(dp00f, 'dp00f ','(a6," =",f10.4," ")')
        call blkinr(ds00,  'ds00  ','(a6," =",f10.4," m")')
        call blkinr(ds00x, 'ds00x ','(a6," =",f10.4," m")')
        call blkinr(ds00f, 'ds00f ','(a6," =",f10.4," ")')
      else !dp0k
        dp0k(1) = dp00
        dp00    = -1.0  !signal that dp00 is not input
        do k=2,kdm
          call blkinr(dp0k(k), 'dp0k  ','(a6," =",f10.4," m")')
!
          if      (k.gt.nhybrd .and. dp0k(k).ne.0.0) then
            if (mnproc.eq.1) then
            write(lp,'(/ a,i3 /)')  &
              'error - dp0k must be zero for k>nhybrd'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif !k>nhybrd&dp0k(k)!=0
        enddo !k
        do k=1,nsigma
          call blkinr(ds0k(k), 'ds0k  ','(a6," =",f10.4," m")')
        enddo !k
      endif !dp00:dp0k
      call blkinr(dp00i, 'dp00i ','(a6," =",f10.4," m")')
      call blkinr(isotop,'isotop','(a6," =",f10.4," m")')
!
! --- isopycnal (MICOM-like) iff nhybrd is 0
      isopyc = nhybrd .eq. 0
      hybrid = .not. isopyc
      if (hybrid .and. nsigma.le.1) then
        nsigma=1
        if     (dp00.lt.0.0) then
          ds0k(1) = dp0k(1)
        endif
      endif
      if (isopyc .and. sigver.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
          'error - MICOM-like requires the 7-term eqn. of state'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if     (nhybrd.gt.kdm) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)')  &
          'error - maximum nhybrd is kdm =',kdm
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (nsigma.gt.nhybrd) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3 /)')  &
          'error - maximum nsigma is nhybrd =',nhybrd
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !error
      if     (dp00.ge.0.0) then
        if (isopyc .and. max(dp00,dp00x).ne.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - must have dp00x==dp00==0.0 for isopycnal case'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (dp00f.eq.1.0 .and. dp00.ne.dp00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - must have dp00x==dp00 for dp00f==1.0'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (dp00f.gt.1.0 .and. dp00.ge.dp00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - dp00x must be > dp00 for dp00f>1'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (dp00f.lt.1.0 .and. dp00.le.dp00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - dp00x must be < dp00 for dp00f<1'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (ds00.gt.dp00) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - must have ds00 <= dp00'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (ds00.le.0.0 .and. .not.isopyc) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - must have ds00>0.0'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (ds00f.eq.1.0 .and. ds00.ne.ds00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - must have ds00x==ds00 for ds00f==1.0'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (ds00f.gt.1.0 .and. ds00.ge.ds00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - ds00x must be > ds00 for ds00f>1'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (ds00f.lt.1.0 .and. ds00.le.ds00x) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - ds00x must be < ds00 for ds00f<1'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
        if (isotop.eq.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - isotop cannot be 0.0 (layer 1 never isopycnal)'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !error
      endif !dp00 used
!
! --- 'oneta0' = minimum 1+eta, must be > 0.0
! --- 'saln0'  = initial salinity value (psu), only used for iniflg<2
! --- 'locsig' = locally-referenced potential density for stability (0=F,1=T)
! --- 'kapref' = thermobaric reference state (-1=input,0=none,1,2,3=constant)
! ---              1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
! --- 'thflag' = reference pressure flag (0=Sigma-0, 2=Sigma-2)
! ---            this is a check on the compile-time stmt_funcs.h setup.
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(oneta0,'oneta0','(a6," =",f10.4," ")')
      call blkinr(saln0, 'saln0 ','(a6," =",f10.4," psu")')
      call blkinl(locsig,'locsig')
      call blkini(kapref,'kapref')
      call blkini(thflag,'thflag')
! --- kapnum is number of thermobaric reference states (1 or 2)
      if     (kapref.ne.-1) then
        kapnum=1
      else
        kapnum=2
      endif
!
      if     (oneta0.le.0.0 .or. oneta0.ge.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
          'error - oneta0 must be above 0.0 and below 1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !oneta0 error
      if     (kapref.lt.-1 .or. kapref.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i1 /)') &
          'error - kapref must be between -1 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !kapref error
      if     (thflag.eq.0) then
        if     (sigver.eq.1) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 7-term sigma-0'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.3) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 9-term sigma-0'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.5) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 17-term sigma-0'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.7) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 12-term sigma-0'
          call flush(lp)
          endif !1st tile
        else
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !7-term:9-term:error
      elseif (thflag.eq.2) then
        if     (sigver.eq.2) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 7-term sigma-2'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.4) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 9-term sigma-2'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.6) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 17-term sigma-2'
          call flush(lp)
          endif !1st tile
        elseif (sigver.eq.8) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'equation of state is 12-term sigma-2'
          call flush(lp)
          endif !1st tile
        else
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
            'error - thflag not consistent with sig()'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif !7-term:9-term:error
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
          'error - thflag must be 0 or 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !thflag
      if     (thflag.eq.0) then
        sigfmt = '(a6," =",f10.4," sigma-0")'
      elseif (thflag.eq.2) then
        sigfmt = '(a6," =",f10.4," sigma-2")'
      endif
!
! --- 'thbase' = reference density (sigma units)
      call blkinr(thbase,'thbase',sigfmt)
!
! --- 'vsigma' = spacially varying isopycnal layer target densities (0=F,1=T)
! ---            if true, target densities input from file iso.sigma.[ab]
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinl(vsigma,'vsigma')
!
! --- 'sigma ' = isopycnal layer target densities (sigma units)
      do k=1,kdm
        call blkinr(sigma(k),'sigma ',sigfmt)
!
        if     (k.gt.1) then
          if      (sigma(k).le.sigma(k-1)) then
            if (mnproc.eq.1) then
            write(lp,'(/ a,i3 /)')  &
              'error - sigma is not stabally stratified'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif !sigma(k).le.sigma(k-1)
        endif !k>1
      enddo !k
!
! --- 'iniflg' = initial state flag (0=level,1=zonal,2=climatology)
! --- 'jerlv0' = initial jerlov water type (1 to 5; 0 kpar, -1 chl)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(iniflg,'iniflg')
      call blkini(jerlv0,'jerlv0')
!
      if (iniflg.lt.0 .or. iniflg.gt.3) then  !inicon==3 ok, for old .inputs
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
          'error - iniflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (jerlv0.lt.-1 .or. jerlv0.gt.5) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
          'error - jerlv0 must be -1 or 0 or between 1 and 5'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
! --- 'yrflag' = days in year flag   (0=360,1=366,2=366Jan1,3=actual)
! ---             (-2=366Jan1 with 732-day forcing repeat)
! --- 'sshflg' = diagnostic SSH flag (0=SSH,1=SSH&stericSSH,2=SSH&stericMONTG)
! ---             note that sshflg==1 implies reading relax.ssh.a
! ---              and that sshflg==2 implies reading relax.montg.a
! --- 'dsurfq' = number of days between model diagnostics at the surface
! ---             (-ve to output only the .txt file at itest,jtest)
! --- 'diagfq' = number of days between model diagnostics
! ---             (-ve same as -diagfq but always write archive at end)
! --- 'proffq' = number of days between model diagnostics at selected locations
! --- 'tilefq' = number of days between model diagnostics on selected tiles
! --- 'meanfq' = number of days between model diagnostics (time averaged)
! --- 'rstrfq' = number of days between model restart output
! ---             (-ve same as -rstrfq but no restart at end)
! --- 'bnstfq' = number of days between baro nesting archive input
! ---             (-ve for mean archives spanning -bnstfq days;
! ---              0.0 if lbflag is not 2 or 4)
! --- 'nestfq' = number of days between 3-d  nesting archive input
! ---             (-ve for mean archives spanning -nestfq days;
! ---              0.0 turns off relaxation to 3-d nesting input)
! --- 'cplifq' = number of days (or time steps) between sea ice coupling
! ---             (positive days or negative time steps)
! --- 'baclin' = baroclinic time step (seconds), int. divisor of 86400
! --- 'batrop' = barotropic time step (seconds), int. divisor of baclin/2
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(yrflag,'yrflag')
      call blkini(sshflg,'sshflg')
      call blkinr(dsurfq,'dsurfq','(a6," =",f10.4," days")')
      call blkinr(diagfq,'diagfq','(a6," =",f10.4," days")')
      call blkinr(proffq,'proffq','(a6," =",f10.4," days")')
      call blkinr(tilefq,'tilefq','(a6," =",f10.4," days")')
      call blkinr(meanfq,'meanfq','(a6," =",f10.4," days")')
      call blkinr(rstrfq,'rstrfq','(a6," =",f10.4," days")')
      call blkinr(bnstfq,'bnstfq','(a6," =",f10.4, &
                                &" days (-ve mean span)")')
      call blkinr(nestfq,'nestfq','(a6," =",f10.4, &
                                &" days (-ve mean span)")')
      call blkinr(cplifq,'cplifq','(a6," =",f10.4, &
                                &" days (-ve time steps)")')
      call blkinr(baclin,'baclin','(a6," =",f10.4," sec")')
      call blkinr(batrop,'batrop','(a6," =",f10.4," sec")')
!
      if     (yrflag.eq.-2) then
        yrflag = 2
        wndrep = 732.0d0  !two year forcing repeat period
      elseif (yrflag.eq. 2) then
        wndrep = 366.0d0  !one year forcing repeat period
      elseif (yrflag.eq. 4) then
        wndrep = 365.0d0  !one year forcing repeat period
      endif
#if defined(RELO)
      if     (yrflag.lt.2) then
        natm = 4  !monthly forcing
      else
        natm = 2  !high-frequency forcing
      endif
#else
      if     (yrflag.lt.2) then
        if (natm.lt.4) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'error - natm must be 4 for monthly forcing (yrflag==0,1)'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      else
        if (natm.ne.2) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'warning - natm can be 2 for this forcing (yrflag!=0,1)'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
#endif
!
      dsur1p = dsurfq .lt. 0.0
      if     (dsur1p) then
        dsurfq = -dsurfq
      endif
!
      arcend = diagfq.lt.0.0
      diagfq = abs(diagfq)
!
      if     (cplifq.ge.0.0) then
        icpfrq = nint( cplifq*(86400.0/baclin) )
      else
        icpfrq = nint(-cplifq)
      endif
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,'(a,i10)') 'icpfrq =',icpfrq
      endif !1st tile
!
      if (yrflag.lt.0 .or. yrflag.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - yrflag must be between 0 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (yrflag.le.1) then  ! monthly forcing
        if     (abs(nint(86400.0/baclin)-86400.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') &
           &'error - baclin not an integer divisor of 24 hours'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        else  ! make it exact
!         write(lp,*) 'old baclin = ',baclin
          baclin = 86400.0/nint(86400.0/baclin)
!         write(lp,*) 'new baclin = ',baclin
        endif
      else  ! high frequency forcing
        if     (abs(nint(21600.0/baclin)-21600.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') &
           &'error - baclin not an integer divisor of 6 hours'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        else  ! make it exact
!         write(lp,*) 'old baclin = ',baclin
          baclin = 21600.0/nint(21600.0/baclin)
!         write(lp,*) 'new baclin = ',baclin
        endif
      endif
      if     (abs(nint(0.5*baclin/batrop)- &
                       0.5*baclin/batrop  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - batrop not an integer divisor of baclin/2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old batrop = ',batrop
        batrop = baclin/nint(baclin/batrop)
!       write(lp,*) 'new batrop = ',batrop
      endif
      if     (abs(nint((dsurfq*86400.d0)/baclin)- &
                       (dsurfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
         &'error - ', &
         &'dsurfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old dsurfq = ',dsurfq
        dsurfq = nint((dsurfq*86400.d0)/baclin)*(baclin/86400.d0)
!       write(lp,*) 'new dsurfq = ',dsurfq
      endif
      if     (abs(nint((diagfq*86400.d0)/baclin)- &
                       (diagfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
         &'error - ', &
         &'diagfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old diagfq = ',diagfq
        diagfq = nint((diagfq*86400.d0)/baclin)*(baclin/86400.d0)
!       write(lp,*) 'new diagfq = ',diagfq
      endif
      if     (abs(nint((proffq*86400.d0)/baclin)- &
                       (proffq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
         &'error - ', &
         &'proffq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old proffq = ',proffq
        proffq = nint((proffq*86400.d0)/baclin)*(baclin/86400.d0)
!       write(lp,*) 'new proffq = ',tilefq
      endif
      if     (abs(nint((tilefq*86400.d0)/baclin)- &
                       (tilefq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
         &'error - ', &
         &'tilefq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old tilefq = ',tilefq
        tilefq = nint((tilefq*86400.d0)/baclin)*(baclin/86400.d0)
!       write(lp,*) 'new tilefq = ',tilefq
      endif
      if     (abs(nint((meanfq*86400.d0)/baclin)- &
                       (meanfq*86400.d0)/baclin  ).gt.0.01) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
         &'error - ', &
         &'meanfq is not a whole number of baroclinic time steps'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      else  ! make it exact
!       write(lp,*) 'old meanfq = ',meanfq
        meanfq = nint((meanfq*86400.d0)/baclin)*(baclin/86400.d0)
!       write(lp,*) 'new meanfq = ',meanfq
      endif
#if defined(RELO)
      if (nestfq.ne.0.0) then
        kknest = kdm
      else
        kknest =   1
      endif
#else
      if (kknest.ne.kdm .and. nestfq.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - kknest (dimensions.h) must be kdm when nesting'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif
!
! --- 'incflg' = incremental update flag (0=no, -/+1=yes, -/+2=full-velocity)
!                  -ve to read in difference archive for increment
! --- 'incstp' = no. timesteps for full update (1=full insertion)
! --- 'incupf' = number of days of incremental updating input
! ---             (-ve to write a restart after each update completes)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(incflg, 'incflg')
      call blkini(incstp, 'incstp')
      call blkini(incupf, 'incupf')
!
! --- 'stfflg' = stochastic anomaly forcing flag (0=no, 1=TS, 2=TSV)
! --- 'stfrdt' = stochastic T anomaly forcing e-folding depth (m)
! --- 'stfrds' = stochastic S anomaly forcing e-folding depth (m)
! --- 'stfrdv' = stochastic V anomaly forcing e-folding depth (m)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(stfflg, 'stfflg')
      call blkinr(stfrdt, 'stfrdt','(a6," =",f10.4," m")')
      call blkinr(stfrds, 'stfrds','(a6," =",f10.4," m")')
      call blkinr(stfrdv, 'stfrdv','(a6," =",f10.4," m")')
!
! --- 'ra2fac' = weight for Robert-Asselin time filter
! ---             (old equivalent: wuv2=wts2=0.5*ra2fac,
! ---                              wuv2>wts2 no longer supported)
! --- 'wbaro ' = weight for time smoothing of barotropic fields
! --- 'btrlfr' = leapfrog barotropic time step (0=F,1=T)
! --- 'btrmas' = barotropic is mass conserving (0=F,1=T)
! --- 'hybraf' = HYBGEN:  Robert-Asselin flag  (0=F,1=T)
! ---             (use 'hybraf'=0 to recover pre-2.2.38 behaviour)
! --- 'hybrlx' = HYBGEN: inverse relaxation coefficient (time steps)
!                 (1.0 for no relaxation)
! --- 'hybiso' = HYBGEN: Use PCM if layer is within hybiso of target density
!                 (0.0 for no PCM; large to recover pre-2.2.09 behaviour)
! --- 'hybmap' = HYBGEN:  remapper  flag (0=PCM, 1=PLM,    2=PPM,  3=WENO-like)
! --- 'hybflg' = HYBGEN:  generator flag (0=T&S, 1=th&S,   2=th&T)
! --- 'advflg' = thermal  advection flag (0=T&S, 1=th&S,   2=th&T)
! --- 'advtyp' = scalar   advection type (0=PCM, 1=MPDATA, 2=FCT2, 4=FCT4)
! --- 'momtyp' = momentum advection type (2=2nd; 3=S.QUICK; 4=M.S.QUICK)
! --- 'slip'   = +1 free-slip, +1<slip<-1 partial, -1 no-slip boundary conditions
! --- 'visco2' = deformation-dependent Laplacian  viscosity factor
! --- 'visco4' = deformation-dependent biharmonic viscosity factor
! --- 'facdf4' =       speed-dependent biharmonic viscosity factor
! ---                   (diffusion velocity is facdf4*|v|)
! --- 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissipation
! ---             (negative to input spacially varying diffusion velocity)
! --- 'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissipation
! ---             (negative to input spacially varying diffusion velocity)
! --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
! --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion
! ---             (negative to input spacially varying diffusion velocity)
! --- 'temdf2' = diffusion velocity (m/s) for Laplacian  temp/saln diffusion
! --- 'temdfc' = temp diffusion conservation (0.0,1.0 all dens,temp resp.)
! --- 'vertmx' = diffusion velocity (m/s) for mom.mixing at mix.layr.base
! ---             (vertmx only used in MICOM-like isopycnal mode)
! --- 'cbar'   = rms flow speed     (m/s) for linear bottom friction
! ---             (negative to input spacially varying rms flow speed)
! --- 'cb'     = coefficient of quadratic bottom friction
! ---             (negative to input spacially varying coefficient)
! --- 'drglim' = limiter for explicit friction (1.0 no limiter, 0.0 implicit)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(ra2fac,'ra2fac','(a6," =",f10.4," ")')
      call blkinr(wbaro ,'wbaro ','(a6," =",f10.4," ")')
      call blkinl(btrlfr,'btrlfr')
      call blkinl(btrmas,'btrmas')
      call blkinl(hybraf,'hybraf')
      call blkinr(hybrlx,'hybrlx','(a6," =",f10.4," time steps")')
      call blkinr(hybiso,'hybiso','(a6," =",f10.4," kg/m^3")')
      call blkini(hybmap,'hybmap')
      call blkini(hybflg,'hybflg')
      call blkini(advflg,'advflg')
      call blkini(advtyp,'advtyp')
      call blkini(momtyp,'momtyp')
      call blkinr(slip,  'slip  ', &
           &'(a6," =",f10.4," (-1=no-slip, -1>:>+1=partial, +1=free)")')
      call blkinr(visco2,'visco2','(a6," =",f10.4," ")')
      call blkinr(visco4,'visco4','(a6," =",f10.4," ")')
      call blkinr(facdf4,'facdf4','(a6," =",f10.4," ")')
      call blkinr(veldf2,'veldf2','(a6," =",f10.4, &
                               &" m/s (-ve if variable)")')
      call blkinr(veldf4,'veldf4','(a6," =",f10.4, &
                               &" m/s (-ve if variable)")')
      call blkinr(thkdf2,'thkdf2','(a6," =",f10.4," m/s")')
      call blkinr(thkdf4,'thkdf4','(a6," =",f10.4, &
                               &" m/s (-ve if variable)")')
      call blkinr(temdf2,'temdf2','(a6," =",f10.4," m/s")')
      call blkinr(temdfc,'temdfc', &
         &'(a6," =",f10.4," (0.0,1.0 conserve dens,temp resp.)")')
      call blkinr(vertmx,'vertmx','(a6," =",f10.4," m/s")')
      call blkinr(cbar,  'cbar  ','(a6," =",f10.4, &
                               &" m/s (-ve if variable)")')
      call blkinr(cb,    'cb    ','(a6," =",f10.4, &
                               &"     (-ve if variable)")')
      call blkinr(drglim,'drglim','(a6," =",f10.4," ")')   
!
      if (hybrlx.lt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - hybrlx must be at least 1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      qhybrlx = 1.0/hybrlx
!
      if (btrmas .and. .not.btrlfr) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - btrmas==.true. requires btrlfr==.true.'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      isopcm = hybiso.gt.0.0  !use PCM for isopycnal layers?
      if (hybmap.lt.0 .or. hybmap.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')  &
         &'error - hybmap must be ', &
         &'0 (PCM) or 1 (PLM) or 2 (PPM) or 3 (WENO-like)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (hybflg.lt.0 .or. hybflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - hybflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (hybflg.ne.0 .and. (sigver.eq.5 .or. sigver.eq.6)) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - hybflg must be 0 (T&S) for 17-term eqn. of state'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.lt.0 .or. advflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - advflg must be 0 (T&S) or 1 (th&S) or 2 (th&T)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.ne.0 .and. (sigver.eq.5 .or. sigver.eq.6)) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - advflg must be 0 (T&S) for 17-term eqn. of state'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advflg.eq.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - advflg==2 (th&T) not yet implemented'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (advtyp.ne.0 .and. advtyp.ne.1 .and. &
          advtyp.ne.2 .and. advtyp.ne.4      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')  &
        &'error - advtyp must be 0,1,2,4', &
        &' (PCM,MPDATA,FCT2,FCT4)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (btrmas .and. advtyp.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)') &
        &'error - advtyp must be 2', &
        &' (FCT2) when btrmas==.true.'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (momtyp.lt.2 .or. momtyp.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
        &'error - momtyp must be 2 or 3 or 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (momtyp.eq.2) then
        if (facdf4.ne.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
          &'error - facdf4 must be 0.0 for montyp==2'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
      if (momtyp.eq.3) then
        if (abs(facdf4-0.0625).gt.1.e-5) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') &
          &'error - facdf4 must be 1/16 for momtyp==3'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
      if (momtyp.eq.4) then
        if (facdf4.le.0.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)') &
          &'error - facdf4 must be positive for momtyp==4'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      endif
      if (slip.lt.-1.0 .or. slip.gt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')  &
         &'error - slip must be ', &
         &'between -1.0 (no-slip) and +1.0 (free-slip)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (thkdf2.ne.0.0 .and. thkdf4.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - only one of thkdf2 and thkdf4 can be non-zero'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (isopyc .and. temdfc.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - isopycnal mode must have temdfc=0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (temdfc.lt.0.0 .or. temdfc.gt.1.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - temdfc must be between 0.0 and 1.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (temdfc.ne.1.0 .and. (sigver.eq.5 .or. sigver.eq.6)) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - temdfc must be 1.0 (all T) for 17-term eqn. of state'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
!
! --- 'lstep' = number of barotropic time steps per baroclinic time step.
! ---             lstep   m u s t   be even.
! --- 'dlt'   = barotropic time step
!
      lstep=nint(baclin/batrop)
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      if (mnproc.eq.1) then
      write (lp,'(i4,a/)') &
        lstep,' barotropic steps per baroclinic time step'
      endif !1st tile
#if ! defined(RELO)
!
      if     (btrmas) then
        if (max_nsteps_batrop.lt.lstep) then
          if     (mnproc.eq.1) then
            write(lp,*) 'max_nsteps_batrop = ',max_nsteps_batrop
            write(lp,*) 'max_nsteps_batrop too small, must be',lstep
          endif
          call xcstop('(blkdat)')
                 stop   blkdat
        endif  !max_nsteps_batrop
      endif  !btrmas
#endif
!
! --- 'thkbot' = thickness of bottom boundary layer (m)
! --- 'sigjmp' = minimum density jump across interfaces  (kg/m**3)
! --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
! --- 'sefold' = e-folding time                  for SSS relaxation (days)
! ---             (use 'sefold'=30.0 to recover original behaviour;
! ---              -ve to input, with sefold(i,j)=0.0 for no relax at i,j)
! --- 'thkmls' = reference mixed-layer thickness for SSS relaxation (m)
! --- 'thkmlt' = reference mixed-layer thickness for SST relaxation (m)
! --- 'thkriv' = nominal thickness of river inflow (m)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(thkbot,'thkbot','(a6," =",f10.4," m")')
      call blkinr(sigjmp,'sigjmp','(a6," =",f10.4," kg/m**3")')
      call blkinr(tmljmp,'tmljmp','(a6," =",f10.4," degC")')
      call blkinr(sefold,'sefold','(a6," =",f10.4, &
                               & " days (-ve if variable)")')
      call blkinr(thkmls,'thkmls','(a6," =",f10.4," m")')
      call blkinr(thkmlt,'thkmlt','(a6," =",f10.4," m")')
      call blkinr(thkriv,'thkriv','(a6," =",f10.4," m")')
!
      if (sefold.eq.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - sefold may be positive or negative, but not 0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
! --- 'thkcdw' = thickness for near-surface currents in ice-ocean stress (m)
! --- 'thkfrz' = maximum thickness of near-surface freezing zone (m)
! --- 'iceflg' = sea ice model flag (0=none,1=energy loan,1>coupled/esmf)
! ---             2=2-way,3=no IO stress, 4=3+no ocean currents
! ---             also, icmflg=3 for ENLN relaxed to coupler ice concentration
! --- 'tfrz_0' = ice melting point (degC) at S=0psu
! --- 'tfrz_s' = gradient of ice melting point (degC/psu)
! --- 'frzifq' = e-folding time scale back to tfrz
! ---             (positive days or negative time steps)
! --- 'ticegr' = ENLN: vertical temperature gradient inside ice (deg/m)
! ---                   (0.0 to get ice surface temp. from atmos. surtmp)
! --- 'hicemn' = ENLN: minimum ice thickness (m)
! --- 'hicemx' = ENLN: maximum ice thickness (m)
! --- 'ishelf' = ice shelf flag (0=none,1=ice shelf over ocean)
!
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinr(thkcdw,'thkcdw','(a6," =",f10.4," m")')
      call blkinr(thkfrz,'thkfrz','(a6," =",f10.4," m")')
      call blkini(iceflg,'iceflg')
      icegln = iceflg.eq.1  !ENLN, but see icmflg.eq.3 below
!
      if (iceflg.lt.0 .or. iceflg.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - iceflg must be between 0 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
! --- to recover original ENLN model use:
! --- tfrz_0=-1.8, tfrz_s=0.0,   ticegr=2.0, hicemn=0.5, hicemx=10.0.
! --- for freezing point linear in S, typically use:
! --- tfrz_0=0.0, tfrz_s=-0.054, ticegr=0.0, hicemn=0.5, hicemx=10.0.
!
      call blkinr(tfrz_0,'tfrz_0','(a6," =",f10.4," degC")')
      call blkinr(tfrz_s,'tfrz_s','(a6," =",f10.4," degC/psu")')
      call blkinr(frzifq,'frzifq','(a6," =",f10.4, &
                                &" days (-ve time steps)")')
      call blkinr(ticegr,'ticegr','(a6," =",f10.4," degC/m")')
      call blkinr(hicemn,'hicemn','(a6," =",f10.4," m")')
      call blkinr(hicemx,'hicemx','(a6," =",f10.4," m")')
!
      if     (frzifq.ge.0.0) then
        icefrq = nint( frzifq*(86400.0/baclin) )
      else
        icefrq = nint(-frzifq)
      endif
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,'(a,i10)') 'icefrq =',icefrq
      endif !1st tile
!
      call blkini(ishelf,'ishelf')
!
! --- 'ntracr' = number of tracers (0=none,negative to initialize)
! --- 'trcflg' = tracer flags      (one digit per tracer, most sig. replicated)
! ---              0: passive, 100% at surface
! ---              1: passive, psudo-silicate
! ---              2: passive, temperature
! ---              3: passive
! ---            4-8: unused
! ---              9: default biology (NPZ-3,NPZD-4,Chai-9)
!
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(ntracr,'ntracr')
      trcrin = ntracr.gt.0  ! positive from restart, otherwise initialize
      trcout = ntracr.ne.0
      ntracr = abs(ntracr)
!
      if (ntracr.gt.mxtrcr) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,i3, a /)')  &
         &'error - maximum ntracr is',mxtrcr, &
         &'  (recompile with larger mxtrcr)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      call blkini(trcflg1,'trcflg')
      do k= 1,ntracr
        trcflg(k) = mod(trcflg1,10)  ! least significant decimal digit
        if     (trcflg1.ge.10) then
          trcflg1 = trcflg1/10  ! shift by one decimal digit
        else
          ! replicate last decimal digit across remaining tracers
        endif
        if (mnproc.eq.1) then
        write(lp,'(a,i3,i2)') '    k,trcflg =',k,trcflg(k)
        endif !1st tile
        if     (trcflg(k).gt.3 .and. trcflg(k).lt.9) then !not 0,1,2,3,9
          if (mnproc.eq.1) then
          write(lp,'(/ a,i3 /)')  &
           &'error - unknown tracer type for tracer',k
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
      enddo !k
!
! --- 'tsofrq' = number of time steps between anti-drift offset calcs
! --- 'tofset' = temperature anti-drift offset (degC/century)
! --- 'sofset' = salnity     anti-drift offset  (psu/century)
!
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(tsofrq,'tsofrq')
      call blkinr(tofset,'tofset','(a6," =",f10.4," degC/century")')
      call blkinr(sofset,'sofset','(a6," =",f10.4,"  psu/century")')
!
! --- convert from per century to per second.
      tofset = tofset / (100.d0*365.d0*86400.d0)
      sofset = sofset / (100.d0*365.d0*86400.d0)
!
! --- 'mlflag' = mixed layer flag  (0=none,1=KPP,2-3=KTa-b,4=PWP,5=MY2.5,6=GISS)
! ---    'mxlkpp' = KPP:     activate mixed layer model (mlflag==1)
! ---    'mxlkrt' = KT:      MICOM or HYCOM Kraus-Turner (mlflag==2,3)
! ---    'mxlkta' = KT:      activate    original mixed layer model (mlflag==2)
! ---    'mxlktb' = KT:      activate alternative mixed layer model (mlflag==3)
! ---    'mxlpwp' = PWP:     activate mixed layer model (mlflag==4)
! ---    'mxlmy ' = MY:      activate mixed layer model (mlflag==5)
! ---   'mxlgiss' = GISS:    activate mixed layer model (mlflag==6)
! --- 'pensol' = KT/PWP   activate penetrating solar radiation (0=F,1=T)
! --- 'dtrate' = KT:      maximum permitted m.l. detrainment rate (m/day)
! --- 'thkmin' = KT/PWP:  minimum mixed-layer thickness (m)
! --- 'dypflg' = KT/PWP:  diapycnal mixing flag (0=none,1=KPP,2=explicit)
! ---                      always none (0) for MY/GISS and explicit (2) for PWP
! --- 'mixfrq' = KT/PWP:  number of time steps between diapycnal mixing calcs
! --- 'diapyc' = KT/PWP:  diapycnal diffusivity x buoyancy freq. (m**2/s**2)
! --- 'rigr'   = PWP:     critical gradient richardson number
! --- 'ribc'   = PWP:     critical bulk     richardson number
! --- 'rinfty' = KPP:     maximum  gradient richardson number (shear inst.)
! --- 'ricr'   = KPP:     critical bulk     richardson number
! --- 'bldmin' = KPP:     minimum surface boundary layer thickness (m)
! --- 'bldmax' = KPROF:   maximum surface boundary layer thickness (m)
! --- 'cekman' = KPP/KT:  scale factor for Ekman depth
! --- 'cmonob' = KPP:     scale factor for Monin-Obukov depth
! --- 'bblkpp' = KPP:     activate bottom boundary layer    (0=F,1=T)
! --- 'shinst' = KPP:     activate shear instability mixing (0=F,1=T)
! --- 'dbdiff' = KPP:     activate double diffusion  mixing (0=F,1=T)
! --- 'nonloc' = KPP:     activate nonlocal b. layer mixing (0=F,1=T)
! --- 'botdiw' = KPROF:   activate bot.enhan.int.wav mixing (0=F,1=T)
! --- 'difout' = KPROF:   output visc/diff coffs in archive (0=F,1=T)
! --- 'difsmo' = KPROF:   number of layers with horiz smooth diff coeffs
! --- 'difm0'  = KPP:     max viscosity   due to shear instability (m**2/s)
! --- 'difs0'  = KPP:     max diffusivity due to shear instability (m**2/s)
! --- 'difmiw' = KPP/MY:  background/internal wave viscosity       (m**2/s)
! ---             (negative to input spacially varying viscosity)
! --- 'difsiw' = KPP/MY:  background/internal wave diffusivity     (m**2/s)
! ---             (negative to input spacially varying diffusivity)
! --- 'dsfmax' = KPP:     salt fingering diffusivity factor        (m**2/s)
! --- 'rrho0'  = KPP:     salt fingering rp=(alpha*delT)/(beta*delS)
! --- 'cs'     = KPP:     value for nonlocal flux term
! --- 'cstar'  = KPP:     value for nonlocal flux term
! --- 'cv'     = KPP:     buoyancy frequency ratio (0.0 to use a fn. of N)
! --- 'c11'    = KPP:     value for turb velocity scale
! --- 'hblflg' = KPP:     b. layer interpolation flag (0=con.,1=lin.,2=quad.)
! --- 'niter'  = KPP:     iterations for semi-implicit soln. (2 recomended)
! --- 'langmr' = KPP:     Langmuir flag (0=no;1=Sul;2=Smy;3=Har;4=Tak)
! ---                       0:None
! ---                       1:McWilliams-Sullivan
! ---                       2:Smyth
! ---                       3:McWilliams-Harcourt
! ---                       4:Takaya
!
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(mlflag,'mlflag')
!
      if (mlflag.lt.0 .or. mlflag.gt.6) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - mlflag must be between 0 and 6'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      mxl_no = mlflag.eq.0
!
      if (isopyc .and. mlflag.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - isopycnal mode requires KT mixed layer (mlflag=2)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      mxlkta = mlflag.eq.2 .and. hybrid
      mxlktb = mlflag.eq.3 .and. hybrid
      mxlkrt = mlflag.eq.2 .or. mlflag.eq.3
      call blkinl(pensol,'pensol')
      call blkinr(dtrate,'dtrate','(a6," =",f10.4," m/day")')
      call blkinr(thkmin,'thkmin','(a6," =",f10.4," m")')
      call blkini(dypflg,'dypflg')
      call blkini(mixfrq,'mixfrq')
      call blkinr(diapyc,'diapyc','(a6," =",f10.4," m**2/s**2")')
!
      if (isopyc .and. pensol) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - isopycnal mode not consistent with pensol'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (dypflg.lt.0 .or. dypflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - dypflg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if ((mxlkta .or. mxlktb) .and. &
          dypflg.ne.0 .and. thkriv.gt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - thkriv>0 not consistent with KT unless dypflg=1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      mxlmy = mlflag.eq.5
!
      if (mxlmy) then
#if ! defined(RELO)
        if (kkmy25.ne.kdm) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'error - kkmy25 (dimensions.h) must be kdm when mlflag==5'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
#endif
        if     (dypflg.ne.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'warning - dypflg reset to 0 for M-Y 2.5'
          call flush(lp)
          endif !1st tile
          dypflg = 0
        endif
      endif
!
      mxlgiss = mlflag.eq.6
!
      if (mxlgiss) then
#if ! defined(RELO)
        if (nlgiss.ne.762) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'error - nlgiss (dimensions.h) must be 762 when mlflag==6'
          call flush(lp)
          endif !1st tile
          call xcstop('(blkdat)')
                 stop '(blkdat)'
        endif
#endif
        if     (dypflg.ne.0) then
          if (mnproc.eq.1) then
          write(lp,'(/ a /)')  &
           &'warning - dypflg reset to 0 for GISS'
          call flush(lp)
          endif !1st tile
          dypflg = 0
        endif
      endif
!
      mxlpwp = mlflag.eq.4
      call blkinr(rigc  ,'rigr  ','(a6," =",f10.4," ")')
      call blkinr(ribc  ,'ribc  ','(a6," =",f10.4," ")')
      call blkinr(qrinfy,'rinfty','(a6," =",f10.4," ")')
      qrinfy = 1.0/qrinfy
      call blkinr(ricr  ,'ricr  ','(a6," =",f10.4," ")')
!
      if (mxlpwp .and. dypflg.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'warning - dypflg  reset to 2 for PWP mixed layer'
        call flush(lp)
        endif !1st tile
        dypflg = 2
      endif
!
      if (mxlpwp .and. thkriv.gt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - thkriv>0 not consistent with PWP mixed layer'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      call blkinr(bldmin,'bldmin','(a6," =",f10.4," m")')
      call blkinr(bldmax,'bldmax','(a6," =",f10.4," m")')
!
      call blkinr(cekman,'cekman','(a6," =",f10.4," ")')
      call blkinr(cmonob,'cmonob','(a6," =",f10.4," ")')
!
      mxlkpp = mlflag.eq.1
      call blkinl(bblkpp,'bblkpp')
      call blkinl(shinst,'shinst')
      call blkinl(dbdiff,'dbdiff')
      call blkinl(nonloc,'nonloc')
      call blkinl(botdiw,'botdiw')
      call blkinl(difout,'difout')
      call blkini(difsmo,'difsmo') !now an integer
      call blkinr(difm0 ,'difm0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difs0 ,'difs0 ','(a6," =",f10.4," m**2/s")')
      call blkinr(difmiw,'difmiw','(a6," =",f10.4, &
                               &" m**2/s (-ve if variable)")')
      call blkinr(difsiw,'difsiw','(a6," =",f10.4, &
                               &" m**2/s (-ve if variable)")')
      call blkinr(dsfmax,'dsfmax','(a6," =",f10.4," m**2/s")')
      call blkinr(rrho0 ,'rrho0 ','(a6," =",f10.4," ")')
      call blkinr(cs    ,'cs    ','(a6," =",f10.4," ")')
      call blkinr(cstar ,'cstar ','(a6," =",f10.4," ")')
      call blkinr(cv    ,'cv    ','(a6," =",f10.4," ")')
      call blkinr(c11   ,'c11   ','(a6," =",f10.4," ")')
      call blkini(hblflg,'hblflg')
      call blkini(niter ,'niter ')
      call blkini(lngblk,'langmr')  !lngblk is local for testing 0 <= 4
!
#if defined(STOKES)
      langmr = lngblk
#else
      if     (lngblk.ne.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ 2a /)')  &
         &'error - langmr must be =0 unless the macro STOKES', &
         &' is defined at compile time'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif   
      if (difsiw*difmiw.lt.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ 2a /)')  &
         &'error - difmiw and difsiw must either both be positive or', &
         &' both be negative'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (hblflg.lt.0 .or. hblflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ 2a /)')  &
         &'error - hblflg must be', &
         &' 0 (constant) or 1 (linear) or 2 (quadratic)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (mxlkpp) then
! ---   for KPP, diapyc and vertmx are not used
        dypflg = 0
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlmy) then
! ---   for M-Y, diapyc and vertmx are not used
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlgiss) then
! ---   for GISS, diapyc and vertmx are not used
        diapyc = 0.0
        vertmx = 0.0
      endif
      if (mxlpwp) then
! ---   for PWP, vertmx is not used
        vertmx = 0.0
      endif
      if (mxlkta .or. mxlktb) then
! ---   for HYCOM KT, vertmx is not currently used
        vertmx = 0.0
      endif
!
      if     (dypflg.eq.2) then
        if     (max(tofset,sofset).gt.0.0) then
          if     (diapyc.le.0.0) then
            if (mnproc.eq.1) then
            write(lp,'(/ a /)')  &
             &'error - diapyc must be positive if [ts]ofset is'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          elseif (mixfrq.ne.tsofrq) then
            if (mnproc.eq.1) then
            write(lp,'(/ a /)')  &
             &'error - mixfrq and tsofrq must be equal'
            call flush(lp)
            endif !1st tile
            call xcstop('(blkdat)')
                   stop '(blkdat)'
          endif !error checks
        endif ![ts]ofset>0.0
      endif !dypflg.eq.2
!
! --- 'fltflg' = FLOATS: synthetic float flag (0=no; 1=yes)
! --- 'nfladv' = FLOATS: advect every nfladv bacl. time steps (even, >=4)
! --- 'nflsam' = FLOATS: output (0=every nfladv steps; >0=# of days)
! --- 'intpfl' = FLOATS: horiz. interp. (0=2nd order+bilinear; 1=bilinear)
! --- 'iturbv' = FLOATS: add horiz. turb. advection velocity (0=no; 1=yes)
! --- 'ismpfl' = FLOATS: sample water properties at float (0=no; 1=yes)
! --- 'tbvar'  = FLOATS: horizontal turbulent velocity variance scale
! --- 'tdecri' = FLOATS: inverse decorrelation time scale
!
      call blkini(fltflg,'fltflg')
      synflt = fltflg.ge.1
      call blkini(nfladv,'nfladv')
      call blkini(nflsam,'nflsam')
      call blkini(intpfl,'intpfl')
      call blkini(iturbv,'iturbv')
      turbvel = iturbv.ge.1
      call blkini(ismpfl,'ismpfl')
      samplfl = ismpfl.ge.1
      call blkinr(tbvar ,'tbvar ','(a6," =",f10.4," m**2/s**2")')
      call blkinr(tdecri,'tdecri','(a6," =",f10.4," 1/day")')
!
! --- 'lbflag' = lateral baro. bndy flag (0=none;nest:2=B-K,4=Flather)
! ---             (port: 1=Browning-Kreiss,3=Flather)
! --- 'lbmont' = baro nesting archives have sshflg=2
! ---             sshflg=2 is always ok, but is required if lbmont is set
! --- 'tidflg' = TIDES: tidal forcing flag    (0=no;1=bdy;2=body;3=bdy&body)
! --- 'tidein' = TIDES: tide field input flag (0=no;1=yes;2=sal)
! --- 'tidcon' = TIDES: 1 digit per constituent (Q1K2P1N2O1K1S2M2), 0=off,1=on
! --- 'tidsal' = TIDES: scalar self attraction and loading factor
! ---             (negative to input spacially varying SAL factor)
! --- 'tiddrg' = TIDES: tidal drag flag (0=no;1=scalar;2=tensor)
! --- 'thkdrg' = TIDES: thickness of bottom boundary layer for tidal drag (m)
! ---             (zero to apply tidal drag to barotropic mode)
! --- 'drgscl' = TIDES: scale factor for tidal drag
! --- 'tidgen' = TIDES: generic time (0=F,1=T)
! --- 'tidrmp' = TIDES:            ramp time (days)
! --- 'tid_t0' = TIDES: origin for ramp time (model day)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(lbflag,    'lbflag')
      call blkinl(lbmont,    'lbmont')
      call blkini(tidflg,    'tidflg')
      call blkini(tidein,    'tidein')
      call blkini(tidcon,    'tidcon')
      call blkinr(tidsal,    'tidsal','(a6," =",f10.4, &
                                   &"  (-ve if variable)")')
      call blkini(tiddrg,    'tiddrg')
      call blkinr(thkdrg,    'thkdrg','(a6," =",f10.4, &
                                   &" m (0 if barotropic)")')
      call blkinr(drgscl,    'drgscl','(a6," =",f10.4," ")')
      call blkinl(tidgen,    'tidgen')
      call blkinr(ramp_time ,'tidrmp','(a6," =",f10.4," days")')
      call blkin8(ramp_orig ,'tid_t0','(a6," =",f10.4," model day")')  !real*8
!
      if     (tidflg.lt.0 .or. tidflg.gt.3) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tidflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tidflg.ge.2 .and. .not.tidgen .and. yrflag.ne.3) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tidgen must be .true. for body tide and yrflag.ne.3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tidflg.gt.0) then
        if     (abs(nint(3600.0/baclin)-3600.0/baclin).gt.0.01) then
          if (mnproc.eq.1) then                                     
          write(lp,'(/ a /)')                                        &
           &'error - baclin not an integer divisor of 1 hour'
          call flush(lp)                                     
          endif !1st tile                                    
          call xcstop('(blkdat)')
                 stop '(blkdat)' 
        else  ! make it exact    
!         write(lp,*) 'old baclin = ',baclin
          baclin = 3600.0/nint(3600.0/baclin)
!         write(lp,*) 'new baclin = ',baclin 
        endif                                
      endif                                 
      if     (tidein.lt.0 .or. tidein.gt.2) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tidein must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tidein.gt.0 .and. tidflg.lt.2) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tidein implies input but there is no body tide'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tiddrg.lt.0 .or. tiddrg.gt.2) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tiddrg must be between 0 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tiddrg.eq.0 .and. drgscl.ne.0.0) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tiddrg=0 must have drgscl=0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tiddrg.ne.0 .and. drgscl.le.0.0) then
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tiddrg>0 must have drgscl>0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if     (tiddrg.eq.2 .and. thkdrg.ne.0.0) then
! ---   tensor drag only implemented on barotropic velocity
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - tiddrg=2 only implemented for thkdrg=0.0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
! --- 'clmflg' = climatology frequency flag (6=bimonthly,12=monthly)
! --- 'wndflg' = wind stress input flag (0=none,1=uv-grid,2,3=p-grid,4,5=wnd10m)
! ---             (=3 wind speed from wind stress; =4,5 wind stress from wind)
! ---             (=4 for COARE 3.0; =5 for COREv2 bulk parameterization)
! ---             (4,5 use relative wind U10-Uocn; -4,-5 use absolute wind U10)
! --- 'ustflg' = ustar forcing   flag          (3=input,1,2=wndspd,4=stress)
! --- 'flxflg' = thermal forcing flag   (0=none,3=net_flux,1-2,4-6=sst-based)
! ---             (=1 MICOM bulk parameterization)
! ---             (=2 Kara  bulk parameterization)
! ---             (=4 COARE bulk parameterization, approx.)
! ---             (=5 L&Y   bulk parameterization)
! ---             (=6 COARE bulk parameterization, approx., better pressure)
! --- 'empflg' = E-P     forcing flag   (0=none,3=net_E-P, 1-2,4-6=sst-based_E)
! ---             (flxflg to use model sst, -flxflg to use seatmp)
! --- 'emptgt' = E-P     balance target (mm/week, into ocean)
! --- 'empbal' = E-P     balance flag   (0=none,1=offset,2=scale)
! --- 'dswflg' = diurnal shortwv flag   (0=none,1=daily to diurnal correction)
! --- 'albflg' = ocean albedo    flag   (0=none,1=const,2=L&Y)
! --- 'sssflg' = SSS relaxation  flag   (0=none,1=clim,-1=clim&rmx)
! --- 'sssbal' = SSS rlx balance flag   (0=none,1=offset,2=scale)
! --- 'lwflag' = longwave corr.  flag   (0=none,1=clim,2=nwp,-1=lwdn), sst-based
! --- 'sstflg' = SST relaxation  flag   (0=none,1=clim,2=atmos,3=observed)
! --- 'icmflg' = ice mask        flag   (0=none,1=clim,2=atmos,3=obs/coupled)
! --- 'prsbas' = msl prs is input field + prsbas (Pa)
! --- 'mslprf' = msl prs forcing flag   (0=F,1=T)
! --- 'stroff' = net strs offset flag   (0=F,1=T)
! --- 'flxoff' = net flux offset flag   (0=F,1=T)
! --- 'flxsmo' = smooth surface fluxes  (0=F,1=T)
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkini(clmflg,'clmflg')
      call blkini(wndflg,'wndflg')
      if     (wndflg.lt.0) then
        wndflg = -wndflg
        amoflg = 0  !U10
      else
        amoflg = 1  !U10-Uocn
      endif
      call blkini(ustflg,'ustflg')
      call blkini(flxflg,'flxflg')
      call blkini(empflg,'empflg')
      call blkinr(emptgt,'emptgt','(a6," =",f10.4," mm/wk")')
      call blkini(empbal,'empbal')
      call blkini(dswflg,'dswflg')
      call blkini(albflg,'albflg')
      call blkini(sssflg,'sssflg')
      call blkini(sssbal,'sssbal')
      call blkini(lwflag,'lwflag')
      call blkini(sstflg,'sstflg')
      call blkini(icmflg,'icmflg')
      call blkinr(prsbas,'prsbas','(a6," =",f10.1," Pa")')
      call blkinl(mslprf,'mslprf')
      call blkinl(stroff,'stroff')
      call blkinl(flxoff,'flxoff')
      call blkinl(flxsmo,'flxsmo')
!
      emptgt = emptgt*rhoref/(7.d0*86400.d0*1000.d0)  !m/s kg/m^3 into ocean
!
      if (clmflg.ne.6 .and. clmflg.ne.12) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - clmflg must be 6 or 12'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#if defined(RELO)
      if (lbflag.eq.1 .or. lbflag.eq.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - lbflag 1 and 3 not supported with dynamic memory'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif
      if (lbflag.lt.0 .or. lbflag.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - lbflag must be between 0 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbflag.ne.2 .and. lbflag.ne.4 .and. bnstfq.ne.0.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - bnstfq must be 0.0 unless lbflag is 2 or 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lbmont .and. sshflg.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
          'error - sshflg must be 2 if baro nesting archives have this'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (wndflg.lt.0 .or. wndflg.gt.5) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - wndflg must be between 0 and 5'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (ustflg.lt.1 .or. ustflg.gt.4) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - ustflg must be between 1 and 4'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.lt.0 .or. flxflg.gt.6) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - flxflg must be between 0 and 6'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.gt.0 .and. &
          wndflg.eq.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - wndflg must be non-zero when flxflg>0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.6 .and. &
          wndflg.ne.4      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - wndflg must be 4 when flxflg==6'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#if ! defined(USE_NUOPC_CESMBETA)
      if (flxflg.eq.5 .and. &
          wndflg.ne.5      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - wndflg must be 5 when flxflg==5'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif
      if (flxflg.eq.3 .and.    &
          ustflg.eq.4      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - ustflg must be 1,2,3 when flxflg==3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif    
      if (abs(empflg).lt.0 .or. abs(empflg).gt.6) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - empflg must be between 0 and 6'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.3 .and. (empflg.ne.0 .and. empflg.ne.3)) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - empflg must be 0 or 3 when flxflg==3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (abs(empflg).ne.flxflg .and. empflg.ne.0 &
                                .and. empflg.ne.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - empflg must be 0 or 3 or +\-flxflg'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (empflg.lt.0 .and. sstflg.lt.2 .and. lwflag.ne.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')  &
         &'error - negative empflg requires sst forcing', &
         &' (i.e. sstflg=3 or sstflg=2 or lwflag=2)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.0 .and. &
          empflg.ne.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - empflg must be 0 when flxflg==0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (flxflg.eq.3 .and. &
          dswflg.eq.1      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - dswflg must be 0 when flxflg==3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (albflg.gt.0 .and. lwflag.ne.-1) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - albflg>0 requires lwflag=-1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lwflag.lt.-1 .or. lwflag.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - lwflag must be between -1 and 2'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (lwflag.ne.0 .and. &
          flxflg.eq.0      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - flxflg must non-zero when lwflag is non-zero'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (wndflg.eq.0) then
        stroff = .false. !no wind forcing
      endif
!
      if (flxflg.eq.0) then
        flxoff = .false. !no thermal forcing
      endif
!
      if (iceflg.eq.0) then
        icmflg = 0   !no ice
        ticegr = 2.0 !surtmp not needed for ice temperature
      elseif (iceflg.ge.2) then
        icegln = icmflg.eq.3  !ENLN plus relax to coupler ice concentration
      endif
!
      if (icmflg.lt.0 .or. icmflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - icmflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      elseif (icmflg.eq.1) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - icmflg 1 not yet implemented'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      elseif (icmflg.eq.3 .and. iceflg.eq.1) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
         &'error - icmflg 3 not yet implemented for iceflg==1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (icmflg.eq.3 .and. icefrq.lt.icpfrq) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - icefrq (frzifq) smaller than icpfrq (cplifq)'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (sstflg.lt.0 .or. sstflg.gt.3) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - sstflg must be between 0 and 3'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (sstflg.gt.1 .and. &
          yrflag.lt.2      ) then
        if (mnproc.eq.1) then
        write(lp,'(/ a,a /)')  &
         &'error - yrflag must be >1 (high frequency forcing)', &
         &' when sstflg>1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if (sssflg.lt.-1 .or. sssflg.gt.2) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - sssflg must be 0 or 1 or -1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
! ---  s w i t c h e s    (if set to .true., then...)
! ---  (due to a SGI bug: read in an integer with 0=F,1=T)
! --- windf       use wind stress forcing (wndflg>0)
! --- thermo      use thermodynamic forcing 
! --- pensol      use penetrating solar radiation (input above)
! --- pcipf       use E-P forcing (may be redefined in forfun)
! --- priver      rivers as a precipitation bogas
! --- epmass      treat evap-precip as a mass exchange
!
! --- srelax      activate surface salinity        climatological nudging
! --- trelax      activate surface temperature     climatological nudging
! --- relax       activate lateral boundary T/S/p  climatological nudging
! --- trcrlx      activate lateral boundary tracer climatological nudging
! --- relaxt      input tracer  climatological relaxation fields
! --- relaxf      input T/S/p   climatological relaxation fields
! ---              (note that relaxt implies relaxf)
! --- relaxs      input surface climatological relaxation fields only
!
      windf  = wndflg.ne.0
      thermo = flxflg.ne.0
      pensol = pensol .and. thermo
      pcipf  = empflg.ne.0  ! if .true., might later be set .false. by forfun
!
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
      call blkinl(relax, 'relax ')
      call blkinl(trcrlx,'trcrlx')
      call blkinl(priver,'priver')
      call blkinl(epmass,'epmass')
      if (mnproc.eq.1) then
      write(lp,*)
      endif !1st tile
!
      if     (priver .and. .not.thermo) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - priver must be .false. for flxflg=0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      if     (epmass .and. .not.thermo) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - epmass must be .false. for flxflg=0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!!Alex add condition epmass.and.btrlfr
      if     (epmass .and. .not.btrlfr) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') &
          'error - btrlfr must be .true. for epmass=1'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
!
      srelax = sssflg.eq.1 .or. sssflg.eq.-1
      trelax = sstflg.eq.1
!
      if (.not.(thermo .or. sstflg.gt.0 .or. srelax)) then
        niter=1
      elseif (mxlkta) then
        if (mnproc.eq.1) then
        write(lp,'(/ a / a /)')  &
         &'error - KT mixed layer needs thermal forcing, i.e.', &
         &'        mlflag=2 needs max(sstflg,sstflg,flxflg)>0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif !.not(thermo ...):mxlkta
!
      relaxt = trcout .and. (trcrlx .or. &
                             (.not.trcrin .and. iniflg.eq.2))
      relaxf =       relax &
                .or. srelax &
                .or. trelax &
                .or. lwflag.eq.1 &
                .or. relaxt
      relaxs = .not. relax &
               .and. relaxf
#if defined(STOKES)
!
! ---  s w i t c h e s  for Stokes Drift effects  (0=F, 1=T)
!
! --- 'stdflg'  = STOKES: add Stokes Drift velocities to kinematics and dynamics
! --- 'stdsur'  = STOKES: add Stokes Drift Surface Stresses  to dynamics
! --- 'stdbot'  = STOKES: add Stokes Waves Field bottom friction drag
! --- 'stdarc'  = STOKES: Stokes Drift Velocities in Archive
!
! ---  s w i t c h e s  for Stokes Drift effects  (Integers)
!
! --- 'nsdzi '  = STOKES: number of interfaces in Stokes Drift input
!
! --- since they are at the end of blkdat.input they can optionally be
! --- present even if the macro /* STOKES */ is not defined at compile time.
!
      call blkinl(stdflg,'stdflg')
      call blkinl(stdsur,'stdsur') 
      call blkinl(stdbot,'stdbot')
      call blkinl(stdarc,'stdarc')
      call blkini(nsdzi ,'nsdzi ')
!
      if (nsdzi .lt.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
       &'error - nsdzi  must be 0 (no Stokes drift) or positive'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (nsdzi .eq.0 .and. stdflg) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
       &'error - nsdzi  must be > 0 with stdflg'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (nsdzi .eq.0 .and. stdbot) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
       &'error - nsdzi  must be > 0 with stdbot'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
      if (nsdzi .eq.0 .and. langmr.ne.0) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
       &'error - nsdzi  must be > 0 with langmr >0'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif   
!
      close(unit=uoff+99)  !file='blkdat.input'
!
! --- initialize from climatology (update relaxf and relaxs)?
      if     (iniflg.eq.2) then
#if defined (USE_NUOPC_CESMBETA)
          day1 = hycom_start_dtg
#else
        open(unit=uoff+99,file=trim(flnminp)//'limits')  !on all nodes
        read(uoff+99,*) day1
        close(unit=uoff+99) !file='limits'
#endif /* USE_NUOPC_CESMBETA */
        if     (day1.le.0.0) then
          relaxf = .true.
          relaxs = .false.
        endif !initialize from climatology
      endif
!
#if defined(RELO)
      if (relaxf .and. .not. relaxs) then
        kkwall = kdm
      elseif (relaxt) then
        kkwall = kdm  !needed for tracers
      else
        kkwall =   1
      endif
      if (mnproc.eq.1) then
      write(lp,'(a,i10)') 'kkwall =',kkwall
      write(lp,*)
      endif !1st tile
#else
      if (kkwall.ne.kdm .and. relaxf .and. .not. relaxs) then
        if (mnproc.eq.1) then
        write(lp,'(/ a /)')  &
         &'error - kkwall (dimensions.h) must be kdm for 3-d clim'
        call flush(lp)
        endif !1st tile
        call xcstop('(blkdat)')
               stop '(blkdat)'
      endif
#endif
#if defined(OCEANS2)
      if     (nocean.eq.2) then
! ---   slave HYCOM works from ./OCEAN2
!
! ---   i/o file names
!
        flnmshlf = './OCEAN2/regional.iceshelf'
        flnmdep  = './OCEAN2/regional.depth'
        flnmgrd  = './OCEAN2/regional.grid'
        flnmarc  = './OCEAN2/archv.0000_000_00'
        flnmarcs = './OCEAN2/archs.0000_000_00'
        flnmarcm = './OCEAN2/archm.0000_000_00'
        flnmarct = './OCEAN2/archt.0000_000_00'
        flnmovr  = './OCEAN2/ovrtn_out'
        flnmflx  = './OCEAN2/flxdp_out'
        flnmrsi  = './OCEAN2/restart_in'
        flnmrso  = './OCEAN2/restart_out'
!
! --- i/o directory names
!
        flnmfor  = './OCEAN2/'
        flnmforw = './OCEAN2/'
      else
! ---   master HYCOM
!
! ---   i/o file names
!
        flnmshlf = 'regional.iceshelf'
        flnmdep  = 'regional.depth'
        flnmgrd  = 'regional.grid'
        flnmarc  = 'archv.0000_000_00'
        flnmarcs = 'archs.0000_000_00'
        flnmarcm = 'archm.0000_000_00'
        flnmarct = 'archt.0000_000_00'
        flnmovr  = 'ovrtn_out'
        flnmflx  = 'flxdp_out'
        flnmrsi  = 'restart_in'
        flnmrso  = 'restart_out'
!
! --- i/o directory names
!
        flnmfor  = './'
        flnmforw = './'
      endif
#else
!
! --- i/o file names
!
      flnmshlf = 'regional.iceshelf'
      flnmdep  = 'regional.depth'
      flnmgrd  = 'regional.grid'
      flnmarc  = 'archv.0000_000_00'
      flnmarcs = 'archs.0000_000_00'
      flnmarcm = 'archm.0000_000_00'
      flnmarct = 'archt.0000_000_00'
      flnmovr  = 'ovrtn_out'
      flnmflx  = 'flxdp_out'
      flnmrsi  = 'restart_in'
      flnmrso  = 'restart_out'
!
! --- i/o directory names
!
      flnmfor  = './'
      flnmforw = './'
#endif /* OCEANS2:else */
      return
      end
      subroutine blkinr(rvar,cvar,cfmt)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      real      rvar
      character cvar*6,cfmt*(*)
!
!     read in one real value
!
      character*6 cvarin
!
      read(uoff+99,*) rvar,cvarin
      if (mnproc.eq.1) then
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
      endif !1st tile
!
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinr - input ',cvarin, &
                            &' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinr)')
               stop
      endif
      return
      end
      subroutine blkinr2(rvar,nvar,cvar1,cfmt1,cvar2,cfmt2)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      real      rvar
      integer   nvar
      character cvar1*6,cvar2*6,cfmt1*(*),cfmt2*(*)
!
!     read in one real value
!     identified as either cvar1 (return nvar=1) or cvar2 (return nvar=2)
!
      character*6 cvarin
!
      read(uoff+99,*) rvar,cvarin
      if     (cvar1.eq.cvarin) then
        nvar = 1
        if (mnproc.eq.1) then
        write(lp,cfmt1) cvarin,rvar
        call flush(lp)
        endif !1st tile
      elseif (cvar2.eq.cvarin) then
        nvar = 2
        if (mnproc.eq.1) then
        write(lp,cfmt2) cvarin,rvar
        call flush(lp)
        endif !1st tile
      else
        if (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in blkinr2 - input ',cvarin, &
                           &' but should be ',cvar1,' or ',cvar2
        write(lp,*)
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinr2)')
               stop
      endif
      return
      end
      subroutine blkin8(rvar,cvar,cfmt)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      real*8    rvar
      character cvar*6,cfmt*(*)
!
!     read in one real*8 value
!
      character*6 cvarin
!
      read(uoff+99,*) rvar,cvarin
      if (mnproc.eq.1) then
      write(lp,cfmt) cvarin,rvar
      call flush(lp)
      endif !1st tile
!
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkin8 - input ',cvarin, &
                            &' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkin8)')
               stop
      endif
      return
      end
      subroutine blkini(ivar,cvar)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer     ivar
      character*6 cvar
!
!     read in one integer value
!
      character*6 cvarin
!
      read(uoff+99,*) ivar,cvarin
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,ivar
      call flush(lp)
      endif !1st tile
!
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkini - input ',cvarin, &
                            &' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkini)')
               stop
      endif
      return
 6000 format(a6,' =',i10)
      end
      subroutine blkinl(lvar,cvar)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      logical     lvar
      character*6 cvar
!
!     read in one logical value
!     due to a SGI bug for logical I/O: read in an integer 0=F,1=T
!
      character*6 cvarin
      integer     ivar
!
      read(uoff+99,*) ivar,cvarin
      lvar = ivar .ne. 0
      if (mnproc.eq.1) then
      write(lp,6000) cvarin,lvar
      call flush(lp)
      endif !1st tile
!
      if     (cvar.ne.cvarin) then
        if (mnproc.eq.1) then
        write(lp,*) 
        write(lp,*) 'error in blkinl - input ',cvarin, &
                            &' but should be ',cvar
        write(lp,*) 
        call flush(lp)
        endif !1st tile
        call xcstop('(blkinl)')
               stop
      endif
      return
 6000 format(a6,' =',l10)
      end
!>
!> Revision history
!>
!> Oct. 1999 - added variables to control penetrating solar radiation
!> Oct. 1999 - added switch to select mixed layer model
!> Oct. 1999 - dp00 for hybgen is now set here
!> Dec. 1999 - multiple heat flux transfer coefficients (cts1, cts2, ctl)
!>             replace old coefficient ct
!> Jan. 2000 - changed to subroutine with run-time input
!> May. 2000 - conversion to SI units
!> Nov. 2000 - added kapflg,thflag,hybflg,advflg,wndflg
!> Dec. 2000 - added flxflg
!> Aug. 2001 - added bnstfq,nestfq
!> Aug. 2001 - added mapflg==4 for an f-plane
!> Aug. 2001 - betard and betabl inverted to replace division by multiply
!> May  2002 - map projection via regional.grid file (see geopar.f)
!> May  2002 - vertical coordinate via d[ps]00*
!> May  2002 - diffusion variable names include 2/4 for Laplacian/biharmonic
!> May  2002 - added PWP and MY2.5 mixed layer options
!> May  2002 - added thkmlr, distinct from thkmin
!> Aug. 2002 - added ntracr and trcflg to control tracers
!> Nov. 2002 - added jerlv0=0 to use kpar-based turbidity
!> Nov. 2002 - split thkmlr into thkmls and thkmlt
!> Apr. 2003 - added dp00i, vsigma, and priver
!> May  2003 - added bldmin, bldmax, flxsmo, and icmflg
!> June 2003 - added locsig and removed thflag=4
!> Oct. 2003 - thkdf4 negative now signals spacial variation
!> Nov. 2003 - added advtyp
!> Jan. 2004 - added latdiw
!> Jan. 2004 - added bblkpp
!> Jan. 2004 - added hblflg and cv=0.0 option
!> Feb. 2004 - added botdiw (and latdiw) for GISS
!> Feb. 2004 - added temdfc
!> Mar. 2004 - added thkriv and epmass
!> Mar. 2004 - added isotop
!> Mar. 2005 - added tfrz_0,tfrz_s,ticegr,hicemn,hicemx
!> Mar. 2005 - added tsofrq,tofset,sofset
!> Mar. 2005 - added empflg
!> Mar. 2005 - added ustflg, reordered thermal forcing flags
!> Mar. 2005 - added flxoff
!> Apr. 2005 - replaced kapflg with kapref
!> June 2005 - added hybrlx
!> June 2006 - added cplifq, thkfrz and negative dsurfq option
!> Nov. 2006 - version 2.2
!> Nov. 2006 - added incflg,incstp,incupf
!> Nov. 2006 - added FLOATS (fltflg,...)
!> Nov. 2006 - added TIDES  (tidflg,tidrmp,tid_t0,lbflag==3)
!> Mar. 2007 - added drgscl
!> Mar. 2007 - added drglim
!> Mar. 2007 - added tidcon,tidsal,tidgen
!> Apr. 2007 - implemented meanfq 
!> Apr. 2007 - added btrlfr and btrmas (latter not yet implemented)
!> May  2007 - added wbaro
!> June 2007 - moved h1 to momtum
!> June 2007 - added momtyp and facdf4
!> Sep. 2007 - added hybmap and hybiso
!> Feb. 2008 - added thkdrg
!> Feb. 2008 - added sshflg
!> June 2008 - added tilefq
!> July 2008 - added hybmap=3 (WENO-like)
!> Oct. 2008 - added dp0k and ds0k input option
!> Oct. 2008 - added dswflg
!> Dec. 2008 - difsmo is now an integer number of layers
!> Jan. 2009 - added -ve diagfq (arcend)
!> June 2009 - added sssrmx
!> Mar. 2010 - added negative bnstfq and nestfq for mean nesting archives
!> Apr. 2010 - added sstflg=-1 for sssrmx array input (no sssrmx blkdat entry)
!> Apr. 2010 - removed latdiw
!> Apr. 2010 - added proffq
!> Oct. 2010 - added support for 17-term and 12-term equation of state
!> Nov. 2010 - added yrflag=-2
!> Nov. 2010 - empflg negative to use seatmp instead of hycom sst
!> Apr. 2011 - added negative cbar 
!> July 2011 - added negative tidsal
!> Aug. 2011 - added ra2fac, removed global wts[12] and wuv[12]
!> Aug. 2011 - added hybraf
!> Sep. 2011 - added negative cb
!> Nov. 2011 - iniflg=2 now active for yrflag=3
!> Nov. 2011 - added frzifq
!> Jan. 2012 - added thkcdw
!> Jan. 2012 - added lbflag=4
!> Mar. 2012 - new terrain following method, based on dpns and dsns
!> Apr. 2012 - added negative tidflg, to set tidef
!> May  2012 - added tidein, replaces negative tidflg
!> Nov. 2012 - added /* OCEANS2 */ macro for master/slave HYCOM
!> Nov. 2012 - added wndlfg=4 for 10m wind component input
!> Nov. 2012 - added stroff, usualy used with wndflg=4
!> Jan. 2013 - added tiddrg, set to 0 or 1 for backwards compatibility
!> Jan. 2013 - added zero thkdrg to apply tidal drag to barotropic mode
!> June 2013 - added lbflag=6
!> July 2013 - added negative (spacially varying) difmiw and difsiw
!> July 2013 - added 3 Stokes Drift Flags: stkflg,stksur,stkbot
!> Aug. 2013 - added langmr flag for KPP and Stokes Drift
!> Sep. 2013 - added number of interfaces in Stokes Drift Files 
!> Oct. 2013 - added jerlv0=-1 to use chlorophyll-based turbidity
!> Oct. 2013 - fixed srelax used before defined bug
!> Nov. 2013 - added wndlfg=5 for 10m wind component input, COREv2 stress
!> Nov. 2013 - added flxlfg=5 for COREv2 heat flux
!> Nov. 2013 - added lwflag=-1 for radflx=Qlwdn
!> Nov. 2013 - added albflg for ocean albedo, >0 for shwflx=Qswdn
!> Jan. 2014 - added prsbas and mslprf for sea level pressure forcing
!> Jan. 2014 - /consts/ replaced by parameters
!> Jan. 2014 - kdm    is defined here when macro /* RELO */ is set
!> Jan. 2014 - natm   is defined here when macro /* RELO */ is set
!> Jan. 2014 - kkwall is defined here when macro /* RELO */ is set
!> Jan. 2014 - kknest is defined here when macro /* RELO */ is set
!> Jan. 2014 - lbflag 1 and 3 not supported with dynamic memory allocation
!> Apr. 2014 - added ishelf and flnmshlf
!> May  2014 - removed lbflag=6
!> Sep. 2014 - langmr can now be 0 to 4
!> Oct. 2014 - added flxlfg=6, sea level pressure is input even if .not. mslprf
!> Aug. 2015 - added stdarc
!> Aug. 2015 - added mxl_no
!> Nov. 2015 - added iniflg<0 for reading a difference archive
!> July 2017 - added momtyp=3 for Split QUICK momentum advection
!> July 2017 - momtyp=4 now allows non-zero visco2 and veldf2
!> Aug. 2017 - incupf -ve to write a restart after each update completes
!> Sep. 2017 - added stfflg and stfrd[tsv]
!> Apr. 2018 - thkmls is -ve for region-wide balanced relaxation (now sssbal)
!> Apr. 2018 - wndflg=5 and flxflg=5 must be used together
!> Aug. 2018 - allow partial-slip (-1<slip<+1)
!> Aug. 2018 - allow btrmas to be true
!> Aug. 2018 - added nsteps_batrop
!> Nov. 2018 - added empbal and sssbal
!> Nov. 2018 - added sefold
!> Nov. 2018 - added emptgt
!> Dec. 2018 - add /* USE_NUOPC_CESMBETA */ macro for coupled simulation
!> Dec. 2018 - add yrflag=4 for 365 days no-leap calendar (CESM)
!> Feb. 2019 - add sshflg=2 for steric Montg. Potential
!> Mar. 2019 - updated iversn to 23
!> Sep. 2019 - added oneta0
!> Oct. 2019 - added lbmont
!> Nov. 2019 - added wndflg=-4,-5 and amoflg
