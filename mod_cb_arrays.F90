      module mod_cb_arrays
      use mod_dimensions
      implicit none
      public ! everything is public

#if defined (USE_NUOPC_CESMBETA) || (ESPC_COUPLE)
#define USE_NUOPC_GENERIC 1 
#endif

!
!     wrapper for common_blocks.h
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2) ::  &
#endif
       u,v,            & ! velocity components
       dp,dpo,         & ! layer thickness on p-grid
       dpu,dpv,        & ! layer thickness on u-grid and v-grid
       temp,           & ! temperature
       saln,           & ! salinity
       th3d,           & ! potential density
       thstar,         & ! virtual potential density
       montg             ! montgomery potential

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,mxtrcr) ::  &
#endif
       tracer         ! inert tracers

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::  &
#endif
       otemp,          & ! temperature,                 time level t-1
       osaln,          & ! salinity,                    time level t-1
       oth3d,          & ! potential density,           time level t-1
       oq2,            & ! tke                          time level t-1
       oq2l              ! tke*turbulent length scale,  time level t-1

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,mxtrcr) ::  &
#endif
       otracer        ! inert tracers, time level t-1

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       dpmold         ! mixed layer depth at old time step

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::  &
#endif
       p,pu,pv        ! interface pressure

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::  &
#endif
       theta,          & ! isopycnal layer target densties - thbase
       diaflx            ! time integral of diapyc.flux

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       corio,          & ! coriolis parameter
       srfhgt,         & ! sea surface height, g*ssh(m)
       steric,         & ! steric sea surface height, g*sssh(m)
       sshgmn,         & !   mean sea surface height, g*mssh(m)
       thmean,         & !   mean depth averaged density
       montg_c,        & ! Montgomery Potential correction (m)
       montg1,         & ! layer 1 montgomery potential
       skap              ! thermobaric scale factor between reference states

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  &
#endif
       psikk,          & ! montg.pot. in bottom layer
       thkk              ! virtual potential density in bottom layer
!                                                                   
#if defined(RELO)
      integer, save, allocatable, dimension(:,:) ::  &
#else
      integer, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       kapi           ! thermobaric reference state index (1 or 3)
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::  &
#endif
       uflx,vflx,      & ! mass fluxes
       uflxav,vflxav,  & ! average fluxes
       dpav              ! average fluxes

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) ::  &
#endif
       ubavg,vbavg,    & ! barotropic velocity
       pbavg             ! barotropic pressure

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  &
#endif
       oneta,          & ! 1+eta
       onetao,         & ! 1+eta, original for Robert-Asselin filter
       onetamas          ! either 1+eta or 1, based on btrmas

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       ubrhs, vbrhs,   & ! rhs of barotropic u,v eqns.
       utotm, vtotm,   & ! total (barotrop.+baroclin.)..
       utotn, vtotn,   & ! ..velocities at 2 time levels
       uflux, vflux,   & ! horizontal mass fluxes
       uflux2,vflux2,  & ! more mass fluxes
       uflux3,vflux3,  & ! more mass fluxes
       onetacnt          ! 1+eta(t+1) after cnuity
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       util1,util2,    & ! arrays for temporary storage
       util3,util4,    & ! arrays for temporary storage
       util5,util6,    & ! arrays for temporary storage
       plon, plat,     & ! lon,lat at p pts
       qlon, qlat,     & ! lon,lat at q pts
       ulon, ulat,     & ! lon,lat at u pts
       vlon, vlat,     & ! lon,lat at v pts
       pang,           & ! angle between xwards and ewards
       scux, scuy,     & ! mesh size at u pts in x,y dir.
       scvx, scvy,     & ! mesh size at v pts in x,y dir.
       scpx, scpy,     & ! mesh size at p pts in x,y dir.
       scqx, scqy,     & ! mesh size at q pts in x,y dir.
       scu2, scv2,     & ! grid box area at u,v pts
       scp2, scq2,     & ! grid box area at p,q pts
       scp2i,scq2i,    & ! inverses of scp2,scq2
       scuxi,scvyi,    & ! inverses of scux,scvy
       aspux,aspuy,    & ! u-grid aspect ratios for diffusion
       aspvx,aspvy,    & ! v-grid aspect ratios for diffusion
       veldf2u,        & ! u-grid laplacian  diffusion coefficient
       veldf2v,        & ! v-grid laplacian  diffusion coefficient
       veldf4u,        & ! u-grid biharmonic diffusion coefficient
       veldf4v,        & ! v-grid biharmonic diffusion coefficient
       thkdf4u,        & ! u-grid biharmonic diffusion coefficient
       thkdf4v,        & ! v-grid biharmonic diffusion coefficient
       cbp,            & ! p-grid quadratic bottom friction coefficient
       cbarp,          & ! p-grid rms flow speed for linear bottom friction
       pgfx, pgfy,     & ! horiz. presssure gradient
       gradx,grady,    & ! horiz. presssure gradient
       depthu,depthv,  & ! bottom pres. at u,v points
       pvtrop,         & ! pot.vort. of barotropic flow
       depths,         & ! water depth
       drag,           & ! bottom drag
       salfac,         & ! spatialy varying "scalar" SAL factor
       topiso,         & ! shallowest depth for isopycnal layers (pressure units)
       diws,           & ! spacially varying background/internal wave diffusivity
       diwm,           & ! spacially varying background/internal wave viscosity
       diwbot,         & ! background/internal wave diffusivity at the bottom
       diwqh0,         & ! background/internal wave diffusivity vertical scale
       sssrmx,         & ! maximum SSS difference for relaxation (psu)
       tidepg_mn,      & ! tidal pressure gradient forcing, time mean
       displd_mn,      & ! dissipation from linear    drag, time mean
       dispqd_mn         ! dissipation from quadratic drag, time mean

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, dimension(1:2,1:2,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       drgten         ! tidal bottom drag tensor
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       uja,   ujb,     & ! velocities at lateral ..
       via,   vib,     & !       .. neighbor points
       pbot,           & ! bottom pressure at t=0
       sgain,          & ! salin.changes from diapyc.mix.
       surtx,          & ! surface net x-stress on p-grid
       surty,          & ! surface net y-stress on p-grid
       surflx,         & ! surface net thermal energy flux
       sswflx,         & ! surface swv thermal energy flux
       mixflx,         & ! mixed layer thermal energy flux
       sstflx,         & ! surface thermal flux from sst relax
       sssflx,         & ! surface salt    flux from sss relax
       rivflx,         & ! surface water   flux from rivers
       salflx,         & ! surface salt    flux
       wtrflx,         & ! surface water   flux
       buoflx,         & ! mixed layer buoyancy flux
       bhtflx,         & ! mixed layer buoyancy flux from heat
       wndocn,         & ! magnitude of 10m wind minus ocean current
       ustar,          & ! friction velocity
       ustarb,         & ! bottom friction velocity
       turgen,         & ! turb.kin.energ. generation
       thkice,         & ! grid-cell avg. ice thknss (m)
       covice,         & ! ice coverage (rel.units)
       temice,         & ! ice surface temperature
       flxice,         & ! heat  flux under ice
       fswice,         & ! swv   flux under ice
       wflice,         & ! water flux under ice
       wflfrz,         & ! water flux under ice from freeze/melt (diagnostic, icloan)
       sflice,         & ! salt  flux under ice
         si_c,         & ! ice concentration   on p-grid from coupler
         si_h,         & ! ice thickness       on p-grid from coupler
         si_t,         & ! ice temperature     on p-grid from coupler
         si_u,         & ! ice u-velocity      on p-grid from coupler
         si_v,         & ! ice v-velocity      on p-grid from coupler
        si_tx,         & ! x-stress  under ice on p-grid from coupler
        si_ty            ! y-stesss  under ice on p-grid from coupler
!
#if defined(RELO)
      integer, save, allocatable, dimension(:,:) ::  &
#else
      integer, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       klist,          & ! k-index
       jerlov            ! jerlov water type 1-5, 0 for kpar, -1 for chl
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  &
#endif
       dpmixl,         & ! mixed layer depth
       t1sav,          & ! upper sublayer temperature
       s1sav,          & ! upper sublayer salinity
       tmlb,           & ! temp in lyr. containing mlb.
       smlb              ! saln in lyr. containing mlb

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       hekman,         & ! ekman layer thickness
       hmonob,         & ! monin-obukhov length
       dpbl,           & ! turbulent boundary layer depth
       dpbbl,          & ! bottom turbulent boundary layer depth
       tmix,           & ! mixed layer temperature
       smix,           & ! mixed layer salinity
       thmix,          & ! mixed layer potential density
       umix,  vmix       ! mixed layer velocity

#if defined(RELO)
      real, save, allocatable, dimension(:) ::  &
#else
      real, save, dimension(kdm) ::  &
#endif
       dp0k,           & ! minimum deep    z-layer separation
       ds0k              ! minimum shallow z-layer separation

      real, save :: &
       dpns,           & ! depth to start terrain following
       dsns              ! depth to stop  terrain following
!
#if defined(RELO)
      integer, save, allocatable, dimension(:,:,:) ::  &
#else
      integer, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  &
#endif
       nmlb           ! layer containing mlb.
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension (2*max_nsteps_batrop+1,2) ::  &
#endif
          coeflx                ! used to diagnose the barotropic mass flux
!
! ---  s w i t c h e s    (if set to .true., then...)
! --- btrlfr      leapfrog barotropic time step
! --- btrmas      barotropic is mass conserving
! --- diagno      output model fields and diagnostic messages
! --- thermo      use thermodynamic forcing (flxflg>0)
! --- windf       use wind stress   forcing (wndflg>0)
! --- pcipf       use evap-precip surface salinity flux
! --- epmass      treat evap-precip as a mass exchange
! --- mslprf      use msl presssure forcing
! --- priver      use river precip bogas
! --- rivera      annual-only river precip bogas
! --- kparan      annual-only kpar or chl
! --- lbmont      use sshflg=2 correction on lateral baro. bndy
! --- relax       activate lateral boundary T/S/p  climatological nudging
! --- srelax      activate surface salinity        climatological nudging
! ---              (sssflg==1 or -1)
! --- trelax      activate surface temperature     climatological nudging
! ---              (sstflg==1)
! --- trcrlx      activate lateral boundary tracer climatological nudging
! --- relaxf      input T/S/p   relaxation fields
! --- relaxs      input surface relaxation fields only
! --- relaxt      input tracer  relaxation fields
! --- locsig      use locally-referenced potential density for stability
! --- vsigma      use spacially varying target densities
! --- hybrid      use hybrid vertical coordinates
! --- isopyc      use isopycnic vertical coordinates (MICOM mode)
! --- icegln      use energy loan ice model (iceflg==1)
! --- hybraf      HYBGEN: apply Robert-Asselin filter during hybgen
! --- isopcm      HYBGEN: use PCM to remap isopycnal layers
! --- mxl_no      NO:  deactivate             mixed layer model (mlflag==0)
! --- mxlkta      KT:    activate    original mixed layer model (mlflag==2)
! --- mxlktb      KT:    activate alternative mixed layer model (mlflag==3)
! --- mxlkrt      KT:    activate MICOM or HYCOM Kraus-Turner (mlflag==2,3)
! --- pensol      KT:    activate penetrating solar radiation
! --- mxlkpp      KPP:   activate mixed layer model (mlflag==1)
! --- bblkpp      KPP:   activate bottom boundary layer
! --- shinst      KPP:   activate shear instability mixing
! --- dbdiff      KPP:   activate double diffusion  mixing
! --- nonloc      KPP:   activate nonlocal b. layer mixing
! --- botdiw      KPROF: activate bot.enhan.int.wav mixing
! --- difout      KPROF: output visc/diff coeffs in archive
! --- mxlmy       MY2.5: activate mixed layer model (mlflag==5)
! --- mxlpwp      PWP:   activate mixed layer model (mlflag==4)
! --- mxlgiss     GISS:  activate mixed layer model (mlflag==6)
! --- stroff      add a net wind stress offset
! --- flxoff      add a net heat flux offset
! --- flxsmo      activate smoothing of surface fluxes
! --- trcrin      initialize tracer from restart file
! --- trcout      advect tracer and save results in history/restart file
! --- dsur1p      single point only surface diagnostics
! --- arcend      always write a 3-d archive at the end of the run
!
      logical, save :: &
                    btrlfr,btrmas,diagno,thermo,windf,mslprf, &
                    pcipf,epmass,priver,rivera,kparan,lbmont, &
                    relax,srelax,trelax,trcrlx,relaxf,relaxs,relaxt, &
                    locsig,vsigma,hybrid,isopyc,icegln,hybraf,isopcm, &
                    mxl_no,mxlkta,mxlktb,mxlkrt,pensol, &
                    mxlkpp,bblkpp,shinst,dbdiff,nonloc, &
                    botdiw,difout, &
                    mxlmy,mxlpwp,mxlgiss, &
                    stroff,flxoff,flxsmo,trcrin,trcout, &
                    dsur1p,arcend
!
! ---  t e x t
! ---  ctitle     four lines describing the simulation
!
      character*80, save, dimension(4) :: ctitle
!
! --- atmospheric forcing fields, natm is 2 or 4
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm) ::  &
#endif
       taux,           & ! wind stress in x direction
       tauy,           & ! wind stress in y direction
       wndspd,         & ! wind speed
       ustara,         & ! ustar (tke source)
       mslprs,         & ! mean sea level pressure (anomaly)
       airtmp,         & ! air temperature
       vapmix,         & ! atmosph. vapor mixing ratio
       precip,         & ! precipitation
       radflx,         & ! net solar radiation
       swflx,          & ! net shortwave radiation
       surtmp,         & ! surface temp. used to calculate input lw radiation
       seatmp,         & ! best available SST from observations
       stoc_t,         & ! stochastic temperature anomaly forcing
       stoc_s,         & ! stochastic salinty     anomaly forcing
       stoc_u,         & ! stochastic u-velocity  anomaly forcing
       stoc_v            ! stochastic v-velocity  anomaly forcing
!
! --- monthly atmospheric forcing fields
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::  &
#endif
       akpar,          & ! photosynthetically available radiation coefficent or chlorophyll-a (jerlov=-1)
       rivers            ! river inflow bogused to surface precipitation
 
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       oftaux,         & ! wind stress offset in x direction
       oftauy,         & ! wind stress offset in y direction
       offlux            ! net heat flux offset

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(0:24,-91:91) ::  &
#endif
       diurnl         ! hourly vs latitude shortwave scale factor table
!
! --- surface and sidewall and nestwall boundary fields
! ---  (kkwall and kknest are either kdm or, if inactive, 1).
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4) ::  &
#endif
       pwall,          & ! pressure b.c. at sidewalls
       swall,          & ! salinity b.c. at sidewalls
       twall             ! temp.    b.c. at sidewalls

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:,:) ::  &
#else
      real, save,  &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4, &
                                                         mxtrcr) ::  &
#endif
       trwall         ! tracer   b.c. at sidewalls

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2) ::  &
#endif
       pnest,          & ! pressure b.c. at nestwalls
       snest,          & ! salinity b.c. at nestwalls
       tnest,          & ! temp.    b.c. at nestwalls
       unest,          & ! u-vel.   b.c. at nestwalls
       vnest             ! v-vel.   b.c. at nestwalls

      real, save, dimension(2) ::  &
       un1min,          & ! u-vel.   b.c. at nestwalls, minimum layer 1
       un1max,          & ! u-vel.   b.c. at nestwalls, maximum layer 1
       vn1min,          & ! v-vel.   b.c. at nestwalls, minimum layer 1
       vn1max             ! v-vel.   b.c. at nestwalls, maximum layer 1
 
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::  &
#endif
       ubnest,         & ! barotropic u-velocity at nestwalls
       vbnest,         & ! barotropic v-velocity at nestwalls
       ubpnst,         & ! barotropic u-velocity at nestwalls on p-grid
       vbpnst,         & ! barotropic v-velocity at nestwalls on p-grid
       pbnest            ! barotropic pressure   at nestwalls

#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       rmu,            & ! weights for   s.w.b.c. relax
       rmunp,          & ! weights for p.n.w.b.c. relax
       rmunv,          & ! weights for v.n.w.b.c. relax
       rmunvu,         & ! weights for u.n.w.b.c. relax (u-grid, masked)
       rmunvv,         & ! weights for v.n.w.b.c. relax (v-grid, masked)
       rmutra,         & ! weights for tracr.b.c. relax (maximum of all tracers)
       rmus              ! weights for        sss relax !!Alex


#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,mxtrcr) ::  &
#endif
       rmutr          ! weights for tracr.b.c. relax

#if defined(RELO)
      integer, save, allocatable, dimension(:,:) ::  &
#else
      integer, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::  &
#endif
       maskbc         ! mask for nested barotropic boundary condition
!
! --- pwp variables
      real, save :: &
       rigc            & ! PWP: critical gradient richardson number
      ,ribc              ! PWP: critical bulk richardson number
!
! --- m-y 2.5 variables
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1,2) ::  &
#endif
       q2              & !  tke
      ,q2l               !  tke * turbulent length scale

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, &
            dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1) ::  &
#endif
       difqmy          & !  tke diffusivity
      ,vctymy          & !  viscosity on mellor-yamada vertical grid
      ,diftmy            !  temperature diffusivity on mellor-yamada vertical grid

      real, save :: &
       ghc             & !  constant for calculating tke production
      ,sef             & !  constant for calculating tke production
      ,smll            & !  constant for calculating tke
      ,const1          & !  coefficient for estimating surface and bottom bc's
      ,coef4           & !  coefficient for calculating  viscosity/diffusivity
      ,coef5           & !  coefficient for calculating  viscosity/diffusivity
      ,a1my            & !  coefficient for calculating  viscosity/diffusivity
      ,b1my            & !  coefficient for calculating  viscosity/diffusivity
      ,a2my            & !  coefficient for calculating  viscosity/diffusivity
      ,b2my            & !  coefficient for calculating  viscosity/diffusivity
      ,c1my            & !  coefficient for calculating  viscosity/diffusivity
      ,e1my            & !  coefficient for calculating  viscosity/diffusivity
      ,e2my            & !  coefficient for calculating  viscosity/diffusivity
      ,e3my              !  coefficient for calculating  viscosity/diffusivity
!
! --- kpp variables
#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) ::  &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::  &
#endif
       zgrid           & !  grid levels in meters
      ,vcty            & !  vert. viscosity coefficient
      ,difs            & !  vert. scalar diffusivity
      ,dift            & !  vert. temperature diffusivity
      ,ghats             !  nonlocal transport

      real, save :: &
       vonk            & !  von karman constant
      ,zmin,zmax       & !  zehat limits for table
      ,umin,umax       & !  ustar limits for table
      ,epsilon         & ! vertical coordinate scale factor
      ,cmonob          & ! KPP:    scale factor for Monin-Obukov length
      ,cekman          & ! KPP:    scale factor for Ekman depth
      ,qrinfy          & ! KPP:    1/max grad.rich.no. for shear instability
      ,difm0           & ! KPP:    max viscosity    due to shear instability
      ,difs0           & ! KPP:    max diffusivity  due to shear instability
      ,difmiw          & ! KPP/MY: background/internal wave viscosity   (m^2/s)
      ,difsiw          & ! KPP/MY: background/internal wave diffusivity (m^2/s)
      ,dsfmax          & ! KPP:    salt fingering diffusivity factor    (m^2/s)
      ,rrho0           & ! KPP:    salt fingering rp=(alpha*delT)/(beta*delS)
      ,ricr            & ! KPP:    critical bulk richardson number
      ,cs              & ! KPP:    value for nonlocal flux term
      ,cstar           & ! KPP:    value for nonlocal flux term
      ,cv              & ! KPP:    buoyancy frequency ratio (0.0 to use a fn. of N)
      ,c11             & ! KPP:    value for turb velocity scale
      ,deltaz          & ! delta zehat in table
      ,deltau          & ! delta ustar in table
      ,vtc             & ! constant for estimating background shear in rib calc.
      ,cg              & ! constant for estimating nonlocal flux term of diff. eq.
      ,dp0enh            ! dist. for tapering diff. enhancement at interface nbl-1

      integer, save :: &
       niter           & ! KPP: iterations for semi-implicit soln. (2 recomended)
      ,hblflg            ! KPP: b.layer interpolation flag (0=con.1=lin.,2=quad.)
!
! --- nasa giss variables
!
      integer, save :: &
              nextrtbl0,ifexpabstable,nextrtbl1, &
              nextrtbl,nposapprox,mt0,mt,ntbl, &
              mt_ra_r,n_theta_r_oct,nbig
!
      real, save :: &
              deltheta_r,pidbl,rri
!
      integer, save :: &
       irimax(-762:762) &
      ,nb
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::  &
#else
      real, save, dimension(-762:762,-nlgiss:nlgiss) ::  &
#endif
       slq2b &
      ,smb &
      ,shb &
      ,ssb
!
      real, save :: &
       ribtbl(-762:762) &
      ,ridb(  -762:762) &
      ,dri &
      ,back_ra_r(-39:117) &
      ,sisamax(  -39:117) &
      ,ra_rmax(  -39:117) &
      ,c_y_r0(   -39:117) &
      ,sm_r1(    -39:117) &
      ,sh_r1(    -39:117) &
      ,ss_r1(    -39:117) &
      ,slq2_r1(  -39:117) &
      ,b1,visc_cbu_limit,diff_cbt_limit &
      ,theta_rcrp,theta_rcrn
!
      integer, save :: &
              ifback,ifsali,ifepson2,ifrafgmax, &
              ifsalback,ifchengcon,ifunreal,idefmld, &
              ifpolartablewrite,ifbg_theta_interp
!
      real, save :: &
             back_ph_0,adjust_gargett,back_k_0,back_del_0,back_s2, &
             ri0,ebase,epson2_ref, &
             eps_bot0,scale_bot,       & !for bottom-enhanced
             eplatidepmin,wave_30,     & !and latitude dependent mixing
             deltemld,delrhmld, &
             back_sm2,v_back0,t_back0, &
             s_back0,ri_internal,backfrac,backfact,ako,tpvot0,sgmt, &
             tptot0,tpcot0,ttot0,tcot0,tctot0,tpvot,tptot,tpcot, &
             ttot,tcot,tctot,back_l_0
!
      real*8, save :: area,avgbot,watcum,empcum,wndrep
!
      real, save :: &
                     time,delt1,dlt, &
                      w0, w1, w2, w3,   & ! wind  interp. scale factors
                      wk0,wk1,wk2,wk3,  & ! kpar  interp. scale factors
                      wr0,wr1,wr2,wr3,  & ! river interp. scale factors
                      wc0,wc1,wc2,wc3,  & ! clim. interp. scale factors
                      wn0,wn1,          & ! nest  interp. scale factors
                      wb0,wb1             ! baro. interp. scale factors
!
      integer, save :: &
                      nstep,nstep1,nstep2,lstep, &
                      l0, l1, l2, l3,   & ! wind  indexes
                      lk0,lk1,lk2,lk3,  & ! kpar  indexes
                      lr0,lr1,lr2,lr3,  & ! river indexes
                      lc0,lc1,lc2,lc3,  & ! clim. indexes
                      ln0,ln1,          & ! nest  indexes
                      lb0,lb1             ! baro. indexes
!
! --- 'sigma ' = isopyncnal layer target densities (sigma units)
! --- 'thbase' = reference density (sigma units)
! --- 'saln0'  = initial salinity value
! --- 'baclin' = baroclinic time step
! --- 'batrop' = barotropic time step
! ---'qhybrlx' = HYBGEN: relaxation coefficient (inverse baroclinic time steps)
! --- 'hybiso' = HYBGEN: Use PCM if layer is within hybiso of target density
! --- 'visco2' = deformation-dependent Laplacian  viscosity factor
! --- 'visco4' = deformation-dependent biharmonic viscosity factor
! --- 'facdf4' =       speed-dependent biharmonic viscosity factor
! --- 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissipation
! --- 'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissipation
! --- 'temdf2' = diffusion velocity (m/s) for Laplacian  temp/saln diffusion
! --- 'temdfc' = temp diffusion conservation (0.0 all density, 1.0 all temp)
! --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
! --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion
! --- 'vertmx' = diffusion velocity (m/s) for mom.mixing across mix.layr.base
! --- 'tofset' = temperature anti-drift offset (degC/century)
! --- 'sofset' = salnity     anti-drift offset  (psu/century)
! --- 'diapyc' = KT: diapycnal diffusivity x buoyancy freq. (m**2/s**2)
! --- 'dtrate' = KT: maximum permitted m.l. detrainment rate (m/day)
! --- 'slip'   = +1 for free-slip, -1  for non-slip boundary conditions
! --- 'cb'     = coefficient of quadratic bottom friction
! --- 'cbar'   = rms flow speed (m/s) for linear bottom friction law
! --- 'drglim' = limiter for explicit friction (1.0 no limiter, 0.0 implicit)
! --- 'drgscl' = scale factor for tidal drag   (0.0 for no tidal drag)
! --- 'thkdrg' = thickness of bottom boundary layer for tidal drag (m)
! --- 'dsurfq' = number of days between model diagnostics at the surface
! --- 'diagfq' = number of days between model diagnostics
! --- 'proffq' = number of days between model diagnostics at some locations
! --- 'tilefq' = number of days between model diagnostics on some tiles
! --- 'meanfq' = number of days between model diagnostics (time averaged)
! --- 'rstrfq' = number of days between model restart output
! --- 'bnstfq' = number of days between baro. nesting archive input
! --- 'nestfq' = number of days between 3-d   nesting archive input
! --- 'stfflg' = stochastic anomaly forcing flag (0=no, 1=TS, 2=TSV)
! --- 'stfrdt' = stochastic T anomaly forcing e-folding depth (m)
! --- 'stfrds' = stochastic S anomaly forcing e-folding depth (m)
! --- 'stfrdv' = stochastic V anomaly forcing e-folding depth (m)
! --- 'ra2fac' = weight for Robert-Asselin time filter
! --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
! --- 'sefold' = e-folding time                  for SSS relaxation (days)
! --- 'thkmls' = reference mixed-layer thickness for SSS relaxation (m)
! --- 'thkmlt' = reference mixed-layer thickness for SST relaxation (m)
! --- 'thkriv' = nominal thickness of river inflow (m)
! --- 'thkcdw' = thickness for near-surface currents in ice-ocean stress
! --- 'thkfrz' = maximum thickness of near surface freezing zone (m)
! --- 'tfrz_0' = ENLN: ice melting point (degC) at S=0psu
! --- 'tfrz_s' = ENLN: gradient of ice melting point (degC/psu)
! --- 'ticegr' = ENLN: vertical temperature gradient inside ice (deg/m)
! ---                   (0.0 to get ice surface temp. from atmos. surtmp)
! --- 'hicemn' = ENLN: minimum ice thickness (m)
! --- 'hicemx' = ENLN: maximum ice thickness (m)
! --- 'thkmin' = KT/PWP: minimum mixed-layer thickness (m)
! --- 'bldmin' = KPP:    minimum surface boundary layer thickness (m)
! --- 'bldmax' = KPP:    maximum surface boundary layer thickness (m)
! --- 'thkbot' = thickness of bottom boundary layer (m)
! --- 'sigjmp' = minimum density jump across interfaces   (theta units)
! --- 'tmljmp' = equivalent temperature jump across the mixed layer (degC)
! --- 'prsbas' = msl pressure is input field + prsbas (Pa)
! --- 'salmin' = minimum salinity allowed in an isopycnic layer (psu)
! --- 'dp00'   = deep    z-level spacing minimum thickness (m)
! --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
! --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
! --- 'ds00'   = shallow z-level spacing minimum thickness (m)
! --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
! --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
! --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
! --- 'isotop' = shallowest depth for isopycnal layers     (m), <0 from file
! --- 'oneta0' = minimum 1+eta, must be > 0.0
! --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
! --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
! --- 'hybmap' = HYBGEN:  remapper  flag (0=PCM,1=PLM,2=PPM,-ve:isoPCM)
! --- 'hybflg' = HYBGEN:  generator flag (0=T&S,1=th&S,2=th&T)
! --- 'advflg' = thermal  advection flag (0=T&S,1=th&S,2=th&T)
! --- 'advtyp' = scalar   advection type (0=PCM,1=MPDATA,2=FCT2,4=FCT4)
! --- 'momtyp' = momentum advection type (2=2nd order, 4=4th order)
! --- 'kapref' = thermobaric reference state (-1=input,0=none,1,2,3=constant)
! --- 'kapnum' = number of thermobaric reference states (1 or 2)
! --- 'tsofrq' = number of time steps between anti-drift offset calcs
! --- 'mixfrq' = KT: number of time steps between diapycnal mixing calcs
! --- 'icefrq' = relax to tfrz with e-folding time of icefrq time steps
! --- 'icpfrq' = number of time steps between sea-ice updates
! --- 'ntracr' = number of tracers (<=mxtrcr)
! --- 'trcflg' = tracer type flag (one per tracer)
! --- 'clmflg' = climatology frequency flag (6=bimonthly,12=monthly)
! --- 'dypflg' = KT: diapycnal mixing flag (0=none,1=KPP,2=explicit)
! --- 'iniflg' = initial state flag (0=level,1=zonal,2=climatology)
! --- 'lbflag' = lateral baro. bndy flag (0=none;nest:2=B-K,4=Flather,6=clamped)
! ---             (port: 1=Browning-Kreiss,3=Flather)
! --- 'mapflg' = map flag (0=mercator,2=uniform,3=beta-plane,4=input)
! --- 'yrflag' = days in year flag   (0=360,1=366,2=366Jan1,3=actual)
! --- 'sshflg' = diagnostic SSH flag (0=SSH,1=SSH&stericSSH,2=SSH&stericMONTG)
! --- 'iversn' = hycom version number x10
! --- 'iexpt'  = experiment number x10
! --- 'jerlv0' = initial jerlov water type (1 to 5; 0 for kpar, -1 for chl)
! --- 'iceflg' = sea ice model flag   (0=none,1=energy loan,2=coupled/esmf)
! --- 'ishelf' = ice shelf flag       (0=none,1=ice shelf over ocean)
! --- 'wndflg' = wind stress input flag (0=none,1=uv-grid,2,3=p-grid,4,5=wnd10m)
! --- 'amoflg' = relative wind     flag (0=U10:wndflg=-4,-5;1=U10-Uocn:wndflg=4,5)
! --- 'ustflg' = ustar forcing flag          (3=input,1=wndspd,2=stress)
! --- 'flxflg' = thermal forcing flag (0=none,3=net-flux,1,2,4=sst-based)
! --- 'empflg' = E-P     forcing flag (0=none,3=net_E-P, 1,2,4=sst-based_E)
! --- 'empbal' = E-P     balance flag (0=none,1=offset,2=scale)
! --- 'dswflg' = diurnal shortwv flag (0=none,1=daily to diurnal correction)
! --- 'albflg' = ocean albedo    flag (0=none,1=.06,2=L&Y)
! --- 'sssflg' = SSS relaxation  flag (0=none,1=clim,-1=clim+rmx)
! --- 'sssbal' = SSS rlx balance flag (0=none,1=offset,2=scale)
! --- 'lwflag' = longwave corr.  flag (0=none,1=clim,2=atmos), sst-based
! --- 'sstflg' = SST relaxation  flag (0=none,1=clim,2=atmos,3=obs)
! --- 'icmflg' = ice mask        flag (0=none,1=clim,2=atmos,3=obs)
! --- 'difsmo' = KPROF: number of layers with horiz smooth diff coeffs
!
#if defined(RELO)
      real, save, allocatable, dimension(:) ::  &
#else
      real, save, dimension(kdm) ::  &
#endif
        sigma, &
        salmin
!
      real, save :: &
                     thbase,saln0,baclin,batrop, &
                     qhybrlx,hybiso, &
                     visco2,visco4,veldf2,veldf4,facdf4, &
                     temdf2,temdfc,thkdf2,thkdf4,vertmx,diapyc, &
                     tofset,sofset,dtrate,slip,cb,cbar, &
                     drglim,drgscl,thkdrg, &
                     dsurfq,diagfq,proffq,tilefq,meanfq, &
                     rstrfq,bnstfq,nestfq, &
                     stfrdt,stfrds,stfrdv,ra2fac,wbaro, &
                     sefold, &
                     thkmls,thkmlt,thkriv,thkmin,bldmin,bldmax,thkbot, &
                     thkcdw,thkfrz,tfrz_0,tfrz_s,ticegr,hicemn,hicemx, &
                     dp00,dp00f,dp00x,ds00,ds00f,ds00x,dp00i,isotop, &
                     oneta0,sigjmp,tmljmp,prsbas,emptgt
!
      integer, save :: &
                     tsofrq,mixfrq,icefrq,icpfrq,nhybrd,nsigma, &
                     hybmap,hybflg,advflg,advtyp,momtyp,stfflg, &
                     kapref,kapnum, &
                     ntracr,trcflg(mxtrcr), &
                     clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,sshflg, &
                     iversn,iexpt,jerlv0, &
                     iceflg,ishelf,icmflg,wndflg,amoflg,ustflg, &
                     flxflg,empflg,dswflg,albflg,lwflag,sstflg,sssflg, &
                     empbal,sssbal, &
                     difsmo,disp_count
!
      real,    parameter :: &
        g      =   9.806,       & !gravitational acceleration (m/s**2)
        rhoref =   1000.0,      & !reference value of potential density (kg/m**3)  
        svref  =   1.0/rhoref,  & !reference value of specific volume (m**3/kg)
        onem   = g/svref,       & !one meter in pressure units (Pa) 
                                  !(/svref instead of *rhoref for compatibility with older version)
        tenm   = onem*10.0, &
        tencm  = onem* 0.1, &
        onecm  = onem* 0.01, &
        onemm  = onem* 0.001, &
        spcifh = 3990.0,       & !specific heat of sea water (j/kg/deg)
        epsil  = 1.0d-11,      & !small nonzero to prevent division by zero
        hugel  = 2.0**100,     & !large number used to indicate land points
        pi     = 3.14159265358979323846d0, &
        radian = pi/180.0, &
        qonem  = 1.0/onem 
!
! --- grid point where detailed diagnostics are desired:
      integer, save :: itest,jtest,ittest,jttest
!
! --- filenames.
      character*80, save :: &
                    flnmdep,flnmgrd,flnmshlf, &
                    flnmrsi,flnmrso, flnmflx, &
                    flnmarc,flnmovr,flnmfor,flnmforw,flnminp, &
                    flnmarcm, &
                    flnmarcs, &
                    flnmarct
!
! --- CCSM3 variables (should be in CCSM3 modules)
      logical, save :: dosstbud,doovtn,chk_ovtn,dodump,dohist,dorestart
!
      integer, save :: istrt_mo,icurrent_mo,istrt_yr,icurrent_yr
!
! --- diurnal cycle factor for short wave heat flux
      integer, save :: nsteps_per_day,nsteps_today

#if ! defined (ESPC_COUPLE)
! --- needed for restart if coupled
#  if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
#  else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
#  endif
       dhde,dhdn,             & ! sea-surface height slope for CICE
       umxl,vmxl,             & ! surface u and v for CICE
       tml,sml                  ! surface T and S for CICE
#endif

#if defined (USE_NUOPC_GENERIC)
! ---  import from atm
      logical cpl_taux, cpl_tauy, cpl_wndspd, cpl_ustara, &
       cpl_airtmp, cpl_vapmix, cpl_precip, cpl_surtmp, cpl_seatmp

! ---  import from ice
      logical cpl_sic, cpl_sitx, cpl_sity, cpl_siqs, cpl_sifh, &
              cpl_sifs, cpl_sifw, cpl_sit, cpl_sih, cpl_siu, &
              cpl_siv

!
#  if defined(RELO)
      real, target, allocatable,dimension (:,:) :: &
             sic_import, & !Sea Ice Concentration
            sitx_import, & !Sea Ice X-Stress
            sity_import, & !Sea Ice Y-Stress
            siqs_import, & !Solar Heat Flux thru Ice to Ocean
            sifh_import, & !Ice Freezing/Melting Heat Flux
            sifs_import, & !Ice Freezing/Melting Salt Flux
            sifw_import, & !Ice Net Water Flux
             sit_import, & !Sea Ice Temperature
             sih_import, & !Sea Ice Thickness
             siu_import, & !Sea Ice X-Velocity
             siv_import    !Sea Ice Y-Velocity

#  else
      real, target, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
             sic_import, & !Sea Ice Concentration
            sitx_import, & !Sea Ice X-Stress
            sity_import, & !Sea Ice Y-Stress
            siqs_import, & !Solar Heat Flux thru Ice to Ocean
            sifh_import, & !Ice Freezing/Melting Heat Flux
            sifs_import, & !Ice Freezing/Melting Salt Flux
            sifw_import, & !Ice Net Water Flux
             sit_import, & !Sea Ice Temperature
             sih_import, & !Sea Ice Thickness
             siu_import, & !Sea Ice X-Velocity
             siv_import    !Sea Ice Y-Velocity
#  endif

#endif /* USE_NUOPC_GENERIC */

#if defined (USE_NUOPC_CESMBETA)
! --- Average array for export
#  if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
#  else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
#  endif
       sshm,               & ! sea-surface height averaged over 1 coupling sequence
       um,vm,              & ! surface u and v
       ubm,vbm,            & ! surface ubaro and vbaro
       tavgm,savgm,        & ! surface T and S
       frzh                  ! Freezing potential  flux    (W/m2)

!
! --- NUOPC glue code structures
! --- tripolar grid
      logical ltripolar
! --- import from coupler
      real ocn_cpl_frq
! --- precipitation factor for coupled simulation
      real pcp_fact  ! always 1. : no precipiation adjustment
! ---  import from atm
      real nstep1_cpl,nstep2_cpl
      logical cpl_swflx, cpl_lwmdnflx, cpl_lwmupflx, &
       cpl_latflx, cpl_sensflx, &
       cpl_orivers,cpl_irivers

      real cpl_w2, cpl_w3
      logical cpl_implicit

#  if defined(RELO)
      real, target, allocatable,dimension (:,:,:) :: &
       imp_taux, imp_tauy, imp_taue, imp_taun, imp_wndspd, imp_ustara, &
       imp_airtmp, imp_vapmix, imp_swflx, imp_lwdflx, imp_lwuflx, &
       imp_latflx, imp_sensflx, &
       imp_precip, imp_surtmp, imp_seatmp, &
       imp_orivers,imp_irivers

#  else
      real, target, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) :: &
       imp_taux, imp_tauy, imp_taue, imp_taun, imp_wndspd, imp_ustara, &
       imp_airtmp, imp_vapmix, imp_swflx, imp_lwdflx, imp_lwuflx, &
       imp_latflx, imp_sensflx, &
       imp_precip, imp_surtmp, imp_seatmp, &
       imp_orivers,imp_irivers
#  endif
#endif /* USE_NUOPC_CESMBETA */


!
#if defined (ESPC_COUPLE)
! ---  import from atm
      logical cpl_swflx_net, cpl_lwflx_net, &
       cpl_swflx_net2down,cpl_lwflx_net2down, &
       cpl_swflxd,cpl_lwflxd, &
       cpl_sbhflx,cpl_lthflx, &
       cpl_u10,cpl_v10,cpl_mslprs


#  if defined(RELO)
      real, target, allocatable,dimension (:,:) :: &
        exp_sbhflx, &
        exp_lthflx

#  else
      real, target, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
        exp_sbhflx, &
        exp_lthflx
#  endif
#endif /* ESPC_COUPLE */

      contains

      subroutine cb_allocate
!
! --- Allocate saved arrays
!
      call set_r_init
!
#if defined(RELO)
      call gindex_allocate   !from mod_dimensions
!
      allocate( &
                    u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                    v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                   dp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                  dpo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                  dpu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                  dpv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                 temp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                 saln(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                 th3d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
               thstar(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2), &
                montg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2) )
      call mem_stat_add( 11*(idm+2*nbdy)*(jdm+2*nbdy)*kdm*2 )
#endif
                    u = r_init
                    v = r_init
                   dp = r_init
                  dpo = r_init
                  dpu = r_init
                  dpv = r_init
                 temp = r_init
                 saln = r_init
                 th3d = r_init
               thstar = r_init
                montg = r_init
!
      if     (ntracr.gt.0) then
#if defined(RELO)
        allocate( &
              tracer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,ntracr) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*kdm*2*ntracr )
#endif
              tracer = r_init
      endif
!
#if defined(RELO)
      allocate( &
                otemp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                osaln(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                oth3d(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) )
      call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
#endif
                otemp = r_init
                osaln = r_init
                oth3d = r_init
!
      if (mxlmy) then
#if defined(RELO)
        allocate( &
                    oq2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1), &
                   oq2l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1) )
        call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*(kkmy25+2) )
#endif
                    oq2 = r_init
                   oq2l = r_init
      endif !mxlmy
!
      if     (ntracr.gt.0) then
#if defined(RELO)
        allocate( &
             otracer(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,ntracr) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*kdm*ntracr )
#endif
             otracer = r_init
      endif
!
#if defined(RELO)
      allocate( &
              dpmold(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 1*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
              dpmold = r_init
!
#if defined(RELO)
      allocate( &
                    p(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                   pu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                   pv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) )
      call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy)*(kdm+1) )
#endif
                    p = r_init
                   pu = r_init
                   pv = r_init
!
#if defined(RELO)
      allocate( &
               theta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
              diaflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) )
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
#endif
               theta = r_init
              diaflx = r_init
!
#if defined(RELO)
      allocate( &
              corio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             srfhgt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             steric(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             sshgmn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             thmean(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
            montg_c(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             montg1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               skap(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 8*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
              corio = r_init
             srfhgt = r_init
             steric = r_init
             sshgmn = r_init
             thmean = r_init
            montg_c = 0.0 
             montg1 = r_init
               skap = r_init
!
#if defined(RELO)
      allocate( &
               psikk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                thkk(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
#endif
               psikk = r_init
                thkk = r_init
!
#if defined(RELO)
      allocate( &
               kapi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)/2 )  !real=2*int
#endif
               kapi = -99
!
#if defined(RELO)
      allocate( &
                uflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                vflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
              uflxav(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
              vflxav(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                dpav(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
#endif
                uflx = r_init
                vflx = r_init
              uflxav = r_init
              vflxav = r_init
                dpav = r_init
!
#if defined(RELO)
      allocate( &
               ubavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3), &
               vbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3), &
               pbavg(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) )
      call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy)*3 )
#endif
               ubavg = r_init
               vbavg = r_init
               pbavg = r_init
!
#if defined(RELO)
      allocate( &
               oneta(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
              onetao(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
            onetamas(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
#endif
               oneta = 1.0
              onetao = 1.0
            onetamas = 1.0
!
#if defined(RELO)
      allocate( &
               ubrhs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               vbrhs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               utotm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               vtotm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               utotn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               vtotn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) , &
               uflux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               vflux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              uflux2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              vflux2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              uflux3(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              vflux3(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
            onetacnt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 13*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
               ubrhs = r_init
               vbrhs = r_init
               utotm = r_init
               vtotm = r_init
               utotn = r_init
               vtotn = r_init
               uflux = r_init
               vflux = r_init
              uflux2 = r_init
              vflux2 = r_init
              uflux3 = r_init
              vflux3 = r_init
            onetacnt = 1.0
!
#if defined(RELO)
      allocate( &
               util1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               util2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               util3(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               util4(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               util5(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               util6(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                plon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                plat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                qlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                qlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                ulon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                ulat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vlon(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vlat(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                pang(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scuy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scvx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scpx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scpy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scqx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scqy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scu2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scv2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scp2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                scq2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               scp2i(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               scq2i(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               scuxi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               scvyi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               aspux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               aspuy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               aspvx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               aspvy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             veldf2u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             veldf2v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             veldf4u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             veldf4v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             thkdf4u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             thkdf4v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 cbp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               cbarp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                pgfx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                pgfy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               gradx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               grady(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              depthu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              depthv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              pvtrop(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              depths(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                drag(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              salfac(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              topiso(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                diws(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                diwm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              diwbot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              diwqh0(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              sssrmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           tidepg_mn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           displd_mn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           dispqd_mn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 62*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
               util1 = r_init
               util2 = r_init
               util3 = r_init
               util4 = r_init
               util5 = r_init
               util6 = r_init
                plon = r_init
                plat = r_init
                qlon = r_init
                qlat = r_init
                ulon = r_init
                ulat = r_init
                vlon = r_init
                vlat = r_init
                pang = r_init
                scux = r_init
                scuy = r_init
                scvx = r_init
                scvy = r_init
                scpx = r_init
                scpy = r_init
                scqx = r_init
                scqy = r_init
                scu2 = r_init
                scv2 = r_init
                scp2 = r_init
                scq2 = r_init
               scp2i = r_init
               scq2i = r_init
               scuxi = r_init
               scvyi = r_init
               aspux = r_init
               aspuy = r_init
               aspvx = r_init
               aspvy = r_init
             veldf2u = r_init
             veldf2v = r_init
             veldf4u = r_init
             veldf4v = r_init
             thkdf4u = r_init
             thkdf4v = r_init
                 cbp = r_init
               cbarp = r_init
                pgfx = r_init
                pgfy = r_init
               gradx = r_init
               grady = r_init
              depthu = r_init
              depthv = r_init
              pvtrop = r_init
              depths = r_init
                drag = r_init
              salfac = r_init
              topiso = r_init
                diws = r_init
                diwm = r_init
              diwbot = r_init
              diwqh0 = r_init
              sssrmx = r_init
           tidepg_mn = r_init
           displd_mn = r_init
           dispqd_mn = r_init
!
#if defined(RELO)
      allocate( &
              drgten(1:2,1:2,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
#endif
              drgten = r_init
!
#if defined(RELO)
      allocate( &
                 uja(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 ujb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 via(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 vib(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                pbot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sgain(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               surtx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               surty(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              surflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              sswflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              mixflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              sstflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              sssflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              rivflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              salflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              wtrflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              buoflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              bhtflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               ustar(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              ustarb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              turgen(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              thkice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              covice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              temice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              flxice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              fswice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              wflice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              wflfrz(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              sflice(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                si_c(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                si_h(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                si_t(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                si_u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                si_v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               si_tx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               si_ty(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 36*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
                 uja = r_init
                 ujb = r_init
                 via = r_init
                 vib = r_init
                pbot = r_init
               sgain = r_init
               surtx = r_init 
               surty = r_init
              surflx = r_init
              sswflx = r_init
              mixflx = r_init
              sstflx = r_init
              sssflx = r_init
              rivflx = r_init
              salflx = r_init
              wtrflx = 0.0 
              buoflx = r_init
              bhtflx = r_init
               ustar = r_init
              ustarb = r_init
              turgen = r_init
              thkice = r_init
              covice = r_init
              temice = r_init
              flxice = r_init
              fswice = r_init
              wflice = r_init
              wflfrz = 0.0     !diagnostic, icloan only
              sflice = r_init
                si_c = r_init
                si_h = r_init
                si_t = r_init
                si_u = r_init
                si_v = r_init
               si_tx = r_init
               si_ty = r_init
!
      if     (flxflg.eq.6) then
#if defined(RELO)
        allocate(wndocn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
#endif
                 wndocn = r_init
      endif !flxflg.eq.6
!
#if defined(RELO)
      allocate( &
              klist(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             jerlov(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy) ) !real=2*int
#endif
              klist = -99
             jerlov = -99
!
#if defined(RELO)
      allocate( &
              dpmixl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
               t1sav(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
               s1sav(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                tmlb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                smlb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
#endif
              dpmixl = r_init
               t1sav = r_init
               s1sav = r_init
                tmlb = r_init
                smlb = r_init
!
#if defined(RELO)
      allocate( &
              hekman(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              hmonob(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                dpbl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               dpbbl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                tmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                smix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               thmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                umix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 9*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
              hekman = r_init
              hmonob = r_init
                dpbl = r_init
               dpbbl = r_init
                tmix = r_init
                smix = r_init
               thmix = r_init
                umix = r_init
                vmix = r_init
!
#if defined(RELO)
      allocate( &
              nmlb(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy) ) !real=2*int
#endif
              nmlb = -99
!
#if defined(RELO)
      allocate( &
              coeflx(2*lstep+1,2) )
      call mem_stat_add( 2*(2*lstep+1) )
#endif
              coeflx = 0.0
!
#if defined(RELO)
      allocate( &
                 taux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
                 tauy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               wndspd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               ustara(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               mslprs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               airtmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               vapmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               radflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
                swflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               surtmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               seatmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               stoc_t(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               stoc_s(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               stoc_u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm), &
               stoc_v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,natm) )
      call mem_stat_add( 16*(idm+2*nbdy)*(jdm+2*nbdy)*natm )
      allocate( &
                akpar(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4), &
               rivers(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) )
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*4 )
#endif
                 taux = 0.0 ! ESPC r_init 
                 tauy = 0.0 ! ESPC r_init
               wndspd = 0.0
               ustara = 0.0
               mslprs = 0.0 ! ESPC r_init
               airtmp = 0.0
               vapmix = 0.0
               precip = 0.0
               radflx = 0.0
                swflx = 0.0
               surtmp = 0.0
               seatmp = 0.0
               stoc_t = 0.0
               stoc_s = 0.0
               stoc_u = 0.0
               stoc_v = 0.0
                akpar = 0.0 
               rivers = 0.0
!
#if defined(RELO)
      allocate( &
           oftaux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           oftauy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           offlux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
           offlux =  0.0  
           oftaux =  0.0 ! ESPC r_init
           oftauy =  0.0 ! ESPC r_init
!
#if defined(RELO)
      allocate( &
               swall(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4), &
               twall(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4) )
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*kkwall*4 )
#endif
               swall = r_init
               twall = r_init
!
#if defined(RELO)
      allocate( &
               diurnl(0:24,-91:91) )
      call mem_stat_add( 25*193 )
#endif
               diurnl = r_init
!
      if     (relaxt .or. kkwall.eq.kdm) then
#if defined(RELO)
        allocate( &
               pwall(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,4) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*kdm*4 )
#endif
               pwall = r_init
      endif
!
      if     (ntracr.gt.0) then
#if defined(RELO)
        allocate( &
              trwall(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4,ntracr) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*kkwall*4*ntracr )
#endif
              trwall = r_init
      endif
!
#if defined(RELO)
      allocate( &
               pnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2), &
               snest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2), &
               tnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2), &
               unest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2), &
               vnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2) )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*kknest*2 )
#endif
               pnest = r_init
               snest = r_init
               tnest = r_init
               unest = r_init
               vnest = r_init
!
#if defined(RELO)
      allocate( &
              ubnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
              vbnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
              ubpnst(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
              vbpnst(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
              pbnest(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*2 )
#endif
              ubnest = r_init
              vbnest = r_init
              ubpnst = r_init
              vbpnst = r_init
              pbnest = r_init
!
#if defined(RELO)
      allocate( &
                rmu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              rmunp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              rmunv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rmunvu(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rmunvv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
             rmutra(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               rmus(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 7*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif
                rmu = r_init
              rmunp = r_init
              rmunv = r_init
             rmunvv = r_init
             rmutra = r_init
             rmutra = r_init
               rmus = r_init  !!Alex
!
      if     (ntracr.gt.0) then
#if defined(RELO)
        allocate( &
               rmutr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ntracr) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*ntracr )
#endif
               rmutr = r_init
      endif
!
#if defined(RELO)
      allocate( &
           maskbc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)/2 ) !real=2*int
#endif
           maskbc = -99
!
      if (mxlmy) then
#if defined(RELO)
        kkmy25 = kk
        allocate( &
                     q2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1,2), &
                    q2l(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1,2) )
        call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*(kkmy25+2)*2 )
!
        allocate( &
                 difqmy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1), &
                 vctymy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1), &
                 diftmy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1) )
        call mem_stat_add( 3*(idm+2*nbdy)*(jdm+2*nbdy)*(kkmy25+2) )
#endif
                     q2 = r_init
                    q2l = r_init
                 difqmy = r_init
                 vctymy = r_init
                 diftmy = r_init
#if defined(RELO)
      else
        kkmy25 = -1
#endif
      endif !mxlmy
!
#if defined(RELO)
      allocate( &
                zgrid(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                 vcty(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                 difs(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                 dift(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1), &
                ghats(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*(kdm+1) )
#endif
                zgrid = r_init
                 vcty = r_init
                 difs = r_init
                 dift = r_init
                ghats = r_init
!
      if     (mxlgiss) then
#if defined(RELO)
        allocate( &
          slq2b(-762:762,-762:762), &
            smb(-762:762,-762:762), &
            shb(-762:762,-762:762), &
            ssb(-762:762,-762:762) )
        call mem_stat_add( 4*1525*1525 )
#endif
          slq2b = r_init
            smb = r_init
            shb = r_init
            ssb = r_init
      endif !mxlgiss
!

#if defined (USE_NUOPC_GENERIC)
#  if defined(RELO)
      allocate( &
                sic_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sitx_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sity_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               siqs_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sifh_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sifs_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               sifw_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                sit_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                sih_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                siu_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                siv_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 11*(idm+2*nbdy)*(jdm+2*nbdy) )
#  endif

               sic_import = 0.d0
              sitx_import = 0.d0
              sity_import = 0.d0
              siqs_import = 0.d0
              sifh_import = 0.d0
              sifs_import = 0.d0
              sifw_import = 0.d0
               sit_import = 0.d0
               sih_import = 0.d0
               siu_import = 0.d0
               siv_import = 0.d0

#endif /* USE_NUOPC_GENERIC */

#if ! defined (ESPC_COUPLE)
#if defined(RELO)
      allocate( &
                   dhde(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   dhdn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   umxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   vmxl(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                    tml(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                    sml(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)) 
      call mem_stat_add( 6*(idm+2*nbdy)*(jdm+2*nbdy) )
#endif

                dhde = 0.d0
                dhdn = 0.d0
                umxl = 0.d0
                vmxl = 0.d0
                 tml = 0.d0
                 sml = 0.d0
#endif

#if defined (USE_NUOPC_CESMBETA)
#  if defined(RELO)
      allocate( &
                  sshm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                    um(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                    vm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   ubm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   vbm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 tavgm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 savgm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                  frzh(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
      call mem_stat_add( 8*(idm+2*nbdy)*(jdm+2*nbdy) )
#  endif

                sshm = 0.d0
                  um = 0.d0
                  vm = 0.d0
                 ubm = 0.d0
                 vbm = 0.d0
               tavgm = 0.d0
               savgm = 0.d0
                frzh = 0.d0

#  if defined(RELO)
      allocate( &
                   imp_taux(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                   imp_tauy(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                   imp_taue(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                   imp_taun(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_wndspd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_ustara(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_airtmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_vapmix(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                  imp_swflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_lwdflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_lwuflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_latflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                imp_sensflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_precip(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                imp_irivers(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                imp_orivers(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_surtmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2), &
                 imp_seatmp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) )
      call mem_stat_add( 16*(idm+2*nbdy)*(jdm+2*nbdy)*(2) )
#  endif

                  imp_taux = 0.d0
                  imp_tauy = 0.d0
                  imp_taue = 0.d0
                  imp_taun = 0.d0
                imp_wndspd = 0.d0
                imp_ustara = 0.d0
                imp_airtmp = 0.d0
                imp_vapmix = 0.d0
                 imp_swflx = 0.d0
                imp_lwdflx = 0.d0
                imp_lwuflx = 0.d0
                imp_latflx = 0.d0
               imp_sensflx = 0.d0
                imp_precip = 0.d0
                imp_surtmp = 0.d0
                imp_seatmp = 0.d0
               imp_irivers = 0.d0
               imp_orivers = 0.d0
#endif /* USE_NUOPC_CESMBETA */

#if defined (ESPC_COUPLE)
      allocate( &
                   exp_sbhflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   exp_lthflx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
      call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
       exp_sbhflx(:,:) = 0.0 
       exp_lthflx(:,:) = 0.0 

#  endif

      end subroutine cb_allocate

      end module mod_cb_arrays
!
!> Revision history:
!>
!> Feb. 2001 - added halo and converted to f90 declarations
!> Aug. 2001 - added bnstfq,nestfq
!> Jan. 2002 - added curvilinear grid arrays and deleted /pivot/
!> May  2002 - added d[ps]00*
!> May  2002 - added PWP and MY2.5 mixed layer models
!> Aug. 2002 - added ntracr and trcflg
!> Nov. 2002 - added thkmls and thkmlt
!> Apr. 2003 - added dp00i, vsigma, and priver
!> May  2003 - added bldmin, bldmax, flxsmo, and icmflg
!> June 2003 - added locsig
!> Nov. 2003 - added advtyp
!> Jan. 2004 - added latdiw
!> Jan. 2004 - added bblkpp
!> Jan. 2004 - added hblflg
!> Feb. 2004 - added botdiw
!> Feb. 2004 - added temdfc
!> Mar. 2004 - added thkriv and epmass
!> Mar. 2004 - added isotop
!> Mar. 2005 - added tfrz_0, tfrz_s, ticegr, hicemn, and hicemx
!> Mar. 2005 - added tsofrq, tofset, and sofset
!> Mar. 2005 - added empflg
!> May  2005 - added kapref and kapnum
!> June 2006 - added surtx,surty
!> June 2006 - added icefrq,txice,tyice,uice,vice,flxice,fswice,sflice
!> June 2006 - added thkfrz
!> Jan. 2007 - added si_t; renamed si_[chuv] and si_t[xy]
!> Feb. 2007 - added CCSM3-only variables
!> Apr. 2007 - added dragrh,drglim,drgscl
!> Apr. 2007 - added srfhgt,montg1
!> Apr. 2007 - added btrlfr and btrmas
!> June 2007 - added momtyp and facdf4.
!> Sep. 2007 - added hybmap, hybiso and isopcm.
!> Feb. 2008 - added thkdrg.
!> Feb. 2008 - added sshflg and steric,sshgmn,thmean.
!> June 2008 - added tilefq.
!> Mar. 2008 - added dswflg and diurnl.
!> Dec. 2008 - difsmo is now an integer number of layers.
!> Jan. 2009 - added arcend.
!> June 2009 - added sssrmx.
!> Mar. 2010 - added diwbot.
!> Apr. 2010 - changed sssrmx to an array.
!> Apr. 2010 - added diwqh0; removed diwlat.
!> Apr. 2010 - added proffq.
!> Nov. 2010 - added wndrep.
!> Apr. 2011 - added flnmarcs.
!> Apr. 2011 - added cbarp.
!> Apr. 2011 - replaced huge with hugel, to avoid clash with intrinsic huge.
!> July 2011 - added salfac.
!> Aug. 2011 - replaced dpold and dpoldm with dpo
!> Aug. 2011 - added ra2fac, removed wts[12] and wuv[12]
!> Aug. 2011 - added hybraf
!> Sep. 2011 - added cbp.
!> Nov. 2011 - added icpfrq
!> Jan. 2012 - added thkcdw
!> Mar. 2012 - replaced dssk with dpns and dsns
!> Nov. 2012 - added stroff and oftaux,oftauy
!> Jan. 2013 - replaced dragrh with drgten
!> Apr. 2013 - added displd_mn and dispqd_mn and tidepg_mn
!> July 2013 - added diws and diwm
!> Oct. 2013 - jerlov constants now in swfrac_ij (thermf.f)
!> Nov. 2013 - added rivflx
!> Nov. 2013 - added albflg
!> Jan. 2014 - added prsbas, mslprf and mslprs
!> Jan. 2014 - /consts/ replaced by parameters
!> Jan. 2014 - /dicycr/ deleted
!> Jan. 2014 - added natm
!> Jan. 2014 - replaced common_blocks.h with mod_cb_arrays.F
!> Jan. 2014 - macro /* RELO */ for relocatable (region independent) version
!> Apr. 2014 - added ishelf and flnmshlf
!> Oct. 2014 - added wndocn
!> Feb. 2015 - added pang for coupled cases
!> Aug. 2015 - added mxl_no
!> Sep. 2017 - added stoc_[tsv]
!> Sep. 2017 - added stfflg and stfr[tsv]
!> Aug. 2018 - added sflfrz (now wflfrz)
!> Aug. 2018 - added "o" arrays for the asselin filter
!> Aug. 2018 - added oneta arrays
!> Nov. 2018 - virtual salt flux fields replaced with water and actual salt flux
!> Nov. 2018 - added sefold and rmus
!> Nov. 2018 - added epmbal and sssbal
!> Nov. 2018 - added emptgt
!> Dec. 2018 - added /* USE_NUOPC_CESMBETA */ macro for coupled simulation
!> Dec. 2018 - added /* USE_NUOPC_GENERIC */ and /* ESPC_COUPLE */ macros
!> Feb. 2019 - added montg_c
!> Feb. 2019 - removed onetai 
!> Sep. 2019 - five arrays moved to momtum_init
!> Sep. 2019 - added oneta0
!> Oct. 2019 - added lbmont
!> Oct. 2019 - added rmunvu and rmunvv, and layer 1 nested velocity ranges
!> Nov. 2019 - added amoflg
