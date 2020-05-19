#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif
      module mod_hycom
#if defined(USE_ESMF4)
      use ESMF_Mod       ! ESMF  Framework
#endif
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
      use mod_pipe       ! HYCOM debugging interface
      use mod_incupd     ! HYCOM incremental update (for data assimilation)
      use mod_floats     ! HYCOM synthetic floats, drifters and moorings
      use mod_tides      ! HYCOM tides
      use mod_archiv     ! HYCOM archives
      use mod_mean       ! HYCOM mean archives
      use mod_momtum     ! HYCOM momentum
      use mod_tsadvc     ! HYCOM scalar advection
      use mod_barotp     ! HYCOM barotropic
      use mod_asselin    ! HYCOM Asselin filter
      use mod_restart    ! HYCOM restart
#if defined(STOKES)
      use mod_stokes     ! HYCOM Stokes drift
#endif

!
! --- -----------------------------------------
! --- MICOM-based HYbrid Coordinate Ocean Model
! ---               H Y C O M
! ---           v e r s i o n  2.2
! --- -----------------------------------------
!
      implicit none
!
#if defined(USE_ESMF4)
      public HYCOM_SetServices
#else
      public HYCOM_Init, HYCOM_Run, HYCOM_Final
#endif
!
      logical, save, public  :: put_export   !set in main program
      logical, save, public  :: get_import   !set in main program
      logical, save, public  :: end_of_run   !set in HYCOM_Run
      integer, save, public  :: nts_day      !set in HYCOM_init, timesteps/day
      integer, save, public  :: nts_ice      !set in HYCOM_init, timesteps/ice
!
      integer, save, private :: &
               m,n
      real*8,  save, private :: &
              d1,d2,d3,d4,d2a,d3a,d4a, &
              ddsurf,ddiagf,dproff,dtilef,drstrf,dmeanf, &
              dske,dskea,dsmr,dsmra,dsms,dsmsa,dsmt,dsmta,dsum,dsuma, &
              dbot, &
#if defined(STOKES)
              dskes,dskesa, &
#endif
              dsumtr(mxtrcr), &
              dtime,dtime0,dbimon,dmonth,dyear,dyear0, &
              dsmall,dsmall2
      real,    save, private :: &
              smr,sms,smsa,smt,smx,sum,smin,smax, tsur, &
              coord,day1,day2,x,x1,time0,timav,cold,utotp,vtotp
      real,    save, private, allocatable :: &
              sminy(:),smaxy(:)
      integer, save, private :: &
              nstep0,nod, &
              jja, &
              lt,ma0,ma1,ma2,ma3,mc0,mc1,mc2,mc3, &
                 mk0,mk1,mk2,mk3,mr0,mr1,mr2,mr3,mnth,iofl, &
              jday,ihour,iyear
      logical, save, private :: &
               linit,diagsv,hisurf,histry,hiprof,hitile,histmn, &
               restrt,diag_tide
      character, save, private :: &
               intvl*3,c_ydh*14
!
      real*8     hours1,days1,days6
      private    hours1,days1,days6
      parameter (hours1=1.d0/24.d0,days1=1.d0,days6=6.d0)
!
! --- tfrz_n = nominal ice melting point (degC) for ice mask
      real       tfrz_n
      private    tfrz_n
      parameter (tfrz_n=-1.79)  !slightly above -1.8
!
#if defined (USE_NUOPC_CESMBETA)
! --- ntavg number of time step between coupling sequence (CESM)
      integer :: ntavg, nstep2_cpl
#endif
#if defined (ESPC_COUPLE)
! ESPC -- add
      real, save ::  nstep1_cpl,nstep2_cpl
#endif
      logical, save, public  :: end_of_run_cpl !set in HYCOM_Run for coupling

#if defined(USE_ESMF4)
!
! --- Data types for Import/Export array pointers
      type ArrayPtrReal2D
        real(ESMF_KIND_R4), dimension(:,:), pointer :: p
      end type ArrayPtrReal2D
!
! --- Attribute names for fields
      character(ESMF_MAXSTR), save :: &
          attNameLongName = "long_name", &
          attNameStdName  = "standard_name", &
          attNameUnits    = "units", &
          attNameSclFac   = "scale_factor", &
          attNameAddOff   = "add_offset"
!
! --- Import Fields
      integer, parameter :: numImpFields=11
      character(ESMF_MAXSTR), save :: impFieldName(    numImpFields), &
                                      impFieldLongName(numImpFields), &
                                      impFieldStdName( numImpFields), &
                                      impFieldUnits(   numImpFields)
      real(ESMF_KIND_R4),     save :: impFieldSclFac(  numImpFields), &
                                      impFieldAddOff(  numImpFields)
      integer,                save :: impFieldHalo(    numImpFields)
!
! --- Export Fields
      integer, parameter :: numExpFields=7
      character(ESMF_MAXSTR), save :: expFieldName(    numExpFields), &
                                      expFieldLongName(numExpFields), &
                                      expFieldStdName( numExpFields), &
                                      expFieldUnits(   numExpFields)
      real(ESMF_KIND_R4),     save :: expFieldSclFac(  numExpFields), &
                                      expFieldAddOff(  numExpFields)
      integer,                save :: expFieldHalo(    numExpFields)
!
! --- ESMF related variables
      type(ESMF_FieldBundle), save :: expBundle, &
                                      impBundle
      type(ESMF_Field),       save :: expField(numExpFields), &
                                      impField(numImpFields)
      type(ArrayPtrReal2D),   save :: expData( numExpFields), &
                                      impData( numImpFields)
!
      type(ESMF_Clock),       save :: intClock
      type(ESMF_VM),          save :: vm
      type(ESMF_DELayout),    save :: deLayout
      integer,                save :: petCount, localPet, &
                                      mpiCommunicator
      type(ESMF_Grid),        save :: grid2D
      type(ESMF_DistGrid),    save :: distgrid2D
      type(ESMF_ArraySpec),   save :: arraySpec2Dr
!
      real, save, allocatable, dimension (:,:) :: &
             sic_import  & !Sea Ice Concentration
      ,     sitx_import  & !Sea Ice X-Stress
      ,     sity_import  & !Sea Ice Y-Stress
      ,     siqs_import  & !Solar Heat Flux thru Ice to Ocean
      ,     sifh_import  & !Ice Freezing/Melting Heat Flux
      ,     sifs_import  & !Ice Freezing/Melting Salt Flux
      ,     sifw_import  & !Ice Net Water Flux
      ,      sit_import  & !Sea Ice Temperature
      ,      sih_import  & !Sea Ice Thickness
      ,      siu_import  & !Sea Ice X-Velocity
      ,      siv_import  & !Sea Ice Y-Velocity
      ,      ocn_mask      !Ocean Currents Mask
      logical, save :: &
             ocn_mask_init
#endif

      contains

#if defined(USE_ESMF4)
      subroutine HYCOM_SetServices(gridComp, rc)
!
      type(ESMF_GridComp)  :: gridComp
      integer, intent(out) :: rc
!
      call ESMF_GridCompSetEntryPoint( &
           gridComp, &
           ESMF_SETINIT, &
           HYCOM_Init, &
           ESMF_SINGLEPHASE, &
           rc=rc)
      call ESMF_GridCompSetEntryPoint( &
           gridComp, &
           ESMF_SETRUN, &
           HYCOM_Run, &
           ESMF_SINGLEPHASE, &
           rc=rc)
      call ESMF_GridCompSetEntryPoint( &
           gridComp, &
           ESMF_SETFINAL, &
           HYCOM_Final, &
           ESMF_SINGLEPHASE, &
           rc=rc)
!
      end subroutine HYCOM_SetServices

      subroutine Setup_ESMF(gridComp, impState, expState, extClock, rc)
!
! --- Calling parameters
      type(ESMF_GridComp)  :: gridComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
!
! --- set up ESMF data structures for HYCOM.
!
      real(ESMF_KIND_R4),pointer :: Xcoord(:,:),Ycoord(:,:)
      integer,           pointer :: mask_ptr(:,:)
      integer                    :: i,j,rc2
      integer                    :: lbnd(2),ubnd(2)
      character(10)              :: dimNames(2),dimUnits(2)
      type(ESMF_Logical)         :: periodic(2)
      integer(ESMF_KIND_I4)      :: year,month,day,hour,minute
      integer(ESMF_KIND_I4)      :: sec,msec,usec,nsec
      real(8)                    :: dsec,dmsec,dusec,dnsec
      type(ESMF_TimeInterval)    :: timeStep, runDuration
      type(ESMF_Time)            :: startTime
      character(ESMF_MAXSTR)     :: msg, gridName
!
! --- Report
      call ESMF_LogWrite("HYCOM Setup routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
!  Attributes for import fields, identical to CICE export fields
      impFieldAddOff(:) = 0.0     !default is no offset
      impFieldSclFac(:) = 1.0     !default is no scale factor
      impFieldHalo(  :) = halo_ps !default is scalar p-grid
!
      impFieldName(     1) = "sic"
      impFieldLongName( 1) = "Sea Ice Concentration"
      impFieldStdName(  1) = "sea_ice_area_fraction"
      impFieldUnits(    1) = "1"
      impFieldName(     2) = "sitx"
      impFieldLongName( 2) = "Sea Ice X-Stress"
      impFieldStdName(  2) = "downward_x_stress_at_sea_ice_base"
      impFieldSclFac(   2) = -1.0  !field is upward
      impFieldUnits(    2) = "Pa"
      impFieldHalo(     2) = halo_pv !vector p-grid
      impFieldName(     3) = "sity"
      impFieldLongName( 3) = "Sea Ice Y-Stress"
      impFieldStdName(  3) = "downward_y_stress_at_sea_ice_base"
      impFieldSclFac(   3) = -1.0  !field is upward
      impFieldUnits(    3) = "Pa"
      impFieldHalo(     3) = halo_pv !vector p-grid
      impFieldName(     4) = "siqs"
      impFieldLongName( 4) = "Solar Heat Flux thru Ice to Ocean"
      impFieldStdName(  4) = "downward_sea_ice_basal_solar_heat_flux"
      impFieldUnits(    4) = "W m-2"
      impFieldName(     5) = "sifh"
      impFieldLongName( 5) = "Ice Freezing/Melting Heat Flux"
      impFieldStdName(  5) = "upward_sea_ice_basal_heat_flux"
      impFieldSclFac(   5) = -1.0  !field is downward
      impFieldUnits(    5) = "W m-2"
      impFieldName(     6) = "sifs"
      impFieldLongName( 6) = "Ice Freezing/Melting Salt Flux"
      impFieldStdName(  6) = "downward_sea_ice_basal_salt_flux"
      impFieldUnits(    6) = "kg m-2 s-1"
      impFieldName(     7) = "sifw"
      impFieldLongName( 7) = "Ice Net Water Flux"
      impFieldStdName(  7) = "downward_sea_ice_basal_water_flux"
      impFieldUnits(    7) = "kg m-2 s-1"
      impFieldName(     8) = "sit" !diagnostic
      impFieldLongName( 8) = "Sea Ice Temperature"
      impFieldStdName(  8) = "sea_ice_temperature"
      impFieldAddOff(   8) = +273.15 !field is in degC
      impFieldUnits(    8) = "K"
      impFieldName(     9) = "sih" !diagnostic
      impFieldLongName( 9) = "Sea Ice Thickness"
      impFieldStdName(  9) = "sea_ice_thickness"
      impFieldUnits(    9) = "m"
      impFieldName(    10) = "siu" !diagnostic
      impFieldLongName(10) = "Sea Ice X-Velocity"
      impFieldStdName( 10) = "sea_ice_x_velocity"
      impFieldUnits(   10) = "m s-1"
      impFieldHalo(    10) = halo_pv !vector p-grid
      impFieldName(    11) = "siv" !diagnostic
      impFieldLongName(11) = "Sea Ice Y-Velocity"
      impFieldStdName( 11) = "sea_ice_y_velocity"
      impFieldUnits(   11) = "m s-1"
      impFieldHalo(    11) = halo_pv !vector p-grid
!
!     impFieldName(    12) = "patm"
!     impFieldLongName(12) = "Surface Air Pressure"
!     impFieldStdName( 12) = "surface_air_pressure"
!     impFieldUnits(   12) = "Pa"
!     impFieldName(    13) = "xwnd"
!     impFieldLongName(13) = "X-Wind"
!     impFieldStdName( 13) = "x_wind"
!     impFieldUnits(   13) = "m s-1"
!     impFieldHalo(    13) = halo_pv !vector p-grid
!     impFieldName(    14) = "ywnd"
!     impFieldLongName(14) = "Y-Wind"
!     impFieldStdName( 14) = "y_wind"
!     impFieldUnits(   14) = "m s-1"
!     impFieldHalo(    14) = halo_pv !vector p-grid
!
!  Attributes for export fields, identical to CICE import fields
      expFieldAddOff(:) = 0.0 !default is no offset
      expFieldSclFac(:) = 1.0 !default is no scale factor
      expFieldHalo(  :) = halo_ps !default is scalar p-grid
!
      expFieldName(     1) = "sst"
      expFieldLongName( 1) = "Sea Surface Temperature"
      expFieldStdName(  1) = "sea_surface_temperature"
      expFieldAddOff(   1) = +273.15 !field is in degC
      expFieldUnits(    1) = "K"
      expFieldName(     2) = "sss"
      expFieldLongName( 2) = "Sea Surface Salinity"
      expFieldStdName(  2) = "sea_surface_salinity"
      expFieldUnits(    2) = "1e-3"
      expFieldName(     3) = "ssu"
      expFieldLongName( 3) = "Sea Surface X-Current"
      expFieldStdName(  3) = "sea_water_x_velocity"
      expFieldUnits(    3) = "m s-1"
      expFieldHalo(     3) = halo_pv !vector p-grid
      expFieldName(     4) = "ssv"
      expFieldLongName( 4) = "Sea Surface Y-Current"
      expFieldStdName(  4) = "sea_water_y_velocity"
      expFieldUnits(    4) = "m s-1"
      expFieldHalo(     4) = halo_pv !vector p-grid
      expFieldName(     5) = "ssh"
      expFieldLongName( 5) = "Sea Surface Height"
      expFieldStdName(  5) = "sea_surface_height_above_sea_level"
      expFieldUnits(    5) = "m"
      expFieldName(     6) = "ssfi"
      expFieldLongName( 6) = "Oceanic Heat Flux Available to Sea Ice"
      expFieldStdName(  6) = "upward_sea_ice_basal_available_heat_flux"
      expFieldSclFac(   6) = -1.0  !field is downward
      expFieldUnits(    6) = "W m-2"
      expFieldName(     7) = "mlt"  !diagnostic
      expFieldLongName( 7) = "Ocean Mixed Layer Thickness"
      expFieldStdName(  7) = "ocean_mixed_layer_thickness"
      expFieldUnits(    7) = "m"
!
!  Create a DE layout to match HYCOM layout
      deLayout = ESMF_DELayoutCreate(vm, rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "Setup_ESMF: DELayoutCreate failed", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
!  Create array specifications
      call ESMF_ArraySpecSet(arraySpec2Dr, &
           rank=2, &
           typekind=ESMF_TYPEKIND_R4, &
           rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "Setup_ESMF: ArraySpecSet failed", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
!  Create an ESMF grid that matches the HYCOM 2D grid
      dimNames(1)="longitude";    dimNames(2)="latitude";
      dimUnits(1)="degrees_east"; dimUnits(2)="degrees_north";
      periodic(1)=ESMF_TRUE
      periodic(2)=ESMF_FALSE
!
! make a dist grid object using the deBlockList defined in mod_xc_mp.h
#if defined(ARCTIC)
! --- Arctic (tripole) domain, top row is replicated (ignore it)
      distgrid2D=ESMF_DistGridCreate(minIndex=(/   1,     1/), &
                                     maxIndex=(/itdm,jtdm-1/), &
                                   indexflag=ESMF_INDEX_GLOBAL, &
                                 deBlockList=deBlockList(:,:,1:ijpr), &
                                          rc=rc)
#else
      distgrid2D=ESMF_DistGridCreate(minIndex=(/   1,   1/), &
                                     maxIndex=(/itdm,jtdm/), &
                                   indexflag=ESMF_INDEX_GLOBAL, &
                                 deBlockList=deBlockList(:,:,1:ijpr), &
                                          rc=rc)
#endif
      if (ESMF_LogMsgFoundError(rc, &
         "Setup_ESMF: ESMF_DistGridCreate", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
! Now create a 2D grid
      grid2D=ESMF_GridCreate(distGrid=distGrid2D, &
                             coordTypeKind=ESMF_TYPEKIND_R4, &
                             coordDimCount=(/2,2/), &
                             indexflag=ESMF_INDEX_GLOBAL, &
                             rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridCreate failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
! Add Grid Coordinates
      call ESMF_GridAddCoord(grid=grid2D, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridAddCoord failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
       call ESMF_GridGetCoord(grid=grid2D, &
                              CoordDim=1, &
                              localDe=0, &
                              staggerloc=ESMF_STAGGERLOC_CENTER, &
                              computationalLbound=lbnd, &
                              computationalUbound=ubnd, &
                              fptr=Xcoord, &
                              rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridGetCoord-1 failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
#if defined(ARCTIC)
! --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif
      do j= 1,jja
        do i= 1,ii
          Xcoord(i+i0,j+j0) = plon(i,j)
        enddo
      enddo
      call ESMF_GridGetCoord(grid=grid2D, &
                             CoordDim=2, &
                             localDe=0, &
                             staggerloc=ESMF_STAGGERLOC_CENTER, &
                             computationalLbound=lbnd, &
                             computationalUbound=ubnd, &
                             fptr=Ycoord, &
                             rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridGetCoord-2 failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
      do j= 1,jja
        do i= 1,ii
          Ycoord(i+i0,j+j0) = plat(i,j)
        enddo
      enddo
      CALL ESMF_GridAddItem(grid=grid2D, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            item=ESMF_GRIDITEM_MASK, &
                            rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridAddItem failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
      CALL ESMF_GridGetItem(grid=grid2D, &
                            localDE=0, &
                            staggerloc=ESMF_STAGGERLOC_CENTER, &
                            item=ESMF_GRIDITEM_MASK, &
                            fptr=mask_ptr, &
                            rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "GridGetItem failed",rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
      mask_ptr(:,:) = 0  !all land, outside active tile
      do j= 1,jja
        do i= 1,ii
          mask_ptr(i+i0,j+j0) = ishlf(i,j)
        enddo
      enddo
!
!  Associate grid with ESMF gridded component
      call ESMF_GridCompSet(gridComp, grid=grid2D, rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "Setup_ESMF: GridCompSet", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
!  Setup export fields, bundles & state
      do i=1,numExpFields
        expField(i)=ESMF_FieldCreate(grid=grid2D, &
                                     arrayspec=arraySpec2Dr, &
                                     indexflag=ESMF_INDEX_GLOBAL, &
                                     staggerLoc=ESMF_STAGGERLOC_CENTER, &
                                     name=trim(expFieldName(i)), &
                                     rc=rc)
        call ESMF_FieldGet(expField(i),0,expData(i)%p,rc=rc)
        expData(i)%p(:,:) = 0.0
      enddo
!
!  Create bundle from list of fields
       expBundle=ESMF_FieldBundleCreate(numExpFields, &
                        expField(:), name='HYCOM Export', &
                        rc=rc)
!
!  Add bundle to the export state
       call ESMF_StateAdd(expState,expBundle,rc=rc)
!
!  Setup import fields, bundles & state
      do i = 1,numImpFields
        impField(i)=ESMF_FieldCreate(grid2D,arraySpec2Dr, &
                                     StaggerLoc=ESMF_STAGGERLOC_CENTER, &
                                     name=trim(impFieldName(i)), &
                                     rc=rc)
        call ESMF_FieldGet(impField(i),0,impData(i)%p,rc=rc)
        impData(i)%p(:,:)=0.0 ! Initialize fields
      enddo
!
      allocate(  sic_import(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
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
!
       sic_import(:,:) = 0.0 !Sea Ice Concentration
      sitx_import(:,:) = 0.0 !Sea Ice X-Stress
      sity_import(:,:) = 0.0 !Sea Ice Y-Stress
      siqs_import(:,:) = 0.0 !Solar Heat Flux thru Ice to Ocean
      sifh_import(:,:) = 0.0 !Ice Freezing/Melting Heat Flux
      sifs_import(:,:) = 0.0 !Ice Freezing/Melting Salt Flux
      sifw_import(:,:) = 0.0 !Ice Net Water Flux
       sit_import(:,:) = 0.0 !Sea Ice Temperature
       sih_import(:,:) = 0.0 !Sea Ice Thickness
       siu_import(:,:) = 0.0 !Sea Ice X-Velocity
       siv_import(:,:) = 0.0 !Sea Ice Y-Velocity
!
!  Create bundle from list of fields
      impBundle=ESMF_FieldBundleCreate(numImpFields, &
                                       impField(:), &
                                       name='HYCOM Import', &
                                       rc=rc)
!
!  Add bundle to the import state
      call ESMF_StateAdd(impState,impBundle,rc=rc)
!
      ocn_mask_init = .true.  !still need to initialize ocn_mask
!
      end subroutine Setup_ESMF

      subroutine Export_ESMF()
!
! --- Fill export state.
! --- Calculate ssfi "in place"
!
      integer i,j,k,rc
      real    ssh2m
      real    tmxl,smxl,umxl,vmxl,hfrz,tfrz,t2f,ssfi
      real    dp1,usur1,vsur1,psur1,dp2,usur2,vsur2,psur2,thksur
!
! --- Report
      call ESMF_LogWrite("HYCOM Export routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
      if     (ocn_mask_init) then  !very 1st call to this routine
        ocn_mask_init = .false.
!
        allocate( ocn_mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy) )
!
        if     (iceflg.eq.4) then
          ocn_mask(:,:) = 0.0  !export ocean currents nowhere
        elseif (nestfq.ne.0.0) then
!         export ocean currents away from open boundaries
          do j= 1,jj
            do i= 1,ii
              if     (rmunv(i,j).ne.0.0) then
                ocn_mask(i,j) = 0.0
              else
                ocn_mask(i,j) = 1.0
              endif
            enddo !i
          enddo !j
          do i= 1,10
            call psmooth(ocn_mask,0,0, ip, util1)  !not efficient, but only done once
          enddo !i
        else
          ocn_mask(:,:) = 1.0  !export ocean currents everywhere
        endif
      endif !ocn_mask_init
!
! --- Assume Export State is as defined in Setup_ESMF
! --- Average two time levels since (the coupling frequency) icpfrq >> 2
!
      call xctilr(u(      1-nbdy,1-nbdy,1,1),1,2*kk, 1,1, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,1),1,2*kk, 1,1, halo_vv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,  1),1,   2, 1,1, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,  1),1,   2, 1,1, halo_vv)
!
      ssh2m = 1.0/g
      do j= 1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
! ---       quantities for available freeze/melt heat flux
! ---       relax to tfrz with e-folding time of icefrq time steps
! ---       assuming the effective surface layer thickness is hfrz
! ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
            hfrz = min( thkfrz*onem, dpbl(i,j) )
            t2f  = (spcifh*hfrz)/(baclin*icefrq*g)
! ---       average both available time steps, to avoid time splitting.
            smxl = 0.5*(saln(i,j,1,n)+saln(i,j,1,m))
            tmxl = 0.5*(temp(i,j,1,n)+temp(i,j,1,m))
            tfrz = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
            ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
! ---       average currents over top thkcdw meters
            thksur = onem*min( thkcdw, depths(i,j) ) 
            usur1  = 0.0
            vsur1  = 0.0
            psur1  = 0.0
            usur2  = 0.0
            vsur2  = 0.0
            psur2  = 0.0
            do k= 1,kk
              dp1   = min( dp(i,j,k,1), max( 0.0, thksur-psur1 ) )
              usur1 = usur1 + dp1*(u(i,j,k,1)+u(i+1,j,k,1))
              vsur1 = vsur1 + dp1*(v(i,j,k,1)+v(i,j+1,k,1))
#if defined(STOKES)
              usur1 = usur1 + dp1*(usd(i,j,k)+usd(i+1,j,k))
              vsur1 = vsur1 + dp1*(vsd(i,j,k)+vsd(i,j+1,k))
#endif
              psur1 = psur1 + dp1
              dp2   = min( dp(i,j,k,2), max( 0.0, thksur-psur2 ) )
              usur2 = usur2 + dp2*(u(i,j,k,2)+u(i+1,j,k,2))
              vsur2 = vsur2 + dp2*(v(i,j,k,2)+v(i,j+1,k,2))
#if defined(STOKES)
              usur2 = usur2 + dp2*(usd(i,j,k)+usd(i+1,j,k))
              vsur2 = vsur2 + dp2*(vsd(i,j,k)+vsd(i,j+1,k))
#endif
              psur2 = psur2 + dp2
              if     (min(psur1,psur2).ge.thksur) then
                exit
              endif
            enddo
            umxl  = 0.25*( usur1/psur1 + ubavg(i,  j,1) + &
                                         ubavg(i+1,j,1) + &
                           usur2/psur2 + ubavg(i,  j,2) + &
                                         ubavg(i+1,j,2)  )
            vmxl  = 0.25*( vsur1/psur1 + vbavg(i,j,  1) + &
                                         vbavg(i,j+1,1) + &
                           vsur2/psur2 + vbavg(i,j,  2) + &
                                         vbavg(i,j+1,2)  )
            util2(i,j)              = umxl
            util3(i,j)              = vmxl
            expData(1)%p(i+i0,j+j0) = tmxl
            expData(2)%p(i+i0,j+j0) = smxl
            expData(5)%p(i+i0,j+j0) = ssh2m*srfhgt(i,j)  !ssh in m
            expData(6)%p(i+i0,j+j0) = max(-1000.0,min(1000.0,ssfi))  !as in CICE
            expData(7)%p(i+i0,j+j0) = dpbl(i,j)*qonem
          else
            util2(i,j)              = 0.0
            util3(i,j)              = 0.0
          endif !ishlf:else
        enddo !i
      enddo !j
!
! --- Smooth surface ocean velocity fields
#if defined(ARCTIC)
      call xctila( util2,1,1,halo_pv)
      call xctila( util3,1,1,halo_pv)
#endif
      call psmooth(util2,0,0, ishlf, util1)
      call psmooth(util3,0,0, ishlf, util1)
      do j= 1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
            expData(3)%p(i+i0,j+j0) = util2(i,j)*ocn_mask(i,j)
            expData(4)%p(i+i0,j+j0) = util3(i,j)*ocn_mask(i,j)
          endif !ishlf
        enddo !i
      enddo !j
!
! --- Report
      call ESMF_LogWrite("HYCOM Export routine returned", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
      end subroutine Export_ESMF

      subroutine Import_ESMF()
!
! --- Extract import state.
!
      integer i,j,rc
!
! --- Report
      call ESMF_LogWrite("HYCOM Import routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
! --- Assume Import State is as defined in Setup_ESMF
!
#if defined(ARCTIC)
! --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif
      do j= 1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
             sic_import(i,j) = impData( 1)%p(i+i0,j+j0) !Sea Ice Concentration
            sitx_import(i,j) = impData( 2)%p(i+i0,j+j0) !Sea Ice X-Stress
            sity_import(i,j) = impData( 3)%p(i+i0,j+j0) !Sea Ice Y-Stress
            siqs_import(i,j) = impData( 4)%p(i+i0,j+j0) !Solar Heat Flux thru Ice to Ocean
            sifh_import(i,j) = impData( 5)%p(i+i0,j+j0) !Ice Freezing/Melting Heat Flux
            sifs_import(i,j) = impData( 6)%p(i+i0,j+j0) !Ice Freezing/Melting Salt Flux
            sifw_import(i,j) = impData( 7)%p(i+i0,j+j0) !Ice Net Water Flux
             sit_import(i,j) = impData( 8)%p(i+i0,j+j0) !Sea Ice Temperature
             sih_import(i,j) = impData( 9)%p(i+i0,j+j0) !Sea Ice Thickness
             siu_import(i,j) = impData(10)%p(i+i0,j+j0) !Sea Ice X-Velocity
             siv_import(i,j) = impData(11)%p(i+i0,j+j0) !Sea Ice Y-Velocity
            if     (iceflg.ge.2 .and. icmflg.ne.3) then
              covice(i,j) = impData(1)%p(i+i0,j+j0) !Sea Ice Concentration
                si_c(i,j) = impData(1)%p(i+i0,j+j0) !Sea Ice Concentration
              if     (covice(i,j).gt.0.0) then
                 si_tx(i,j) = -impData( 2)%p(i+i0,j+j0) !Sea Ice X-Stress into ocean
                 si_ty(i,j) = -impData( 3)%p(i+i0,j+j0) !Sea Ice Y-Stress into ocean
                fswice(i,j) =  impData( 4)%p(i+i0,j+j0) !Solar Heat Flux thru Ice to Ocean
                flxice(i,j) =  fswice(i,j) + &
                               impData( 5)%p(i+i0,j+j0) !Ice Freezing/Melting Heat Flux
                sflice(i,j) =  impData( 6)%p(i+i0,j+j0)*1.e3
                                                        !Ice Freezing/Melting Salt Flux
                wflice(i,j) =  impData( 7)%p(i+i0,j+j0) !Ice Water Flux
                temice(i,j) =  impData( 8)%p(i+i0,j+j0) !Sea Ice Temperature
                  si_t(i,j) =  impData( 8)%p(i+i0,j+j0) !Sea Ice Temperature
                thkice(i,j) =  impData( 9)%p(i+i0,j+j0) !Sea Ice Thickness
                  si_h(i,j) =  impData( 9)%p(i+i0,j+j0) !Sea Ice Thickness
                  si_u(i,j) =  impData(10)%p(i+i0,j+j0) !Sea Ice X-Velocity
                  si_v(i,j) =  impData(11)%p(i+i0,j+j0) !Sea Ice Y-Velocity
              else
                 si_tx(i,j) = 0.0
                 si_ty(i,j) = 0.0
                fswice(i,j) = 0.0
                flxice(i,j) = 0.0
                sflice(i,j) = 0.0
                wflice(i,j) = 0.0
                temice(i,j) = 0.0
                  si_t(i,j) = 0.0
                thkice(i,j) = 0.0
                  si_h(i,j) = 0.0
                  si_u(i,j) = 0.0
                  si_v(i,j) = 0.0
              endif !covice
            elseif (iceflg.ge.2 .and. icmflg.eq.3) then
                  si_c(i,j) =  impData( 1)%p(i+i0,j+j0) !Sea Ice Concentration
              if     (si_c(i,j).gt.0.0) then
                 si_tx(i,j) = -impData( 2)%p(i+i0,j+j0) !Sea Ice X-Stress into ocean
                 si_ty(i,j) = -impData( 3)%p(i+i0,j+j0) !Sea Ice Y-Stress into ocean
                  si_h(i,j) =  impData( 9)%p(i+i0,j+j0) !Sea Ice Thickness
                  si_t(i,j) =  impData( 8)%p(i+i0,j+j0) !Sea Ice Temperature
                  si_u(i,j) =  impData(10)%p(i+i0,j+j0) !Sea Ice X-Velocity
                  si_v(i,j) =  impData(11)%p(i+i0,j+j0) !Sea Ice Y-Velocity
              else
                 si_tx(i,j) = 0.0
                 si_ty(i,j) = 0.0
                  si_h(i,j) = 0.0
                  si_t(i,j) = 0.0
                  si_u(i,j) = 0.0
                  si_v(i,j) = 0.0
              endif !covice
            endif !iceflg>=2 (icmflg)
          endif !ishlf
        enddo !i
      enddo !j
#if defined(ARCTIC)
!
! --- update last active row of array
      call xctila( sic_import,1,1,halo_ps)  !Sea Ice Concentration
      call xctila(sitx_import,1,1,halo_pv)  !Sea Ice X-Stress
      call xctila(sity_import,1,1,halo_pv)  !Sea Ice Y-Stress
      call xctila(siqs_import,1,1,halo_ps)  !Solar Heat Flux thru Ice to Ocean
      call xctila(sifh_import,1,1,halo_ps)  !Ice Freezing/Melting Heat Flux
      call xctila(sifs_import,1,1,halo_ps)  !Ice Freezing/Melting Salt Flux
      call xctila(sifw_import,1,1,halo_ps)  !Ice Net Water Flux
      call xctila( sit_import,1,1,halo_ps)  !Sea Ice Temperature
      call xctila( sih_import,1,1,halo_ps)  !Sea Ice Thickness
      call xctila( siu_import,1,1,halo_pv)  !Sea Ice X-Velocity
      call xctila( siv_import,1,1,halo_pv)  !Sea Ice Y-Velocity
      if     (iceflg.ge.2 .and. icmflg.ne.3) then
        call xctila(covice,1,1,halo_ps)  !Sea Ice Concentration
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(fswice,1,1,halo_ps)  !Solar Heat Flux thru Ice to Ocean
        call xctila(flxice,1,1,halo_ps)  !Sea Ice Freezing/Melting Heat Flux
        call xctila(sflice,1,1,halo_ps)  !Sea Ice Freezing/Melting Salt Flux
        call xctila(wflice,1,1,halo_ps)  !Sea Ice Water Flux
        call xctila(temice,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(thkice,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      elseif (iceflg.ge.2 .and. icmflg.eq.3) then
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      endif
#endif
!
! --- Smooth Sea Ice velocity fields
      call psmooth(si_u,0,0, ishlf, util1)
      call psmooth(si_v,0,0, ishlf, util1)
#if defined(ARCTIC)
      call xctila(si_u,1,1,halo_pv)
      call xctila(si_v,1,1,halo_pv)
#endif
! --- copy back from si_ to impData for Archive_ESMF
      do j= 1,jja
        do i= 1,ii
          if     (si_c(i,j).gt.0.0) then
            impData(10)%p(i+i0,j+j0) = si_u(i,j) !Sea Ice X-Velocity
            impData(11)%p(i+i0,j+j0) = si_v(i,j) !Sea Ice X-Velocity
          endif !si_c
        enddo !i
      enddo !j
!
! --- Report
      call ESMF_LogWrite("HYCOM Import routine returned", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
      end subroutine Import_ESMF

      subroutine Archive_ESMF(iyear,jday,ihour)
      integer   iyear,jday,ihour
!
! --- Create a HYCOM "archive-like" file from Import/Export state.
! --- Import state may not be at the same time as Export.
! --- Ice Drift has been smoothed since import.
!
      logical      hycom_isnaninf  !function to detect NaN and Inf
!
      character*8  cname
      character*80 cfile
      logical      lexist
      integer      i,j,k,nop,nopa,rc
      real         coord,xmin,xmax,sumssu,sumssv,sumsiu,sumsiv
!
! --- Report
      call ESMF_LogWrite("HYCOM Archive routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
      write(cfile,'(a,i4.4,a1,i3.3,a1,i2.2)') &
         'arche.',iyear,'_',jday,'_',ihour
      nopa=13
      nop =13+uoff
!
! --- Only write out one archive per hour
!
      inquire(file=trim(cfile)//'.a',exist=lexist)
      if     (lexist) then
!
! ---   Report
        call ESMF_LogWrite("HYCOM Archive routine returned early", &
             ESMF_LOG_INFO, rc=rc)
        call ESMF_LogFlush(rc=rc)
!
        if (mnproc.eq.1) then
        write(lp,*) 'skip: ',trim(cfile)
        call flush(lp)
        endif !1st tile
        return
      else
        if (mnproc.eq.1) then
        write(lp,*) 'open: ',trim(cfile)
        call flush(lp)
        endif !1st tile
      endif
      call xcsync(flush_lp)
!
      call zaiopf(trim(cfile)//'.a', 'new', nopa)
      if     (mnproc.eq.1) then
      open (unit=nop,file=trim(cfile)//'.b',status='new') !uoff+13
      write(nop,116) ctitle,iversn,iexpt,yrflag,itdm,jtdm
      call flush(nop)
      endif !1st tile
 116  format (a80/a80/a80/a80/ &
       i5,4x,'''iversn'' = hycom version number x10'/ &
       i5,4x,'''iexpt '' = experiment number x10'/ &
       i5,4x,'''yrflag'' = days in year flag'/ &
       i5,4x,'''idm   '' = longitudinal array size'/ &
       i5,4x,'''jdm   '' = latitudinal  array size'/ &
       'field       time step  model day', &
       '  k  dens        min              max')
!
! --- surface fields
!
      coord=0.0
      do k= 1,numExpFields
        do j= 1,jja
          do i= 1,ii
            if     (ishlf(i,j).eq.1) then
              util1(i,j) = expData(k)%p(i+i0,j+j0)
            endif !ishlf
          enddo !i
        enddo !j
#if defined(ARCTIC)
        call xctila(util1,1,1,expFieldHalo(k))
#endif
        cname = expFieldName(k)(1:8)
        call zaiowr(util1,ishlf,.true., &
                    xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      enddo !k
      do j= 1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
            util2(i,j) = impData(1)%p(i+i0,j+j0)  !ice concentration
          else
            util2(i,j) = 0.0
          endif !ishlf
        enddo !i
      enddo !j
#if defined(ARCTIC)
      call xctila(util2,1,1,halo_ps)
#endif
      cname = impFieldName(1)(1:8)
      call zaiowr(util2,ishlf,.true., &
                  xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      do k= 2,3 !si_tx,si_ty
        do j= 1,jja
          do i= 1,ii
            if     (util2(i,j).ne.0.0) then
              util1(i,j) = -impData(k)%p(i+i0,j+j0)  !into ocean
            else
              util1(i,j) = 0.0
            endif !ice:no-ice
          enddo !i
        enddo !j
#if defined(ARCTIC)
        call xctila(util1,1,1,impFieldHalo(k))
#endif
        cname = impFieldName(k)(1:8)
        call zaiowr(util1,ishlf,.true., &
                    xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) cname(1:4)//'down', &
                        nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      enddo !k
      do k= 4,7 !fluxes
        do j= 1,jja
          do i= 1,ii
            if     (util2(i,j).ne.0.0) then
              util1(i,j) = util2(i,j)*impData(k)%p(i+i0,j+j0)
            else
              util1(i,j) = hugel !mask where there is no ice
            endif !ice:no-ice
          enddo !i
        enddo !j
#if defined(ARCTIC)
        vland = hugel
        call xctila(util1,1,1,impFieldHalo(k))
        vland = 0.0
#endif
        cname = impFieldName(k)(1:8)
        call zaiowr(util1,ishlf,.false.,   & !mask on ice
                    xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      enddo !k
      do k= 8,numImpFields
        do j= 1,jja
          do i= 1,ii
            if     (util2(i,j).ne.0.0) then
              util1(i,j) = impData(k)%p(i+i0,j+j0)
            else
              util1(i,j) = hugel !mask where there is no ice
            endif !ice:no-ice
          enddo !i
        enddo !j
#if defined(ARCTIC)
        vland = hugel
        call xctila(util1,1,1,impFieldHalo(k))
        vland = 0.0
#endif
        cname = impFieldName(k)(1:8)
        call zaiowr(util1,ishlf,.false.,   & !mask on ice
                    xmin,xmax, nopa, .false.)
        if     (mnproc.eq.1) then
        write (nop,117) cname,nstep,time,0,coord,xmin,xmax
        call flush(nop)
        endif !1st tile
      enddo !k
      do j= 1,jja
        do i= 1,ii
          if     (util2(i,j).ne.0.0) then
            util1(i,j) = -impData( 7)%p(i+i0,j+j0)  !water flux
          else
            util1(i,j) = hugel !mask where there is no ice
          endif !ice:no-ice
        enddo !i
      enddo !j
#if defined(ARCTIC)
      call xctila(util1,1,1,halo_ps)
#endif
      cname = 'surtx   '
      call zaiowr(surtx,ishlf,.true., &
                  xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      cname = 'surty   '
      call zaiowr(surty,ishlf,.true., &
                  xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
      cname = 'wflice  '
      call zaiowr(util1,ishlf,.false.,   & !mask on ice
                  xmin,xmax, nopa, .false.)
      if     (mnproc.eq.1) then
      write (nop,117) cname,nstep,time,0,coord,xmin,xmax
      call flush(nop)
      endif !1st tile
 117  format (a8,' =',i11,f11.3,i3,f7.3,1p2e16.7)
!
      close (unit=nop)
      call zaiocl(nopa)
!
! --- local-tile test of velocity fields for NaNs
! --- sum should be ok unless NaNs or Infs are present
!
      sumssu = 0.0
      sumssv = 0.0
      sumsiu = 0.0
      sumsiv = 0.0
      do j= 1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
            sumssu = sumssu + expData( 3)%p(i+i0,j+j0)
            sumssv = sumssv + expData( 4)%p(i+i0,j+j0)
          endif !ishlf
          if     (util2(i,j).ne.0.0) then
            sumsiu = sumsiu + impData(10)%p(i+i0,j+j0)
            sumsiv = sumsiv + impData(11)%p(i+i0,j+j0)
          endif !ice
        enddo !i
      enddo !j
      if     (hycom_isnaninf(sumssu) .or. &
              hycom_isnaninf(sumssv) .or. &
              hycom_isnaninf(sumsiu) .or. &
              hycom_isnaninf(sumsiv)     ) then
        call xchalt('Archive_ESMF: NaN or Inf detected')
               stop 'Archive_ESMF: NaN or Inf detected'
      endif !NaN
!
! ---   Report
        call ESMF_LogWrite("HYCOM Archive routine returned", &
             ESMF_LOG_INFO, rc=rc)
        call ESMF_LogFlush(rc=rc)
!
      end subroutine Archive_ESMF
#endif /* USE_ESMF4 */

#if defined(USE_ESMF4)
      subroutine HYCOM_Init  &
                  (gridComp, impState, expState, extClock, rc)
!
! --- Calling parameters
      type(ESMF_GridComp)  :: gridComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
#else 
      subroutine HYCOM_Init &
                  (mpiCommunicator,hycom_start_dtg,hycom_end_dtg, &
                   pointer_filename)
!
      integer,      intent(in), optional :: mpiCommunicator
      real(8),         intent(in), optional :: hycom_start_dtg
      real(8),         intent(in), optional :: hycom_end_dtg
      character*80, intent(in), optional :: pointer_filename
      logical       restart_cpl
      integer :: mpi_comm_ocean,istat
#endif /* USE_ESMF4:USE_NUOPC_CESMBETA:else */
#if defined(USE_NUOPC_CESMBETA)
      real :: ssh_n,ssh_s,ssh_e,ssh_w,dhdx,dhdy
      real :: maskn,masks,maske,maskw
      real :: dp1,usur1,vsur1,psur1,dp2,usur2,vsur2,psur2,thksur, &
              utot,vtot
#endif /*USE_NUOPC_CESMBETA */
!
! --- Initialize (before the 1st time step).
!
      integer      i,j,k,nm,margin
      character*80 flnm,flnmra,flnmrb
!
# include "stmt_fns.h"
!
#if defined(USE_ESMF4)
      integer                :: rc2
      character(ESMF_MAXSTR) :: msg
!
! --- Report
      call ESMF_LogWrite("HYCOM initialize routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
! --- Get VM from gridComp
      call ESMF_GridCompGet(gridComp, vm=vm, rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "HYCOM_Init: GridCompGet failed", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
!
! --- Get VM info
      call ESMF_VMGet(vm, &
           petCount=petCount, localPET=localPet, &
           mpiCommunicator=mpiCommunicator, rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "HYCOM_Init: VMGet failed", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
      write(msg,'(a,i4)') "HYCOM_Init: petCount = ",petCount
      call ESMF_LogWrite(msg, ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
! --- initialize hycom message passing.
      call xcspmd(mpiCommunicator)
#else
! --- initialize hycom message passing.
      if (present(mpiCommunicator)) then
         call MPI_Comm_Dup(mpiCommunicator,mpi_comm_ocean,istat)
         call xcspmd(mpi_comm_ocean)
      else
! --- initialize SPMD processsing
         call xcspmd
! avoid computer warnings not used
          mpi_comm_ocean = 0
          istat = 0 
      endif
#endif
! --- initialize timer names.
!
      call xctmrn(40,'cnuity')
      call xctmrn(41,'tsadvc')
      call xctmrn(42,'momtum')
      call xctmrn(43,'barotp')
      call xctmrn(44,'thermf')
      call xctmrn(45,'ic****')
      call xctmrn(46,'mx****')
      call xctmrn(47,'conv**')
      call xctmrn(48,'diapf*')
      call xctmrn(49,'hybgen')
      call xctmrn(50,'trcupd')
      call xctmrn(51,'restrt')
      call xctmrn(52,'overtn')
      call xctmrn(53,'archiv')
      call xctmrn(54,'incupd')
      call xctmrn(55,'aslsav')
      call xctmrn(56,'asseln')
#if defined(USE_ESMF4) || defined (ESPC_COUPLE)
      call xctmrn(78,'HY_Ini')
      call xctmrn(79,'HY_Out')
      call xctmrn(80,'HY_Run')
      call xctmrs(78,huge(i),1)  !time all instances
      call xctmrs(79,huge(i),1)  !time all instances
      call xctmrs(80,huge(i),1)  !time all instances
#endif
!
! --- machine-specific initialization
      call machine

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

! --- model is to be integrated from time step 'nstep1' to 'nstep2'
! either get time from coupler or read from limits
      if ((present(hycom_start_dtg)) .and. (present(hycom_end_dtg))) then
         day1 = hycom_start_dtg
         day2 = hycom_end_dtg
      elseif ((present(hycom_start_dtg)) .or. (present(hycom_end_dtg))) then
        if (mnproc.eq.1) then
         write(lp,*) 'error in hycom - Both hycom_start_dtg and &
               &''hycom_end_dtg must be present if one is present'
         call flush(lp)
        endif !1st tile
        call xcstop('(hycom)')
               stop '(hycom)'   !won't get here
      else
         if (mnproc.eq.1) then
         write(lp,*)  trim(flnminp)//'limits'  
         call flush(lp)
        endif !1st tile
        open( unit=uoff+99,file=trim(flnminp)//'limits')
        read(      uoff+99,*) day1,day2
        close(unit=uoff+99)
      endif
! --- non-positive day1 indicates a new initialization, or
! --- the start of a yrflag==3 case.
      linit =day1.le.0.0
      day1  =abs(day1)
      dtime=day1
!
! --- initialize scalars
      call blkdat(linit)  !must call before zaiost
!
! --- initialize array i/o.
      call zaiost
      if (mnproc.eq.1) then
        call mem_stat_print('  zaiost  1st:')
      endif !1st tile
      call xcsync(flush_lp)
      if (mnproc.eq.ijpr) then
        call mem_stat_print('  zaiost last:')
      endif !last tile
      call xcsync(flush_lp)
!
! --- initialize common variables
!
      call cb_allocate
      if (mnproc.eq.1) then
        call mem_stat_print('  cb_allocate:')
      endif !1st tile
!
! --- initiate named-pipe comparison utility
      call pipe_init
!
      nstep1=nint(dtime*(86400.0d0/baclin))
      nts_ice = icpfrq                     !no. time steps between ice coupling
      nts_day = nint(86400.0d0/baclin)     !no. time steps per day
      dsmall  = baclin/86400.0d0 * 0.25d0  !1/4 of a time step in days
      dsmall2 = dsmall*2.0d0
      dtime=(nstep1/nts_day)+mod(nstep1,nts_day)*(baclin/86400.0d0)             
      day1 =dtime 
      if     (dsurfq.ge.1.0) then
        ddsurf = dsurfq
!       write(6,'("Case 1 ddsurf =",G25.17)')ddsurf
      elseif (dsurfq.ne.0.0) then
        ddsurf = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*dsurfq)/baclin))
!       write(6,'("Case 1 ddsurf =",G25.17)')ddsurf
      else !no surface archives
        ddsurf = (baclin/86400.0d0)*0.99d0*huge(i)
      endif
      if     (diagfq.ge.1.0) then
        ddiagf = diagfq
      elseif (diagfq.ne.0.0) then
        ddiagf = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*diagfq)/baclin))
      else !no 3-d archives
        ddiagf = (baclin/86400.0d0)*0.99d0*huge(i)
      endif
      if     (proffq.ge.1.0) then
        dproff = proffq
      elseif (proffq.ne.0.0) then
        dproff = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*proffq)/baclin))
      else !no 3-d profiles at selection locations
        dproff = (baclin/86400.0d0)*0.99d0*huge(i)
      endif
      if     (tilefq.ge.1.0) then
        dtilef = tilefq
      elseif (tilefq.ne.0.0) then
        dtilef = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*tilefq)/baclin))
      else !no tiled 3-d archives
        dtilef = (baclin/86400.0d0)*0.99d0*huge(i)
      endif
      if     (meanfq.ge.1.0) then
        dmeanf = meanfq
      elseif (meanfq.ne.0.0) then
        dmeanf = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*meanfq)/baclin))
      else !no mean archives
        dmeanf = (baclin/86400.0d0)*0.99d0*huge(i)
      endif
      if     (rstrfq.eq.0.0) then  ! no restart
        drstrf = rstrfq
      elseif (rstrfq.lt.0.0) then  ! no restart at end of run
        drstrf = -rstrfq
      elseif (rstrfq.ge.1.0) then
        drstrf = rstrfq
      else
        drstrf = (baclin/86400.0d0)* &
                 max(1,nint((86400.0d0*rstrfq)/baclin))
      endif
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,*) 'ddsurf = ',ddsurf,nint((86400.0d0*ddsurf)/baclin)
      write(lp,*) 'ddiagf = ',ddiagf,nint((86400.0d0*ddiagf)/baclin)
      write(lp,*) 'dproff = ',dproff,nint((86400.0d0*dproff)/baclin)
      write(lp,*) 'dtilef = ',dtilef,nint((86400.0d0*dtilef)/baclin)
      write(lp,*) 'dmeanf = ',dmeanf,nint((86400.0d0*dmeanf)/baclin)
      write(lp,*) 'drstrf = ',drstrf,nint((86400.0d0*drstrf)/baclin)
      write(lp,*)
      write (lp,101) thkdf2,temdf2, &
                     thkdf4, &
                     veldf2,visco2, &
                     veldf4,visco4, &
                     diapyc,vertmx
 101  format ( &
        ' turb. flux parameters:',1p/ &
        ' thkdf2,temdf2 =',2e9.2/ &
        ' thkdf4        =', e9.2/ &
        ' veldf2,visco2 =',2e9.2/ &
        ' veldf4,visco4 =',2e9.2/ &
        ' diapyc,vertmx =',2e9.2/)
      endif !1st tile
!
! --- days in year.
!
      if     (yrflag.eq.0) then
! ---   360 days, starting Jan 16
        dmonth =  30.0d0
        dbimon =  60.0d0
        dyear  = 360.0d0
        dyear0 =   0.0d0
      elseif (yrflag.eq.1) then
! ---   366 days, starting Jan 16
        dmonth =  30.5d0
        dbimon =  61.0d0
        dyear  = 366.0d0
        dyear0 =   0.0d0
      elseif (yrflag.eq.2) then
! ---   366 days, starting Jan 1
! ---   also implies high frequency atmospheric forcing
        dmonth =  30.5d0
        dbimon =  61.0d0
        dyear  = 366.0d0
        dyear0 = -15.0d0+dyear
      elseif (yrflag.eq.3) then
! ---   model day is calendar days since 01/01/1901
! ---   also implies high frequency atmospheric forcing
        dyear  = 365.25d0
        dmonth = dyear/12.d0
        dbimon = dyear/ 6.d0
        dyear0 = -15.0d0+dyear
      elseif (yrflag.eq.4) then
! ---   365 days, starting Jan 1 - CESM: NO-LEAP calendar
! ---   also implies high frequency atmospheric forcing
        dyear  = 365.0d0
        dmonth = dyear/12.d0
        dbimon = dyear/ 6.d0
        dyear0 = -15.0d0+dyear
      else
        if (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in hycom - unsupported yrflag value'
        write(lp,*)
        call flush(lp)
        endif !1st tile
        call xcstop('(hycom)')
               stop '(hycom)'   !won't get here
      endif
!
! --- number of baroclinic time steps per day...
      nsteps_per_day = nint(86400.0/baclin)
!
! --- initialize barotp coeflx
      call barotp_init  !!added by Alex
!
! --- set up parameters defining the geographic environment
!
      call geopar
#if defined(USE_ESMF4)
!
! --- set up ESMF data structures
!
      call Setup_ESMF(gridComp, impState, expState, extClock, rc=rc)
      if (ESMF_LogMsgFoundError(rc, &
         "HYCOM_Init: Setup_ESMF failed", rcToReturn=rc2)) &
         call ESMF_Finalize(rc=rc)
#endif

!
! --- set up forcing functions
!
      if (yrflag.lt.2) then
        call forfuna  ! monthly atmospheric forcing
      endif
      if     (jerlv0.eq. 0) then
        call forfunk  !  annual/monthly kpar
      elseif (jerlv0.eq.-1) then
        call forfunc  !  annual/monthly chl
      endif
      call forfunp  !    annual/monthly rivers
      call forfunr  ! bimonthly/monthly climatology
      watcum=0.
      empcum=0.
!
! --- set minimum salinity for each isopycnic layer
      if     (isopyc) then
        do k=2,kk
          cold=-3.0
          salmin(k)=sofsig(sigma(k),cold)
        enddo
      endif
!
! --- layer specific volume is defined as (1-theta)*svref
! --- subtract constant 'thbase' from theta to reduce roundoff errors
!
      if     (vsigma) then
        call forfunv  ! spacially varying isopycnal target densities
      else
        do k=1,kk
          theta(:,:,k)=sigma(k)-thbase
        enddo
      endif
!
! --- minimum depth of isopycnmal layers (pressure units).
!
      if     (isotop.lt.0.0) then
        call forfunt    !spacially varying minimum depths
      else
        topiso(:,:)=onem*isotop  !constant minimum depth
      endif
!
! --- tidal drag roughness (m/s)
!
      if     (drgscl.ne.0.0) then
        call forfund(tiddrg)    !tidal drag scalar or tensor
      else
        drgten(:,:,:,:)=0.0
      endif
!
! --- "scalar" tidal SAL factor
!
      if     (tidflg.eq.0) then
        salfac(:,:)=0.0     !not used, set it for safety
      elseif (tidsal.lt.0.0) then
        call forfuns        !varying tidal SAL factor
      else
        salfac(:,:)=tidsal  !scalar  tidal SAL factor
      endif
!
! --- veldf2, veldf4 and thkdf4 may be spacially varying
!
      call forfundf
!
!
#if defined (ESPC_COUPLE)
      nstep1_cpl=nstep1
      nstep2_cpl=nstep1
#endif

      dtime=day2
      nstep2=nint(dtime*(86400.0d0/baclin))
      dtime=(nstep2/nts_day)+mod(nstep2,nts_day)*(baclin/86400.0d0)
      day2 =dtime
!
      if     (mxlkpp) then
! ---   initialize kpp mixing
        call inikpp
      elseif (mxlmy) then
! ---   initialize m-y 2.5 mixing
        call inimy
      elseif (mxlgiss) then
! ---   initialize nasa giss mixing
        call inigiss
      endif
!
      if (linit) then
!
! ---   set up initial conditions
!
        nstep0=nstep1
        dtime0=(nstep0/nts_day)+mod(nstep0,nts_day)*(baclin/86400.0d0)
        time0=dtime0
        delt1=baclin
        if     (clmflg.eq.12) then
          mnth=   1+nint(mod(dtime0+dyear0,dyear)/dmonth)
        elseif (clmflg.eq.6) then
          mnth=2*(1+nint(mod(dtime0+dyear0,dyear)/dbimon))-1
        endif
        call inicon(mnth)
        trcrin = .false.
        call initrc(mnth)
!
! ---   setup parameters defining tidal body forces
!
        if     (tidflg.gt.0) then
          time_8=dtime0      !'baroclinic' time for body force tides
          call tides_set(0)
        endif
!
! ---   output to archive file
!
        m=mod(nstep0  ,2)+1
        n=mod(nstep0+1,2)+1
        nstep=nstep0
        time=dtime0
        call forday(dtime0,yrflag, iyear,jday,ihour)
!
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
                tmix(i,j) = temp(i,j,1,n)
                smix(i,j) = saln(i,j,1,n)
               thmix(i,j) = th3d(i,j,1,n)
              surflx(i,j) = 0.0
              mixflx(i,j) = 0.0
              buoflx(i,j) = 0.0
              bhtflx(i,j) = 0.0
              salflx(i,j) = 0.0
              wtrflx(i,j) = 0.0
            endif !ip
          enddo !i
          if (isopyc .or. mxlkrt .or. mxl_no) then
            do i=1,ii
              if (SEA_U) then
                umix(i,j)=u(i,j,1,n)
              endif !iu
              if (SEA_V) then
                vmix(i,j)=v(i,j,1,n)
              endif !iv
            enddo !i
          endif !isopyc.or.mxlkrt
        enddo !j
!$OMP   END PARALLEL DO
!
        if (mnproc.eq.1) then
        write (intvl,'(i3.3)') 0
        endif !1st tile
        if     (rstrfq.ne.0.0) then  !don't write if benchmarking (no restart)
          call archiv(n, kk, iyear,jday,ihour, intvl)
        endif
!
      else
!
! ---   start from restart file
!
        restart_cpl = .false.
        if (present(pointer_filename)) then
          open(1,file=trim(pointer_filename),form='formatted', &
                 status='old')
          read(1,'(a)') flnmrsi
          close(1)
          restart_cpl = .true.
        endif
        flnmra = trim(flnmrsi)//'.a'
        flnmrb = trim(flnmrsi)//'.b'
        call restart_in(nstep0,dtime0, flnmra,flnmrb,restart_cpl)
        surflx(:,:) = 0.0
        salflx(:,:) = 0.0
        wtrflx(:,:) = 0.0
        nstep0=nint(dtime0*(86400.0d0/baclin))
        dtime0=(nstep0/nts_day)+mod(nstep0,nts_day)*(baclin/86400.0d0)
        time0=dtime0
        delt1=baclin+baclin
        if (mnproc.eq.1) then
        write (lp,'(a,f8.1,a,i9,a, a,f8.1,a,i9,a)') &
          'restart on day',time0,' (step',nstep0,')', &
          ', wanted day',   day1,' (step',nstep1,')'
        endif !1st tile
        if (nstep0.ne.nstep1) then
          if (mnproc.eq.1) then
          write(lp,'(/a/a,f8.1/a,f8.1/)') &
            'error in hycom - wrong restart (or limits) file', &
            'restart file day is ',time0, &
            'limits start day is ',day1
          endif !1st tile
          call xcstop('(hycom)')
                 stop '(hycom)'   !won't get here
        endif !nstep0.ne.nstep1
!
        if     (clmflg.eq.12) then
          mnth=   1+nint(mod(dtime0+dyear0,dyear)/dmonth)
        elseif (clmflg.eq.6) then
          mnth=2*(1+nint(mod(dtime0+dyear0,dyear)/dbimon))-1
        endif
        call initrc(mnth)
!
! ---   setup parameters defining tidal body forces
!
        if     (tidflg.gt.0) then
          time_8=dtime0      !'baroclinic' time for body force tides
          call tides_set(0)
        endif
!
        if     (trcout .and. .not.trcrin) then
!
! ---     new tracers, so output to archive file
!
          m=mod(nstep0  ,2)+1
          n=mod(nstep0+1,2)+1
          l0=1
          l1=2
          l2=3
          l3=4
          w0=0.0
          w1=0.0
          w2=0.0
          w3=0.0
          call momtum_hs(n,m)  !calculate srfhgt
          nstep=nstep0
          time=dtime0
          call forday(dtime0,yrflag, iyear,jday,ihour)
          if (mnproc.eq.1) then
          write (intvl,'(i3.3)') 0
          endif !1st tile
          surflx(:,:) = 0.0
          salflx(:,:) = 0.0
          wtrflx(:,:) = 0.0
          if     (difout) then
            vcty(:,:,:) = 0.0
            dift(:,:,:) = 0.0
            difs(:,:,:) = 0.0
          endif
          call archiv(n, kk, iyear,jday,ihour, intvl)
        endif !archive output
      endif !initial conditions
!
! --- set barotp.pot.vort. and layer thickness (incl.bottom pressure) at
! --- u,v points
!
      call dpthuv
!
      call xctilr(dp(    1-nbdy,1-nbdy,1,1),1,2*kk, nbdy,nbdy, halo_ps)
      call xctilr(dpmixl(1-nbdy,1-nbdy,  1),1,   2, nbdy,nbdy, halo_ps)
      call xctilr(thkk(  1-nbdy,1-nbdy,1  ),1,   2, nbdy,nbdy, halo_ps)
      call xctilr(psikk( 1-nbdy,1-nbdy,1  ),1,   2, nbdy,nbdy, halo_ps)
!
      margin = nbdy
!
      nstep = nstep0+1  !for pipe_compare
      do nm=1,2
!
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          if     (nm.eq.mod(nstep+1,2)+1) then
            do i=1-margin,ii+margin
              if (SEA_P) then
                dpbl( i,j)=dpmixl(i,j,nm)
                dpbbl(i,j)=thkbot*onem
              endif !ip
            enddo !i
          endif !nm
          do i=1-margin,ii+margin
            if (SEA_P) then
              p(i,j,1)=0.0
              do k=1,kk
                p(i,j,k+1)=p(i,j,k)+dp(i,j,k,nm)
              enddo !k
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        call dpudpv(dpu(1-nbdy,1-nbdy,1,nm), &
                    dpv(1-nbdy,1-nbdy,1,nm), &
                    p,depthu,depthv, margin,max(0,margin-1))
!
        if (.false.) then
!
! ---     ISOPYC TO HYBRID RESTART ONLY
          nstep=nstep1
          n=nm
          m=mod(n,2)+1
          call hybgen(m,n, .false.)
          call xctilr(dp(1-nbdy,1-nbdy,1,n),1,kk, nbdy,nbdy, halo_ps)
          margin = nbdy
!
!$OMP     PARALLEL DO PRIVATE(j,k,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
!DIR$       PREFERVECTOR
            do i=1-margin,ii+margin
              if (SEA_P) then
                p(i,j,1)=0.0
                do k=1,kk
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                enddo !k
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                      dpv(1-nbdy,1-nbdy,1,n), &
                      p,depthu,depthv, margin,max(0,margin-1))
             dpo = dp !for initial pipe_comparall
             call pipe_comparall(m,n, 'hybgen, step')
        endif  !isopyc to hybrid restart only
!
      enddo  !nm=1,2
      nstep = nstep-1  !restore
!
! --- surface archive output flags initiialization
!
      call archiv_init
!
! --- multi-location profile initialization.
!
      if     (proffq.ne.0.0) then
        call archiv_prof_init
      endif
!

      nod=14
      nstep=nstep1
      if (mnproc.eq.1) then
      write (lp,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day', &
         time0,', goes to day',time0+day2-day1,'   (steps',nstep1, &
         ' --',nstep2,')'
      open (unit=nod,file=trim(flnminp)//'summary_out',status='unknown')
      write(nod,'(/2(a,f8.1),2(a,i9),a/)') 'model starts at day', &
         time0,', goes to day',time0+day2-day1,'   (steps',nstep1, &
         ' --',nstep2,')'
      endif !1st tile
!
      timav=time0
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
!
           dpo = dp !for initial pipe_comparall
           call pipe_comparall(m,n, 'restrt, step')
!
      if (synflt) then
! ---   initialize synthetic floats/moorings
        call floats_init(m,n,time0)
        margin = nbdy
      endif
!
      if (yrflag.lt.2) then
!
! ---   read in forcing fields for 4 consecutive months
        ma1=1.+mod(dtime0+dyear0,dyear)/dmonth
        ma0=mod(ma1+10,12)+1
        ma2=mod(ma1,   12)+1
        ma3=mod(ma2,   12)+1
        l0=1
        l1=2
        l2=3
        l3=4
        call rdforf(ma0,l0)
        call rdforf(ma1,l1)
        call rdforf(ma2,l2)
        call rdforf(ma3,l3)
      else
!
! ---   initial day of high frequency atmospheric forcing.
! ---   only two fields are used (linear interpolation in time).
        l0=1
        l1=2
        l2=3
        l3=4
        if     (windf) then
          w0=-99.9
          w1=-99.0
          w2=0.0
          w3=0.0
!!Alex no file flux on restart
#if defined (DMI_CICE_COUPLED)
        call forfunh(dtime0)
        if (mnproc.eq.1) print*,'forfunh(dtime0) HYCOMCICE'
#elif defined (USE_NUOPC_CESMBETA)
          if (.not. (cpl_swflx  .and. cpl_lwmdnflx .and. cpl_lwmupflx &
              .and.  cpl_taux   .and. cpl_tauy     .and. cpl_precip)  &
              .and.  linit ) &
              then
                  call forfunh(dtime0)
                 if (mnproc.eq.1) print*,'forfunh(dtime0) HYCOM_IN_CESM'
          else
!!Alex add initialization of jerlv0 here when forfunh is not called
! ---   jerlov used in call to swfrac_ij
            if (jerlv0.gt.0) then
! ---       calculate jerlov water type,
! ---       which governs the penetration depth of shortwave radiation.
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
              do j=1-nbdy,jj+nbdy
                do i=1-nbdy,ii+nbdy
! ---           map shallow depths to high jerlov numbers
                  jerlov(i,j)=6-max(1,min(5,int(depths(i,j)/15.0)))
                  jerlov(i,j)=max(jerlv0,jerlov(i,j))
                enddo
              enddo
!$OMP END PARALLEL DO
            else
! ---     jerlv0= 0 uses an input annual/monthly kpar field
! ---     jerlv0=-1 uses an input annual/monthly chl  field
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
              do j=1-nbdy,jj+nbdy
                do i=1-nbdy,ii+nbdy
                  jerlov(i,j)=jerlv0
                enddo
              enddo
!$OMP END PARALLEL DO
            endif ! jerlv0.le.0
            if (mnproc.eq.1) print*,'coupler forcing HYCOM_IN_CESM'
          endif ! no coupled fields and restart

#if defined (ARCTIC)
          ltripolar=.true.
#else
          ltripolar=.false.
#endif /* ARCTIC:else */
#else
          call forfunh(dtime0)
#endif /* USE_NUOPC_CESMBETA:else */

        elseif (mslprf) then
          w0=-99.9
          w1=-99.0
          w2=0.0
          w3=0.0
          call forfunhz
          call forfunhp(dtime0)
        else
          w0=0.0
          w1=0.0
          w2=0.0
          w3=0.0
          call forfunhz
        endif
      endif
!
#if defined(STOKES)
!
! --- set up fields for Stokes Drift Velocities
!       (set to zero if stdflg==0)
! --- note that stokes_set calls stokes_forfun if necessary
!
      call stokes_set(dtime0)
#endif
!
      if (jerlv0.le.0) then
! ---   read in kpar or chk field for 4 consecutive months
        mk1=1.+mod(dtime0+dyear0,dyear)/dmonth
        mk0=mod(mk1+10,12)+1
        mk2=mod(mk1,   12)+1
        mk3=mod(mk2,   12)+1
        lk0=1
        lk1=2
        lk2=3
        lk3=4
        call rdkpar(mk0,lk0)
        call rdkpar(mk1,lk1)
        call rdkpar(mk2,lk2)
        call rdkpar(mk3,lk3)
      endif
!
      if (priver) then
! ---   read in rivers field for 4 consecutive months
        mr1=1.+mod(dtime0+dyear0,dyear)/dmonth
        mr0=mod(mr1+10,12)+1
        mr2=mod(mr1,   12)+1
        mr3=mod(mr2,   12)+1
        lr0=1
        lr1=2
        lr2=3
        lr3=4
#if defined (USE_NUOPC_CESMBETA)  &&  !defined (DMI_CICE_COUPLED)
         if (.not. (cpl_orivers .and. cpl_irivers) &
             .and.  linit ) then
           call rdrivr(mr0,lr0)
           call rdrivr(mr1,lr1)
           call rdrivr(mr2,lr2)
           call rdrivr(mr3,lr3)
         endif
#else
        call rdrivr(mr0,lr0)
        call rdrivr(mr1,lr1)
        call rdrivr(mr2,lr2)
        call rdrivr(mr3,lr3)
#endif /* USE_NUOPC_CESMBETA:else */
      endif
!
      if     (clmflg.eq.12) then
! ---   read in relaxation climatology fields for 4 consecutive months
        mc1=1.+mod(dtime0+dyear0,dyear)/dmonth
        mc0=mod(mc1+10,12)+1
        mc2=mod(mc1,   12)+1
        mc3=mod(mc2,   12)+1
        lc0=1
        lc1=2
        lc2=3
        lc3=4
        call rdrlax(mc0,lc0)
        call rdrlax(mc1,lc1)
        call rdrlax(mc2,lc2)
        call rdrlax(mc3,lc3)
      elseif (clmflg.eq.6) then
! ---   read in relaxation fields for 4 consecutive bi-months
        mc1=1.+mod(dtime0+dyear0,dyear)/dbimon
        mc0=mod(mc1+4,6)+1
        mc2=mod(mc1,  6)+1
        mc3=mod(mc2,  6)+1
        lc0=1
        lc1=2
        lc2=3
        lc3=4
        call rdrlax(2*mc0-1,lc0)
        call rdrlax(2*mc1-1,lc1)
        call rdrlax(2*mc2-1,lc2)
        call rdrlax(2*mc3-1,lc3)
      else
        if (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error in hycom - unsupported clmflg value'
        call flush(lp)
        endif !1st tile
        call xcstop('(hycom)')
               stop '(hycom)'   !won't get here
      endif
!
      if     (bnstfq.ne.0.0) then  ! initialize barotropic boundary input
        wb0=-99.0
        wb1=-99.0
        call rdbaro(dtime0)
      endif
!
      if     (nestfq.ne.0.0) then  ! initialise 3-d nesting input
        wn0=-99.0
        wn1=-99.0
        call rdnest(dtime0)
      endif
!
! --- initialize incremental update.
!
      if (incflg.ne.0) then
        call incupd_init(dtime0)
!
        if (incstp.eq.1) then
              call xctmr0(54)
           call incupd(1,restrt)
           call incupd(2,restrt)
              call xctmr1(54)
        endif ! full insertion of update
      endif
!
#if defined(USE_ESMF4)
!
! --- Fill the Export State with the initial fields.
!
      m=mod(nstep0  ,2)+1
      n=mod(nstep0+1,2)+1
      call momtum_hs(n,m)
      call Export_ESMF
!
! --- Report
      call ESMF_LogWrite("HYCOM initialize routine returned", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
#else
!
! --- Only here for compatibility with coupled runs.
!
      m=mod(nstep0  ,2)+1
      n=mod(nstep0+1,2)+1
      call momtum_hs(n,m)
#endif
      call momtum_init()

#if defined (USE_NUOPC_CESMBETA)
! --- Initialization of CICE export
!!Alex update halo for calculation of  of seas surface slope for CICE (NUOPC)
      call xctilr(srfhgt(1-nbdy,1-nbdy    ),1,  1,  6,6, halo_ps)
      call xctilr(u(     1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_uv) !! for export
      call xctilr(v(     1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_vv) !! for export
      call xctilr(ubavg (1-nbdy,1-nbdy,1  ),1,   2, 6,6, halo_uv) !! for export
      call xctilr(vbavg (1-nbdy,1-nbdy,1  ),1,   2, 6,6, halo_vv) !! for export

! --- precipitation factor for CESM platform
      pcp_fact = 1.0 ! always 1.: no precipiation adjustment

!!Alex  calculation of seas surface slope for export to CICE (NUOPC)
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
            ssh_e = 0.0
            ssh_w = 0.0
            ssh_n = 0.0
            ssh_s = 0.0
            dhdx  = 0.0
            dhdy  = 0.0
            maskw = 0.0
            maske = 0.0
            maskn = 0.0
            masks = 0.0
            ! calculate the eastward slope
            if (ip(i-1,j).ne.0) then
              ssh_w = (srfhgt(i  ,j) - srfhgt(i-1,j))/(g*scux(i  ,j))
              maskw = 1.0
            endif
            if (ip(i+1,j).ne.0) then
              ssh_e = (srfhgt(i+1,j) - srfhgt(i  ,j))/(g*scux(i+1,j))
              maske = 1.0
            endif

            if ( maskw .eq.1.0 .or. maske .eq. 1.0 ) then
                dhdx=(ssh_e+ssh_w)/(maskw+maske) !! on the p-grid
            endif

            ! calculate the northward slope
           if (ip(i,j-1).ne.0) then
               ssh_s = (srfhgt(i,j  ) - srfhgt(i,j-1))/(g*scvy(i,j  ))
               masks = 1.0
            endif
            if (ip(i,j+1).ne.0) then
               ssh_n = (srfhgt(i,j+1) - srfhgt(i,j  ))/(g*scvy(i,j+1))
               maskn = 1.0
            endif
            if (masks .eq. 1.0 .or. maskn .eq. 1.0) then
               dhdy=(ssh_n+ssh_s)/(maskn+masks) !! on the p-grid
            endif
             ! convert to eastward/northward grid
               dhde(i,j) = dhdx*cos(pang(i,j)) + dhdy*sin(-pang(i,j))
               dhdn(i,j) = dhdy*cos(pang(i,j)) - dhdx*sin(-pang(i,j))
          endif
        enddo
      enddo
!$OMP END PARALLEL DO
!!Alex calculation of u and v surf for export to CICE
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
! ---       average currents over top thkcdw meters
            thksur = onem*min( thkcdw, depths(i,j) )
            usur1  = 0.0
            vsur1  = 0.0
            psur1  = 0.0
            usur2  = 0.0
            vsur2  = 0.0
            psur2  = 0.0
            do k= 1,kk
              dp1   = min( dp(i,j,k,1), max( 0.0, thksur-psur1 ) )
              usur1 = usur1 + dp1*(u(i,j,k,1)+u(i+1,j,k,1))
              vsur1 = vsur1 + dp1*(v(i,j,k,1)+v(i,j+1,k,1))
              psur1 = psur1 + dp1

              dp2   = min( dp(i,j,k,2), max( 0.0, thksur-psur2 ) )
              usur2 = usur2 + dp2*(u(i,j,k,2)+u(i+1,j,k,2))
              vsur2 = vsur2 + dp2*(v(i,j,k,2)+v(i,j+1,k,2))
              psur2 = psur2 + dp2

              if     (min(psur1,psur2).ge.thksur) then
                exit
              endif
            enddo
            utot  = 0.25*( usur1/psur1 + ubavg(i,  j,1) + &
                                         ubavg(i+1,j,1) + &
                           usur2/psur2 + ubavg(i,  j,2) + &
                                         ubavg(i+1,j,2)  )
            vtot  = 0.25*( vsur1/psur1 + vbavg(i,j,  1) + &
                                         vbavg(i,j+1,1) + &
                           vsur2/psur2 + vbavg(i,j,  2) + &
                                         vbavg(i,j+1,2)  )

             uml(i,j)=utot*cos(pang(i,j)) + vtot*sin(-pang(i,j))
             vml(i,j)=vtot*cos(pang(i,j)) - utot*sin(-pang(i,j))
          endif
        enddo
      enddo
!$OMP END PARALLEL DO
      if (linit) then
!!Alex initialization of time averaged export fields
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
          do j=1-nbdy,jj+nbdy
            do i=1-nbdy,ii+nbdy
              if (SEA_P) then
                  tml(i,j) = 0.5*(temp(i,j,1,1)+temp(i,j,1,2))
                  sml(i,j) = 0.5*(saln(i,j,1,1)+saln(i,j,1,2))
                  sshm(i,j) = srfhgt(i,j)
              endif
            enddo
          enddo
!$OMP END PARALLEL DO
      else
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
          do j=1-nbdy,jj+nbdy
            do i=1-nbdy,ii+nbdy
              if (SEA_P) then
                  if (sml(i,j).eq.0.) then ! tml and sml not initialized in the restart
                      tml(i,j) = 0.5*(temp(i,j,1,1)+temp(i,j,1,2))
                      sml(i,j) = 0.5*(saln(i,j,1,1)+saln(i,j,1,2))
                      sshm(i,j) = srfhgt(i,j)
                  endif
              endif
            enddo
          enddo
!$OMP END PARALLEL DO
      endif
      ntavg = 0 ! initialization of the counter for time averaging
#endif /* USE_NUOPC_CESMBETA */

!
! --- mean archive initialization.
!
      if     (meanfq.ne.0.0) then
        call mean_allocate
        if (mnproc.eq.1) then
          call mem_stat_print('mean_allocate:')
        endif !1st tile
        nstep=nstep0
        time =dtime0
        call mean_zero((time))
        call momtum_hs(n,m)  !calculate srfhgt
        call mean_add(n, 0.5)
      endif
!
! --- report initialization time.
! 
      call xctmrp
      call xctmr0(78)  !time since HYCOM_Init
!
      end subroutine HYCOM_Init


#if defined(USE_ESMF4)
      subroutine HYCOM_Run  &
                 (gridComp, impState, expState, extClock, rc)
!
! --- Calling parameters
      type(ESMF_GridComp)  :: gridComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
#else 
      subroutine HYCOM_Run  &
                 (endtime,pointer_filename,restart_write)
!
! --- Calling parameters
      real(8),         intent(in), optional :: endtime
      character*80, intent(in), optional :: pointer_filename
      logical,      intent(in), optional :: restart_write
#endif
#if defined (USE_NUOPC_CESMBETA)
      real :: ssh_n,ssh_s,ssh_e,ssh_w,dhdx,dhdy
      real :: maskn,masks,maske,maskw
      real :: dp1,usur1,vsur1,psur1,dp2,usur2,vsur2,psur2,thksur, &
              utot,vtot
      real :: inv_cplifq
      integer :: ld
#endif
      logical :: restart_cpl
!
! --- -------------------------
! --- execute a single timestep
! --- -------------------------
!
      logical      lfatal
      integer      i,j,k,ktr,nm,margin
      character*80 flnmra,flnmrb
!      real*8       u3max,u3min
!      integer      iumax3,jumax3,iumin3,jumin3
!
      logical hycom_isnaninf  !function to detect NaN and Inf
!
# include "stmt_fns.h"
!
#if defined(USE_ESMF4) || defined (ESPC_COUPLE)
      if     (nstep.eq.nstep0) then
        call xctmr1(78)   !time since HYCOM_Init
      else
        call xctmr1(79)   !time outside HYCOM_Run
      endif
      call xctmr0(80)
!
#endif
! --- initialize end_of_run
      end_of_run      = .false.   ! initialize on every entry
      end_of_run_cpl  = .false.   ! initialize on every entry

! --- letter 'm' refers to mid-time level (example: dp(i,j,k,m) )
! --- letter 'n' refers to old and new time level
!

#if defined (ESPC_COUPLE)
      nstep2_cpl=nint(endtime*(86400.0d0/baclin))

#  if defined (ESPC_IMPLICIT)
      w0 = (nstep2_cpl-nstep)/(nstep2_cpl-nstep1_cpl)
      w1 = 1.0 - w0
      w2 = 0
      w3 = 0
#  else
      w0=1.
      w1=0.
      w2=0.
      w3=0.
#  endif
      if(mnproc.eq.1) write(*,910) w0,w1,nstep,nstep1_cpl,nstep2_cpl
 910  format('HYCOM_Run, start...w0,w1,nstep..',2F6.3,2x,I15,2F15.1)
#endif


      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
!
      nstep=nstep+1
      dtime=(nstep/nts_day)+mod(nstep,nts_day)*(baclin/86400.0d0)
!
      time  =dtime
      time_8=dtime  !'baroclinic' time for body force tides
      if     (tidflg.gt.0 .and. &
              mod(dtime+dsmall,hours1).lt.dsmall2) then
        call tides_detide(n, .true.)  !update 49-hour filter
      endif
      hisurf=mod(dtime+dsmall,ddsurf).lt.dsmall2
      histry=mod(dtime+dsmall,ddiagf).lt.dsmall2 .or. &
              (nstep.ge.nstep2 .and. arcend)
      hiprof=mod(dtime+dsmall,dproff).lt.dsmall2
      hitile=mod(dtime+dsmall,dtilef).lt.dsmall2
      histmn=mod(dtime+dsmall,dmeanf).lt.dsmall2 .or. &
               nstep.ge.nstep2
      if     (rstrfq.eq.0.0) then  ! no restart
        restrt=.false.  ! for benchmark cases only
      elseif (drstrf.gt.dtime0) then  !at most one restart within the run
        if (rstrfq.lt.0.0) then  ! no restart at end of run
          restrt=mod(dtime+dsmall,drstrf).lt.dsmall2
        else
          restrt=mod(dtime+dsmall,drstrf).lt.dsmall2 .or. &
                 nstep.ge.nstep2
        endif
      else
        if (rstrfq.lt.0.0) then  ! no restart at end of run
          restrt=mod(dtime-dtime0+dsmall,drstrf).lt.dsmall2
        else
          restrt=mod(dtime-dtime0+dsmall,drstrf).lt.dsmall2 .or. &
                 nstep.ge.nstep2
        endif
      endif
      diagno=mod(dtime+dsmall,ddiagf).lt.dsmall2 .or. &
             restrt .or. nstep.ge.nstep2
      if (yrflag.lt.2) then
!
! ---   set weights for quasi-hermite time interpolation for
! ---   monthly atmospheric forcing fields
        x=1.+mod(dtime+dyear0,dyear)/dmonth
! ---   keep quadruplet of forcing functions centered on model time
        if (int(x).ne.ma1) then
          ma1=x
          ma0=mod(ma1+10,12)+1
          ma2=mod(ma1,   12)+1
          ma3=mod(ma2,   12)+1
          lt=l0
          l0=l1
          l1=l2
          l2=l3
          l3=lt
! ---     newest set of forcing functions overwrites set no longer needed
          call rdforf(ma3,l3)
        endif
        x=mod(x,1.)
        x1=1.-x
        w1=x1*(1.+x *(1.-1.5*x ))
        w2=x *(1.+x1*(1.-1.5*x1))
        w0=-.5*x *x1*x1
        w3=-.5*x1*x *x
!diag   if (mnproc.eq.1) then
!diag   write (lp,'(i9,'' atmos time interpolation: months'',4i3, &
!diag     '',  weights '',4f6.3)') nstep,l0,l1,l2,l3,w0,w1,w2,w3
!diag   endif !1st tile
      elseif (windf) then
!
! ---   set weights and fields for high frequency atmospheric forcing.
! ---   only two fields are used (linear interpolation in time).
#if defined (USE_NUOPC_CESMBETA)
          if (.not. (cpl_swflx  .and. cpl_lwmdnflx .and. cpl_lwmupflx &
              .and.  cpl_taux   .and. cpl_tauy     .and. cpl_precip)) &
              then
               call forfunh(dtime)
               if (mnproc.eq.1) print*,'not cpl_forcing, forfunh(dtime)'  
          endif
          if (mnproc.eq.1) print*,'cpl_forcing from coupler'
#  if defined(ARCTIC)
! --- update last active row of array
      do ld = 1, 2
        call xctila(   imp_taux,1,ld,halo_pv) ! Ocean+Ice zonal      stress
        call xctila(   imp_tauy,1,ld,halo_pv) ! Ocean+Ice meridional stress
        call xctila(  imp_swflx,1,ld,halo_ps) ! Ocean+Ice surf. solar heat    flux
        call xctila( imp_lwdflx,1,ld,halo_ps) ! Ocean+Ice surf. longwave down flux
        call xctila( imp_lwuflx,1,ld,halo_ps) ! Ocean+Ice surf. longwave up   flux
        call xctila( imp_latflx,1,ld,halo_ps) ! Ocean+Ice surf. latent        flux
        call xctila(imp_sensflx,1,ld,halo_ps) ! Ocean+Ice surf. sensible      flux
        call xctila( imp_precip,1,ld,halo_ps) ! Ocean+Ice surf. precipitation
        call xctila(imp_irivers,1,ld,halo_ps) ! surf. ice   river
        call xctila(imp_orivers,1,ld,halo_ps) ! surf. ocean river
      enddo

      if     (iceflg.ge.2 .and. icmflg.ne.3) then
        call xctila(covice,1,1,halo_ps)  !Sea Ice Concentration
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(fswice,1,1,halo_ps)  !Solar Heat Flux thru Ice to Ocean
        call xctila(flxice,1,1,halo_ps)  !Ice Freezing/Melting Heat Flux
        call xctila(sflice,1,1,halo_ps)  !Ice Virtual Salt Flux
        call xctila(wflice,1,1,halo_ps)  !Ice Freshwater   Flux
        call xctila(temice,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(thkice,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      elseif (iceflg.ge.2 .and. icmflg.eq.3) then
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      endif
#  endif /* ARCTIC */

#  else
          call forfunh(dtime)
#endif /* USE_NUOPC_CESMBETA:else */
      elseif (mslprf) then
! ---   pressure can be the only atmospheric forcing,
! ---   set weights and fields for high frequency pressure forcing.
        call forfunhp(dtime)
      endif
!
! --- set weights for quasi-hermite time interpolation for kpar.
      if     (jerlv0.le.0) then
! ---   monthly fields.
        x=1.+mod(dtime+dyear0,dyear)/dmonth
        if (int(x).ne.mk1) then
          mk1=x
          mk0=mod(mk1+10,12)+1
          mk2=mod(mk1,   12)+1
          mk3=mod(mk2,   12)+1
          lt =lk0
          lk0=lk1
          lk1=lk2
          lk2=lk3
          lk3=lt
          call rdkpar(mk3,lk3)
        endif
        x=mod(x,1.)
        x1=1.-x
        wk1=x1*(1.+x *(1.-1.5*x ))
        wk2=x *(1.+x1*(1.-1.5*x1))
        wk0=-.5*x *x1*x1
        wk3=-.5*x1*x *x
      endif
!
! --- set weights for quasi-hermite time interpolation for rivers.
      if     (priver) then
         if (.not. (cpl_orivers .and. cpl_irivers)) then
! ---   monthly fields.
            x=1.+mod(dtime+dyear0,dyear)/dmonth
            if (int(x).ne.mr1) then
               mr1=x
               mr0=mod(mr1+10,12)+1
               mr2=mod(mr1,   12)+1
               mr3=mod(mr2,   12)+1
               lt =lr0
               lr0=lr1
               lr1=lr2
               lr2=lr3
               lr3=lt
               call rdrivr(mr3,lr3)
            endif
            x=mod(x,1.)
            x1=1.-x
            wr1=x1*(1.+x *(1.-1.5*x ))
            wr2=x *(1.+x1*(1.-1.5*x1))
            wr0=-.5*x *x1*x1
            wr3=-.5*x1*x *x
         endif
      endif
!
! --- set weights for quasi-hermite time interpolation for temperature,
! --- salinity and pressure relaxation fields.
      if     (clmflg.eq.12) then
! ---   monthly fields.
        x=1.+mod(dtime+dyear0,dyear)/dmonth
        if (int(x).ne.mc1) then
          mc1=x
          mc0=mod(mc1+10,12)+1
          mc2=mod(mc1,   12)+1
          mc3=mod(mc2,   12)+1
          lt =lc0
          lc0=lc1
          lc1=lc2
          lc2=lc3
          lc3=lt
          call rdrlax(mc3,lc3)
        endif
        x=mod(x,1.)
        x1=1.-x
        wc1=x1*(1.+x *(1.-1.5*x ))
        wc2=x *(1.+x1*(1.-1.5*x1))
        wc0=-.5*x *x1*x1
        wc3=-.5*x1*x *x
      elseif (clmflg.eq.6) then
! ---   bi-monthly fields.
        x=1.+mod(dtime+dyear0,dyear)/dbimon
        if (int(x).ne.mc1) then
          mc1=x
          mc0=mod(mc1+4,6)+1
          mc2=mod(mc1,  6)+1
          mc3=mod(mc2,  6)+1
          lt =lc0
          lc0=lc1
          lc1=lc2
          lc2=lc3
          lc3=lt
          call rdrlax(2*mc3-1,lc3)
        endif
        x=mod(x,1.)
        x1=1.-x
        wc1=x1*(1.+x *(1.-1.5*x ))
        wc2=x *(1.+x1*(1.-1.5*x1))
        wc0=-.5*x *x1*x1
        wc3=-.5*x1*x *x
      endif
!diag if (mnproc.eq.1) then
!diag write (lp,'(i9,'' relax time interpolation: months'',4i3, &
!diag   '',  weights '',4f6.3)') nstep,lc0,lc1,lc2,lc3,wc0,wc1,wc2,wc3
!diag endif !1st tile
!
#if defined(USE_ESMF4)
      if     (get_import) then     ! new ESMF Import fields
        call Import_ESMF
        if (incflg.ne.0) then
          call incupd_si_c(dtime)
        endif
        if     (nstep.eq.nstep0+1) then
          call momtum_hs(n,m)
          time=dtime0
          call forday(dtime0,yrflag, iyear,jday,ihour)
          call Archive_ESMF(iyear,jday,ihour)
          time=dtime
        endif !initial ESMF Archive
      endif !get_import
#endif
!
      if     (bnstfq.ne.0.0) then  ! new fields for baro nesting
        call rdbaro(dtime)
      endif
!
      if     (nestfq.ne.0.0) then  ! new fields for 3-d  nesting
        call rdnest(dtime)
      endif
!
         call pipe_comparall(n,m, 'ENTERm, step')
         call pipe_comparall(m,n, 'ENTERn, step')
         call xctmr0(55)
      call asselin_save(m,n)
         call xctmr1(55)
         call xctmr0(40)
      call cnuity(m,n)
         call xctmr1(40)
         call pipe_comparall(m,n, 'cnuity, step')
         call xctmr0(42)
      if     (momtyp.eq.2) then
        call momtum(m,n)
      else   !momtyp==3,4
        call momtum4(m,n)
      endif
         call xctmr1(42)
         call pipe_comparall(m,n, 'momtum, step')
         call xctmr0(43)
      call barotp(m,n)
         call xctmr1(43)
         call pipe_comparall(m,n, 'barotp, step')
         call xctmr0(41)
      call tsadvc(m,n)
         call xctmr1(41)
         call pipe_comparall(m,n, 'tsadvc, step')
         call xctmr0(44)
      call thermf(m,n, dtime)
         call xctmr1(44)
         call pipe_comparall(m,n, 'thermf, step')
      if (icegln) then
            call xctmr0(45)
         call icloan(m,n)
         call thermf_oi(m,n)
            call xctmr1(45)
            call pipe_comparall(m,n, 'icloan, step')
      elseif (iceflg.ge.2) then
            call xctmr0(45)
         call thermf_oi(m,n)
            call xctmr1(45)
            call pipe_comparall(m,n, 'icecpl, step')
      else
            call xctmr0(45)
         call thermf_oi(m,n)
            call xctmr1(45)
            call pipe_comparall(m,n, 'thermi, step') !thermf_oi
      endif !icegln:icecpl:else
      if (trcout) then
           call xctmr0(50)
        call trcupd(m,n)
           call xctmr1(50)
           call pipe_comparall(m,n, 'trcupd, step')
      endif !trcout
      if (hybrid) then
         diagsv = diagno
         diagno = diagno .or. nstep.eq.nstep0+1 .or. &
                  histry .or. hisurf .or. hiprof .or. hitile .or. &
                  mod(dtime+dsmall,days6).lt.dsmall2
         if (mxlkpp .or. mxlmy .or. mxlgiss) then
               call xctmr0(46)
            call mxkprf(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkprf, step')
         elseif (mxlpwp) then
               call xctmr0(46)
            call mxpwp(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxpwp,  step')
         elseif (mxlkta) then
               call xctmr0(46)
            call mxkrta(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkrta, step')
         elseif (mxlktb) then
               call xctmr0(46)
            call mxkrtb(m,n)
               call xctmr1(46)
               call pipe_comparall(m,n, 'mxkrtb, step')
         else
               call xctmr0(46)
               call xctmr1(46)
         endif !mixed layer
         diagno = diagsv
         if (mxlkpp .or. mxlmy .or. mxlgiss) then  !tsoff in mxkprf
               call xctmr0(47)
               call xctmr1(47)
               call xctmr0(48)
               call xctmr1(48)
         else  ! mxlpwp has dypflg=2
               call xctmr0(47)
            call convch(m,n)
               call xctmr1(47)
               call pipe_comparall(m,n, 'convch, step')
           if     (dypflg.eq.1) then  ! KPP-like, tsoff in diapf1
                 call xctmr0(48)
              call diapf1(m,n)
                 call xctmr1(48)
                 call pipe_comparall(m,n, 'diapf1, step')
           elseif (dypflg.eq.2) then  ! explicit, tsoff in diapf2
                 call xctmr0(48)
              call diapf2(m,n)
                 call xctmr1(48)
                 call pipe_comparall(m,n, 'diapf2, step')
           else
                 call xctmr0(48)
                 call xctmr1(48)
           endif
        endif !diapycnal mixing
           call xctmr0(54)
        if     (incflg.ne.0 .and. incstp.gt.1) then
        call incupd(n,restrt)
           call pipe_comparall(m,n, 'incupd, step')
        endif ! incremental update
        call stfupd(n)
           call xctmr1(54)
           call pipe_comparall(m,n, 'stfupd, step')
           call xctmr0(56)
        call asselin_filter(m,n)
           call xctmr1(56)
           call pipe_comparall(m,n, 'asseln, step')
           call xctmr0(49)
        call hybgen(m,n, hybraf)
           call xctmr1(49)
           call pipe_comparall(m,n, 'hybgen, step')
      else  ! isopyc
            call xctmr0(46)
         call mxkrtm(m,n)
            call xctmr1(46)
            call pipe_comparall(m,n, 'mxkrtm, step')
            call xctmr0(47)
         call convcm(m,n)
            call xctmr1(47)
            call pipe_comparall(m,n, 'convcm, step')
            call xctmr0(48)
         call diapf3(m,n)
            call xctmr1(48)
            call pipe_comparall(m,n, 'diapf3, step')
            call xctmr0(54)
            call xctmr1(54)
            call xctmr0(56)
         call asselin_filter(m,n)
            call xctmr1(56)
            call pipe_comparall(m,n, 'asseln, step')
            call xctmr0(49)
            call xctmr1(49)
      endif !hybrid:isopyc

#if defined (USE_NUOPC_CESMBETA)
!!Alex update halo for calculation of  of seas surface slope for CICE (NUOPC)
      call xctilr(srfhgt(1-nbdy,1-nbdy    ),1,  1,  6,6, halo_ps)
      call xctilr(u(     1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_uv) !! for export
      call xctilr(v(     1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_vv) !! for export
      call xctilr(ubavg (1-nbdy,1-nbdy,1  ),1,   2, 6,6, halo_uv) !! for export
      call xctilr(vbavg (1-nbdy,1-nbdy,1  ),1,   2, 6,6, halo_vv) !! for export
#endif /* USE_NUOPC_CESMBETA */
!
! --- update floats/moorings
      if     (synflt) then
        nstepfl=nstepfl+1
        if (nstepfl.eq.1) then
          iofl=1
          call floats(m,n,0.0,iofl)
        endif
        if (mod(nstepfl,nfldta).eq.0) then
          iofl=0
          if (mod(nstepfl,nflsam).eq.0) then
            iofl=1
          endif
          call floats(m,n,nstepfl*baclin/86400.0,iofl)
        endif
      endif !synflt
!
! ---------------------------------------------------------------------------
!
! --- output and diagnostic calculations
!
! ---------------------------------------------------------------------------
!
      lfatal    = .false.
      diag_tide = tidflg.gt.0 .and. &
                  mod(dtime+dsmall,hours1).lt.dsmall2  !at least hourly
      if (diagno .or. &
          histry .or. hiprof .or. hitile .or. hisurf .or. histmn .or. &
          nstep.eq.nstep0+1 .or. &
          diag_tide .or. &
          mod(dtime+dsmall,ddsurf).lt.dsmall2 .or. &
          mod(dtime+dsmall, days1).lt.dsmall2     ) then     ! at least daily
!
        call forday(dtime,yrflag, iyear,jday,ihour)
        write(c_ydh,'('' ('',i4.4,''/'',i3.3,1x,i2.2,'')'')') &
          iyear,jday,ihour
!
! ---   diagnose mean sea surface height
        if     (.not.allocated(sminy)) then
          allocate( sminy(1:jj), smaxy(1:jj) )
        endif
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          sminy(j)= huge(sminy(j))
          smaxy(j)=-huge(smaxy(j))
          do i=1,ii
            if (SEA_P) then
! ---         always use mass-conserving diagnostics
              oneta(i,j,n) = max( oneta0, 1.0 + pbavg(i,j,n)/pbot(i,j) )
              if     (tidflg.gt.0) then
                util2(i,j)=(srfhgt(i,j)/g)**2*scp2(i,j)
              endif
                util3(i,j)=srfhgt(i,j)*scp2(i,j)
              if     (sshflg.eq.1) then
                util4(i,j)=steric(i,j)*scp2(i,j)
              else
                util4(i,j)=montg1(i,j)*scp2(i,j)
              endif !sshflg
                util5(i,j)=pbot(i,j)*scp2(i,j)
!
              sminy(j)=min(sminy(j),srfhgt(i,j))
              smaxy(j)=max(smaxy(j),srfhgt(i,j))
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
        smin=minval(sminy(1:jj))
        smax=maxval(smaxy(1:jj))
        call xcminr(smin)
        call xcmaxr(smax)
        call xcsum( dsum, util3,ipa)
        call xcsum( dsmt, util4,ipa)
        call xcsum( dbot, util5,ipa)
        sum=dsum
        smt=dsmt
!        write(6,'("sum,dsum=",G25.17,G25.17)')sum,dsum
!        print *,'srfhgt(1,1)=',srfhgt(1,1)
!        print *,'scp2(1,1)  =',scp2(1,1)
!       write(6,'(I3,I3,G25.17)')iumax3,jumax3,u3max
!       write(6,'(I3,I3,G25.17)')iumin3,jumin3,u3min
        if (mnproc.eq.1) then
        write (lp,'(i9,a, &
                   &'' mean      SSH (mm):'',f8.2, &
                   &''  ('',1pe8.1,'' to '',e8.1,'')'')') &
        nstep,c_ydh, &
        sum/(area*svref*onemm),smin/(svref*onemm),smax/(svref*onemm)
        write(nod,'(i9,a, &
                   &'' mean      SSH (mm):'',f8.2, &
                   &''  ('',1pe8.1,'' to '',e8.1,'')'')') &
        nstep,c_ydh, &
        sum/(area*svref*onemm),smin/(svref*onemm),smax/(svref*onemm)
        if     (sshflg.ne.1) then
          write (lp,'(i9,a, &
                     &'' mean MontgPot (mm):'',f8.2)') &
          nstep,c_ydh, &
          smt/(area*svref*onemm)
          call flush(lp)
          write(nod,'(i9,a, &
                     &'' mean MontgPot (mm):'',f8.2)') &
          nstep,c_ydh, &
          smt/(area*svref*onemm)
        else
          write (lp,'(i9,a, &
                     &'' mean   St-SSH (mm):'',f8.2)') &
          nstep,c_ydh, &
          smt/(area*svref*onemm)
          call flush(lp)
          write(nod,'(i9,a, &
                     &'' mean   St-SSH (mm):'',f8.2)') &
          nstep,c_ydh, &
          smt/(area*svref*onemm)
        endif !sshflg
        call flush(nod)
        endif !1st tile
! ---   NaN detection.
        if     (hycom_isnaninf(smin) .or. &
                hycom_isnaninf(sum)  .or. &
                hycom_isnaninf(smax) .or. &
                hycom_isnaninf(smt)      ) then
          if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - NaN or Inf detected'
          write(lp,*)
          call flush(lp)
          endif !1st tile
          lfatal = .true.  !delay exit to allow archive output
        endif !NaN
!                
        if     (tidflg.gt.0) then
          call xcsum( dsms, util2,ipa)
          sms=dsms/area
!         
          call xctilr(u(    1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
          call xctilr(v(    1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
          call xctilr(ubavg(1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_uv)
          call xctilr(vbavg(1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_vv)
!         
#if defined(STOKES) 
          dskesa=0.0d0
#endif
          dskea =0.0d0
          do k=1,kk
!$OMP       PARALLEL DO PRIVATE(j,i,utotp,vtotp) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  utotp=0.5*( u(i,  j,k,n)+ubavg(i,  j,n) + &
                              u(i+1,j,k,n)+ubavg(i+1,j,n)  )
                  vtotp=0.5*( v(i,j,  k,n)+vbavg(i,j,  n) + &
                              v(i,j+1,k,n)+vbavg(i,j+1,n)  )
                  util4(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                            0.5*(rhoref+th3d(i,j,k,n)+thbase)* &
                                (utotp**2+vtotp**2)
#if defined(STOKES)
                  utotp=utotp+0.5*( usd(i,  j,k) + &
                                    usd(i+1,j,k)  )
                  vtotp=vtotp+0.5*( vsd(i,j,  k) + &
                                    vsd(i,j+1,k)  )
!
                  util6(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                            0.5*(rhoref+th3d(i,j,k,n)+thbase)* &
                                (utotp**2+vtotp**2)
#endif
                endif !ip
              enddo !i
            enddo !j
!$OMP       END PARALLEL DO
#if defined(STOKES)
            call xcsum(dskes, util6,ipa)
            dskesa=dskesa+dskes
#endif
            call xcsum(dske, util4,ipa)
            dskea=dskea+dske
          enddo !k
#if defined(STOKES)
          smx=(dskesa-dskea)/(area*onem)
#endif
          sum=dskea/(area*onem)
          if (mnproc.eq.1) then
#if defined(STOKES)
          write (lp,'(i9,a, &
                       &'' region-wide mean SKE:'',f20.10)') &
              nstep,c_ydh, &
                smx
#endif
          write (lp,'(i9,a, &
                       &'' region-wide mean KE: '',f20.10)') &
              nstep,c_ydh, &
                sum
          write (lp,'(i9,a, &
                       &'' region-wide mean APE:'',f20.10)') &
              nstep,c_ydh, &
                sms*0.5*g*(rhoref+thbase)
          endif
        endif !tidflg
!     else
!       if (mnproc.eq.1) then
!       write (lp,'('' time step ='',i9)') nstep
!       call flush(lp)
!       endif !1st tile
      endif !daily or hourly
!
! --- diagnose heat/salt flux, ice, layer thickness and temperature,
! --- mean temperature and mean kinetic energy
! --- note that mixed-layer fields must be switched on in mxkprf
      if (diagno .or. &
          nstep.eq.nstep0+1 .or. &
          mod(dtime+dsmall,days6).lt.dsmall2) then  ! at least every 6 days
!
        call forday(dtime,yrflag, iyear,jday,ihour)
        write(c_ydh,'('' ('',i4.4,''/'',i3.3,1x,i2.2,'')'')') &
          iyear,jday,ihour
!
        if (thermo .or. sstflg.gt.0 .or. srelax) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                util1(i,j)=buoflx(i,j)*scp2(i,j)
                util2(i,j)=bhtflx(i,j)*scp2(i,j)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(dsum, util1,ipa)
          call xcsum(dsmt, util2,ipa)
          sum= (dsum*1.00D9)/area  ! 1.e9*m**2/sec**3
          smt= (dsmt*1.00D9)/area  ! 1.e9*m**2/sec**3
          if (mnproc.eq.1) then
          write (lp, '(i9,a, &
             &'' mean BFL (m^2/s^3):'',f8.2, &
                            &'' hfl:'',f8.2)') &
            nstep,c_ydh, &
            sum,smt
          call flush(lp)
          write (nod,'(i9,a, &
             &'' mean BFL (m^2/s^3):'',f8.2, &
                            &'' hfl:'',f8.2)') &
            nstep,c_ydh, &
            sum,smt
          call flush(nod)
          endif !1st tile
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                util1(i,j)=surflx(i,j)*scp2(i,j)
                util2(i,j)=mixflx(i,j)*scp2(i,j)
                util3(i,j)=sstflx(i,j)*scp2(i,j)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(dsum, util1,ipa)
          call xcsum(dsmt, util2,ipa)
          call xcsum(d3,   util3,ipa)
          sum= dsum/area
          smt= dsmt/area
          smr= d3  /area
          if (mnproc.eq.1) then
          write (lp, '(i9,a, &
             &'' mean HFLUX (w/m^2):'',f8.2, &
                            &'' sst:'',f8.2, &
                            &''  ml:'',f8.2)') &
            nstep,c_ydh, &
            sum,smr,smt
          call flush(lp)
          write (nod,'(i9,a, &
             &'' mean HFLUX (w/m^2):'',f8.2, &
                            &'' sst:'',f8.2, &
                            &''  ml:'',f8.2)') &
            nstep,c_ydh, &
            sum,smr,smt
          call flush(nod)
          endif !1st tile
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                util1(i,j)=wtrflx(i,j)*scp2(i,j)
                util2(i,j)=sssflx(i,j)*scp2(i,j)/max(saln(i,j,1,n), &
                                                     epsil)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(dsms, util1,ipa)
          call xcsum(dsum, util2,ipa)
          sms= (dsms*svref*7.0D0*8.64D7)/area  ! P-E in mm/week
          smr=-(dsum*svref*7.0D0*8.64D7)/area  ! SSS relax in mm/week
          if (mnproc.eq.1) then
          write (lp, '(i9,a, &
             &'' mean WFLUX (mm/wk):'',f8.2, &
                            &'' sss:'',f8.2)') &
            nstep,c_ydh, &
            sms,smr
          call flush(lp)
          write (nod,'(i9,a, &
             &'' mean WFLUX (mm/wk):'',f8.2, &
                            &'' sss:'',f8.2)') &
            nstep,c_ydh, &
            sms,smr
          call flush(nod)
          endif !1st tile
        endif
!
        if (iceflg.ne.0) then  ! basin-wide ice
! ---     use CICE fields for ice statistics when available
          if     (iceflg.eq.1) then
            si_c(:,:) = covice(:,:)
            si_h(:,:) = thkice(:,:)
            si_t(:,:) = temice(:,:)
          endif
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                if     (si_c(i,j).ne.0.0) then
                  util2(i,j)=si_c(i,j)*          scp2(i,j)
                  util3(i,j)=          si_h(i,j)*scp2(i,j)
                  util4(i,j)=si_c(i,j)*si_t(i,j)*scp2(i,j)
                else
                  util2(i,j)=0.0
                  util3(i,j)=0.0
                  util4(i,j)=0.0
                endif
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(d2, util2,ipa)
          call xcsum(d3, util3,ipa)
          call xcsum(d4, util4,ipa)
          if     (d2.gt.0.0d0) then
            sum=d3/d2   !average ice thickness,   where there is ice
            smt=d4/d2   !average ice temperature, where there is ice
            sms=d2/area * 100.0  !ice coverage, percent of total area
          else
            sum=0.0 
            smt=0.0 
            sms=0.0 
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                   &'' mean  ice thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(lp)
          write(nod,'(i9,a, &
                   &'' mean  ice thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(nod)
          endif !1st tile
        endif  ! iceflg.ne.0
!
        if (nreg.ne.0 .and. iceflg.ne.0) then  ! southern hemisphere ice
          d2a = d2
          d3a = d3
          d4a = d4
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                if     (si_c(i,j).ne.0.0 .and. &
                        plat(i,j).lt.0.0      ) then
                  util2(i,j)=si_c(i,j)*          scp2(i,j)
                  util3(i,j)=          si_h(i,j)*scp2(i,j)
                  util4(i,j)=si_c(i,j)*si_t(i,j)*scp2(i,j)
                else
                  util2(i,j)=0.0
                  util3(i,j)=0.0
                  util4(i,j)=0.0
                endif
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(d2, util2,ipa)
          call xcsum(d3, util3,ipa)
          call xcsum(d4, util4,ipa)
          if     (d2.gt.0.0d0) then
            sum=d3/d2   !average ice thickness,   where there is ice in S.H.
            smt=d4/d2   !average ice temperature, where there is ice in S.H.
            sms=d2/area * 100.0  !S.H. ice coverage, percent of total area
          else
            sum=0.0 
            smt=0.0 
            sms=0.0 
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                   &'' mean SH I thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(lp)
          write(nod,'(i9,a, &
                   &'' mean SH I thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(nod)
          endif !1st tile
!
          d2 = d2a - d2
          d3 = d3a - d3
          d4 = d4a - d4
          if     (d2.gt.0.0d0) then
            sum=d3/d2   !average ice thickness,   where there is ice in N.H.
            smt=d4/d2   !average ice temperature, where there is ice in N.H.
            sms=d2/area * 100.0  !N.H. ice coverage, percent of total area
          else
            sum=0.0 
            smt=0.0 
            sms=0.0 
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                   &'' mean NH I thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(lp)
          write(nod,'(i9,a, &
                   &'' mean NH I thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' pcen:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
          call flush(nod)
          endif !1st tile
        endif  ! iceflg.ne.0 .and. nreg.ne.0
!
        if     (icegln .and. icmflg.eq.3) then
!
! ---     HYCOM covice when relaxing to CICE si_
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                util1(i,j) = scp2(i,j)*covice(i,j)
                if     (plat(i,j).lt.0.0) then
                  util2(i,j) = util1(i,j)
                else
                  util2(i,j) = 0.0
                endif
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(dsmt, util1,ipa)
          call xcsum(dsms, util2,ipa)
          smt=dsmt/area * 100.0
          sms=dsms/area * 100.0
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                   &'' mean   covice pcen:'',f9.3, &
                                 &''  S.H:'',f7.3, &
                                 &''  N.H:'',f7.3)') &
          nstep,c_ydh, &
          smt,sms,smt-sms
          call flush(lp)
          write(nod,'(i9,a, &
                   &'' mean   covice pcen:'',f9.3, &
                                 &''  S.H:'',f7.3, &
                                 &''  N.H:'',f7.3)') &
          nstep,c_ydh, &
          smt,sms,smt-sms
          call flush(nod)
          endif !1st tile
        endif  !HYCOM covice
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              util1(i,j)=dpmixl(i,j,n)*scp2(i,j)
              util2(i,j)=dpmixl(i,j,n)*scp2(i,j)*tmix(i,j)
              util3(i,j)=dpmixl(i,j,n)*scp2(i,j)*smix(i,j)
              util4(i,j)=temp(i,j,1,n)*scp2(i,j)
              util5(i,j)=saln(i,j,1,n)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
        call xcsum(dsum, util1,ipa)
        call xcsum(dsmt, util2,ipa)
        call xcsum(dsms, util3,ipa)
        if     (dsum.ne.0.0d0) then
          sum=dsum/(area*onem)
          smt=dsmt/dsum
          sms=dsms/dsum
        else
          sum=0.0
          smt=0.0
          sms=0.0
        endif
        if (mnproc.eq.1) then
        write (lp,'(i9,a, &
                   &'' mean mixl thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' saln:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
        call flush(lp)
        write(nod,'(i9,a, &
                   &'' mean mixl thk. (m):'',f8.2, &
                                &''  temp:'',f7.3, &
                                 &'' saln:'',f7.3)') &
            nstep,c_ydh, &
            sum,smt,sms
        call flush(nod)
        endif !1st tile
!
        call xcsum(dsmt, util4,ipa)
        call xcsum(dsms, util5,ipa)
        smt=dsmt/area
        sms=dsms/area
        if (mnproc.eq.1) then
        write (lp,'(i9,a, &
                 &'' mean surf thk. (m):'',f8.2, &
                              &''   sst:'',f7.3, &
                               &''  sss:'',f7.3)') &
          nstep,c_ydh, &
          dp00*qonem,smt,sms  !dp00 is max, not mean, surf thk.
        call flush(lp)
        write(nod,'(i9,a, &
                 &'' mean surf thk. (m):'',f8.2, &
                              &''   sst:'',f7.3, &
                               &''  sss:'',f7.3)') &
          nstep,c_ydh, &
          dp00*qonem,smt,sms  !dp00 is max, not mean, surf thk.
        call flush(nod)
        endif !1st tile
! ---   NaN detection.
        if     (hycom_isnaninf(sms) .or. &
                hycom_isnaninf(smt)      ) then
          if (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) 'error - NaN or Inf detected'
          write(lp,*)
          call flush(lp)
          endif !1st tile
          lfatal = .true.  !delay exit to allow archive output
        endif !NaN
!
        if     (relaxf .or. sstflg.ne.0) then
!
! ---     mean surface climatology.
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                if     (relaxf .and. sstflg.le.1) then
                  util1(i,j)=scp2(i,j)* &
                      (twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1 &
                      +twall(i,j,1,lc2)*wc2+twall(i,j,1,lc3)*wc3)
                elseif (natm.eq.2) then !hf synoptic observed sst
                  util1(i,j)=scp2(i,j)* &
                      (seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1)
                else !monthly synoptic observed sst
                  util1(i,j)=scp2(i,j)* &
                      (seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1 &
                      +seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3)
                endif
                if     (relaxf) then !sss
                  util2(i,j)=scp2(i,j)* &
                        (swall(i,j,1,lc0)*wc0+swall(i,j,1,lc1)*wc1 &
                        +swall(i,j,1,lc2)*wc2+swall(i,j,1,lc3)*wc3)
                elseif (natm.eq.2) then !hf synoptic observed surface temp
                  util1(i,j)=scp2(i,j)* &
                      (surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1)
                else !monthly synoptic observed surface temperature
                  util1(i,j)=scp2(i,j)* &
                      (surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1 &
                      +surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3)
                endif
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
          call xcsum(dsmt, util1,ipa)
          call xcsum(dsms, util2,ipa)
          smt=dsmt/area
          sms=dsms/area
          if (mnproc.eq.1) then
          if     (relaxf) then
            write (lp,'(i9,a, &
                     &'' mean clim thk. (m):'',f8.2, &
                                  &''   sst:'',f7.3, &
                                   &''  sss:'',f7.3)') &
              nstep,c_ydh, &
              thkmin,smt,sms
            call flush(lp)
            write(nod,'(i9,a, &
                     &'' mean clim thk. (m):'',f8.2, &
                                  &''   sst:'',f7.3, &
                                   &''  sss:'',f7.3)') &
              nstep,c_ydh, &
              thkmin,smt,sms
            call flush(nod)
          else !.not.relaxf
            write (lp,'(i9,a, &
                     &'' mean clim thk. (m):'',f8.2, &
                                  &''   sst:'',f7.3, &
                                   &'' surt:'',f7.3)') &
              nstep,c_ydh, &
              thkmin,smt,sms
            call flush(lp)
            write(nod,'(i9,a, &
                     &'' mean clim thk. (m):'',f8.2, &
                                  &''   sst:'',f7.3, &
                                   &'' surt:'',f7.3)') &
              nstep,c_ydh, &
              thkmin,smt,sms
            call flush(nod)
          endif !relaxf:else
          endif !1st tile
        endif
!
        call xctilr(u(    1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
        call xctilr(v(    1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
        call xctilr(ubavg(1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_uv)
        call xctilr(vbavg(1-nbdy,1-nbdy,  n),1, 1, 1,1, halo_vv)
!
        dsuma=0.0d0
        dsmta=0.0d0
        dsmsa=0.0d0
        dsmra=0.0d0
        dskea=0.0d0
#if defined(STOKES)
        dskesa=0.0d0
#endif
        do k=1,kk
!$OMP     PARALLEL DO PRIVATE(j,i,utotp,vtotp) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                  utotp=0.5*( u(i,  j,k,n)+ubavg(i,  j,n) + &
                              u(i+1,j,k,n)+ubavg(i+1,j,n)  )
                  vtotp=0.5*( v(i,j,  k,n)+vbavg(i,j,  n) + &
                              v(i,j+1,k,n)+vbavg(i,j+1,n)  )
                util4(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                          0.5*(rhoref+th3d(i,j,k,n)+thbase)* &
                              (utotp**2+vtotp**2)
#if defined(STOKES)
                  utotp=utotp+0.5*( usd(i,  j,k) + &
                                    usd(i+1,j,k)  )
                  vtotp=vtotp+0.5*( vsd(i,j,  k) + &
                                    vsd(i,j+1,k)  )
                  util6(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                            0.5*(rhoref+th3d(i,j,k,n)+thbase)* &
                                (utotp**2+vtotp**2)
#endif
!
                util1(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)
                util2(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                                      temp(i,j,k,n)
                util3(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                                      saln(i,j,k,n)
                util5(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                                      th3d(i,j,k,n)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
#if defined(STOKES)
          call xcsum(dskes, util6,ipa)
          dskesa=dskesa+dskes
#endif
          call xcsum(dsum, util1,ipa)
          call xcsum(dsmt, util2,ipa)
          call xcsum(dsms, util3,ipa)
          call xcsum(dsmr, util5,ipa)
          call xcsum(dske, util4,ipa)
          dsuma=dsuma+dsum
          dsmta=dsmta+dsmt
          dsmsa=dsmsa+dsms
          dsmra=dsmra+dsmr
          dskea=dskea+dske
          if     (dsum.ne.0.0d0) then
            sum=dsum/(area*onem)
            smt=dsmt/dsum
            sms=dsms/dsum
          else
            sum=0.0
            smt=0.0
            sms=0.0
          endif
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                     &'' mean L '',i2,'' thk. (m):'',f8.2, &
                                        &''  temp:'',f7.3, &
                                         &'' saln:'',f7.3)') &
              nstep,c_ydh, &
              k,sum,smt,sms
          call flush(lp)
          write(nod,'(i9,a, &
                     &'' mean L '',i2,'' thk. (m):'',f8.2, &
                                        &''  temp:'',f7.3, &
                                         &'' saln:'',f7.3)') &
              nstep,c_ydh, &
              k,sum,smt,sms
          call flush(nod)
          endif !1st tile
        enddo !k
#if defined(STOKES)
        smx=(dskesa-dskea)/(area*onem)
#endif
        sum =dskea/(area*onem)
        smt =dsmta/dsuma
        sms =dsmsa/dsuma
        smr =dsmra/dsuma
        smsa=(dsmsa - dbot*saln0)/(area*onem)  !salt anomaly
        if (mnproc.eq.1) then
#if defined(STOKES)
        write (lp,'(i9,a, &
                     &'' region-wide mean Stokes K.E.:'',f20.10)') &
            nstep,c_ydh, &
              smx
#endif
        write (lp,'(i9,a, &
                     &'' region-wide mean Kin. Energy:'',f20.10)') &
            nstep,c_ydh, &
              sum
        write (lp,'(i9,a, &
                     &'' region-wide mean Temperature:'',f20.10)') &
            nstep,c_ydh, &
              smt
        write (lp,'(i9,a, &
                     &'' region-wide mean Salinity:   '',f20.10)') &
            nstep,c_ydh, &
              sms
        write (lp,'(i9,a, &
                     &'' region-wide Salt Anomaly:    '',f20.10)') &
            nstep,c_ydh, &
              smsa
        write (lp,'(i9,a, &
                     &'' region-wide mean Density Dev:'',f20.10)') &
            nstep,c_ydh, &
              smr
        call flush(lp)
#if defined(STOKES)
        write (nod,'(i9,a, &
                     &'' region-wide mean Stokes K.E.:'',f20.10)') &
            nstep,c_ydh, &
              smx
#endif
        write(nod,'(i9,a, &
                     &'' region-wide mean Kin. Energy:'',f20.10)') &
            nstep,c_ydh, &
              sum
        write(nod,'(i9,a, &
                     &'' region-wide mean Temperature:'',f20.10)') &
            nstep,c_ydh, &
              smt
        write(nod,'(i9,a, &
                     &'' region-wide mean Salinity:   '',f20.10)') &
            nstep,c_ydh, &
              sms
        write(nod,'(i9,a, &
                     &'' region-wide Salt Anomaly:    '',f20.10)') &
            nstep,c_ydh, &
              smsa
        write(nod,'(i9,a, &
                     &'' region-wide mean Density Dev:'',f20.10)') &
            nstep,c_ydh, &
              smr
        call flush(nod)
        endif !1st tile
!
        do ktr= 1,ntracr
          dsumtr(ktr)=0.0d0
          do k=1,kk
!$OMP       PARALLEL DO PRIVATE(j,i) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  util2(i,j)=oneta(i,j,n)*dp(i,j,k,n)*scp2(i,j)* &
                                      tracer(i,j,k,n,ktr)
                endif !ip
              enddo !i
            enddo !j
!$OMP       END PARALLEL DO
            call xcsum(dsmt, util2,ipa)
            dsumtr(ktr)=dsumtr(ktr)+dsmt
          enddo !k
          smt=dsumtr(ktr)/dsuma  !dsuma still good from K.E. loops
          if (mnproc.eq.1) then
          write (lp,'(i9,a, &
                       &'' region-wide mean tracer'',i3.2, &
                       &'':  '',f20.10)') &
              nstep,c_ydh, &
                ktr,smt
          call flush(lp)
          write(nod,'(i9,a, &
                       &'' region-wide mean tracer'',i3.2, &
                       &'':  '',f20.10)') &
              nstep,c_ydh, &
                ktr,smt
          call flush(nod)
          endif !1st tile
! ---     NaN detection, for each tracer.
          if     (hycom_isnaninf(smt)) then
            if (mnproc.eq.1) then
            write(lp,*)
            write(lp,*) 'error - NaN or Inf detected'
            write(lp,*)
            call flush(lp)
            endif !1st tile
            lfatal = .true.  !delay exit to allow archive output
          endif !NaN
          if     (ktr.ge.3 .and. trcflg(ktr-2).eq.903) then !NPZ
            smt=(dsumtr(ktr)+dsumtr(ktr-1)+dsumtr(ktr-2))/dsuma
            if (mnproc.eq.1) then
            write (lp,'(i9,a, &
                         &'' region-wide mean N+P+Z:      '',f20.10)') &
                nstep,c_ydh, &
                  smt
            call flush(lp)
            write(nod,'(i9,a, &
                         &'' region-wide mean N+P+Z:      '',f20.10)') &
                nstep,c_ydh, &
                  smt
            call flush(nod)
            endif !1st tile
          elseif (ktr.ge.4 .and. trcflg(ktr-3).eq.904) then !NPZD
            smt=(dsumtr(ktr)  +dsumtr(ktr-1)+ &
                 dsumtr(ktr-2)+dsumtr(ktr-3))/dsuma
            if (mnproc.eq.1) then
            write (lp,'(i9,a, &
                         &'' region-wide mean N+P+Z+D:    '',f20.10)') &
                nstep,c_ydh, &
                  smt
            call flush(lp)
            write(nod,'(i9,a, &
                         &'' region-wide mean N+P+Z+D:    '',f20.10)') &
                nstep,c_ydh, &
                  smt
            call flush(nod)
            endif !1st tile
          endif !NPZ:NPZD
        enddo !ktr
      endif  !diagno ...
!
! --- diagnose meridional overturning and heat flux
      if     (mod(dtime+dsmall,dmonth).lt.dsmall2) then
        call xctmr0(52)
        call overtn(dtime,dyear)
        call xctmr1(52)
      elseif (nstep.ge.nstep2) then
        call xctmr0(52)
        call overtn(dtime,dyear)
        call xctmr1(52)
      endif
!
      if     (meanfq.ne.0.0) then
        if     (.not. histmn) then
          call mean_add(n, 1.0)
        else  ! histmn
          call mean_add(n, 0.5)
!
          call xctmr0(53)
!
! ---     output to mean archive file
!
          call mean_end(dtime)
          call forday(time_ave,yrflag, iyear,jday,ihour)
          call mean_archiv(n, iyear,jday,ihour)
!
          call mean_zero(dtime)
          call mean_add(n, 0.5)
!
          call xctmr1(53)
        endif  ! histmn
      endif  !meanfq
!
      if (histry .or. hiprof .or. hitile .or. hisurf .or. lfatal) then
        call xctmr0(53)
!
! ---   output to archive file
!
        call forday(dtime,yrflag, iyear,jday,ihour)
!
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          if (isopyc .or. mxlkrt) then
            do i=1,ii
              if (SEA_U) then
                umix(i,j)=u(i,j,1,n)
              endif !iu
              if (SEA_V) then
                vmix(i,j)=v(i,j,1,n)
              endif !iv
            enddo !i
          endif !isopyc .or. mxlkrt
!
          if (mxl_no) then
            do i=1,ii
              if (SEA_P) then
                  tmix(i,j) = temp(i,j,1,n)
                  smix(i,j) = saln(i,j,1,n)
                 thmix(i,j) = th3d(i,j,1,n)
                mixflx(i,j) = 0.0
                buoflx(i,j) = 0.0
                bhtflx(i,j) = 0.0
              endif !ip
              if (SEA_U) then
                   umix(i,j) =    u(i,j,1,n)
              endif !iu
              if (SEA_V) then
                   vmix(i,j) =    v(i,j,1,n)
              endif !iv
            enddo !i
          endif !mxl_no
!
          if (histry) then
            do k= 1,kk
              do i=1,ii
                if (SEA_P) then
! ---             convert diapycnal thickness changes into 
! ---             actual interface fluxes
                  if (k.gt.1) then
                    diaflx(i,j,k)=diaflx(i,j,k)/(2.*onem) + &
                                  diaflx(i,j,k-1)
                  else
                    diaflx(i,j,k)=diaflx(i,j,k)/(2.*onem)
                  endif
                endif !ip
              enddo !i
            enddo !k
          endif !histry
        enddo !j
!$OMP   END PARALLEL DO
!
        if     (lfatal) then  !write archive and exit
          if (mnproc.eq.1) then
          write (intvl,'(i3.3)') 0
          endif !1st tile
          call archiv(n, kk, iyear,jday,ihour, intvl)
          call xcstop('(hycom)')
                 stop '(hycom)'   !won't get here
        else
          if (mnproc.eq.1) then
          write (intvl,'(i3.3)') int(dtime-timav+dsmall)
          endif !1st tile
          if (hisurf) then
            call archiv(n, 1,  iyear,jday,ihour, intvl)
          endif
          if (histry) then
            call archiv(n, kk, iyear,jday,ihour, intvl)
          endif
          if     (hiprof) then
            call archiv_prof(n, kk, iyear,jday,ihour)
          endif
          if     (hitile) then
            call archiv_tile(n, kk, iyear,jday,ihour)
          endif
        endif  !lfatal:else
!
        if (histry) then
!$OMP     PARALLEL DO PRIVATE(j,k,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do k= 1,kk
              do i=1,ii
                if (SEA_P) then
                  diaflx(i,j,k)=0.0
                endif !ip
              enddo !i
            enddo !k
          enddo !j
!$OMP     END PARALLEL DO
!
          timav=time
        endif
        call xctmr1(53)
      endif  ! histry.or.hiprof.or.hitile.or.hisurf.or.lfatal
!
#if defined(USE_ESMF4)
      if     (put_export) then
!
! ---   fill the ESMF Export State.
!
        call Export_ESMF
        call forday(dtime,yrflag, iyear,jday,ihour)
        call Archive_ESMF(iyear,jday,ihour)
      endif
#endif
!
#if defined (USE_NUOPC_CESMBETA)
! --- end of coupling sequence
      ! stepable end-of-run condition
      ! endtime present means that "end_of_run_cpl" needs to be set correctly
#if defined (DMI_CICE_COUPLED)
      end_of_run_cpl = .true.
#else
      nstep2_cpl=nint(endtime*(86400.0d0/baclin))

      if(mnproc.eq.1) print *,"mod_hycom,nstep,nstep2,nstep2_cpl =", &
          nstep,nstep2,nstep2_cpl
      if(mnproc.eq.1) print *,"endtime =", endtime

      if(nstep.ge.nstep2_cpl) then
        end_of_run_cpl = .true.
      else
        end_of_run_cpl = .false.
      endif

!!Alex
!! add mean export field
      inv_cplifq= 1./icefrq
#endif
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          if (SEA_P) then
             sshm(i,j) = sshm(i,j) + srfhgt(i,j)

             um(i,j)  = um(i,j)  + 0.5*(    u(i,j,1,1)+    u(i,j,1,2)) 
             ubm(i,j) = ubm(i,j) + 0.5*(ubavg(i,j,  1)+ubavg(i,j,  2))

             vm(i,j)  = vm(i,j)  + 0.5*(    v(i,j,1,1)+    v(i,j,1,2)) 
             vbm(i,j) = vbm(i,j) + 0.5*(vbavg(i,j,  1)+vbavg(i,j,  2)) 

             tavgm(i,j) = tavgm(i,j) + 0.5*(temp(i,j,1,1) &
                                           +temp(i,j,1,2))
             savgm(i,j) = savgm(i,j) + 0.5*(saln(i,j,1,1) &
                                           +saln(i,j,1,2))

          endif
        enddo
      enddo
!$OMP END PARALLEL DO
      ntavg = ntavg + 1

      if (end_of_run_cpl) then

        sshm(:,:)=srfhgt(:,:)
!Alex  calculation of seas surface slope for export to CICE (NUOPC)
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
          if (SEA_P) then
           ssh_e = 0.0
           ssh_w = 0.0
           ssh_n = 0.0
           ssh_s = 0.0
           dhdx  = 0.0
           dhdy  = 0.0
           maskw = 0.0
           maske = 0.0
           maskn = 0.0
           masks = 0.0
           ! calculate the eastward slope
           if (ip(i-1,j).ne.0) then
             ssh_w = (srfhgt(i  ,j) - srfhgt(i-1,j))/(g*scux(i  ,j))
             maskw = 1.0
           endif
           if (ip(i+1,j).ne.0) then
             ssh_e = (srfhgt(i+1,j) - srfhgt(i  ,j))/(g*scux(i+1,j))
             maske = 1.0
           endif

           if ( maskw .eq.1.0 .or. maske .eq. 1.0 ) then
              dhdx=(ssh_e+ssh_w)/(maskw+maske) !! on the p-grid
           endif

           ! calculate the northward slope
           if (ip(i,j-1).ne.0) then
              ssh_s = (srfhgt(i,j  ) - srfhgt(i,j-1))/(g*scvy(i,j  ))
              masks = 1.0
           endif
           if (ip(i,j+1).ne.0) then
              ssh_n = (srfhgt(i,j+1) - srfhgt(i,j  ))/(g*scvy(i,j+1))
              maskn = 1.0
           endif
           if (masks .eq. 1.0 .or. maskn .eq. 1.0) then
              dhdy=(ssh_n+ssh_s)/(maskn+masks) !! on the p-grid
           endif
           ! convert to eastward/northward grid
           dhde(i,j) = (dhdx*cos(pang(i,j)) + dhdy*sin(-pang(i,j)))/ntavg
           dhdn(i,j) = (dhdy*cos(pang(i,j)) - dhdx*sin(-pang(i,j)))/ntavg
         endif
         enddo
        enddo
!$OMP END PARALLEL DO
!!Alex calculation of u and v surf for export to CICE
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
           if (SEA_P) then
! ---      average currents over top thkcdw meters
             utot = 0.5*(um(i  ,j) + ubm(i  ,j) + &
                         um(i+1,j) + ubm(i+1,j))
             vtot = 0.5*(vm(i,j  ) + vbm(i,j  ) + &
                         vm(i,j+1) + vbm(i,j+1))

             uml(i,j)=(utot*cos(pang(i,j)) + vtot*sin(-pang(i,j)))/ntavg
             vml(i,j)=(vtot*cos(pang(i,j)) - utot*sin(-pang(i,j)))/ntavg
           endif
          enddo
        enddo
!$OMP END PARALLEL DO

!!Alex calculation of tmxl,smxl
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
           if (SEA_P) then
              tml(i,j) = tavgm(i,j)/ntavg
              sml(i,j) = savgm(i,j)/ntavg
           endif
          enddo
        enddo
!$OMP END PARALLEL DO
!        if (ntavg.ne.icefrq) then
!          if (ntavg.eq.2*icefrq) then !! first coupling cycle ...
!            tml(:,:)=tml(:,:)*(icefrq/ntavg)
!            sml(:,:)=sml(:,:)*(icefrq/ntavg)
!            uml(:,:)=uml(:,:)*(icefrq/ntavg)
!            vml(:,:)=vml(:,:)*(icefrq/ntavg)
!          else
!            if (mnproc.eq.1) then
!              write(lp,*) 
!              write(lp,*) ' icefrq and cplfrq should be the same ....'
!              write(lp,*) ' icefrq, ntavg:',icefrq, ntavg
!              call flush(lp)
!            endif !1st tile
!            call xcstop('(hycom)')
!          endif
!        endif

! --- reset average fields
        ntavg = 0
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            if (SEA_P) then
               sshm(i,j) = 0.0
                um(i,j)  = 0.0
                ubm(i,j) = 0.0
                vm(i,j)  = 0.0
                vbm(i,j) = 0.0
              tavgm(i,j) = 0.0
              savgm(i,j) = 0.0
            endif
          enddo
        enddo
!$OMP END PARALLEL DO

      endif ! end_of_run_cpl
#endif /* USE_NUOPC_CESMBETA */
! --- restart
!TILL 28/3  I am a bit in doubt here whether it is correct
      if (present(restart_write) ) then
          restart_cpl = restart_write .and. end_of_run_cpl
      else
          restart_cpl = .false.
      endif
      if (restrt .or. restart_cpl) then

        call xctmr0(51)
!
! ---   output to restart and flux statitics files
!
        if (mxlkpp) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                dpmixl(i,j,m) = dpmixl(i,j,n)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
        endif
!
        call forday(dtime,yrflag, iyear,jday,ihour)
        if (mnproc.eq.1) then
        write (lp,'(a,i9, 9x,a,i6.4, 9x,a,i5.3, 9x,a,i4.2)') &
         &' time step',nstep, &
         &'y e a r',   iyear, &
         &'d a y',     jday, &
         &'h o u r',   ihour
        call flush(lp)
        endif !1st tile
!
        if (restart_cpl) then
          write(flnmra,'(a,i4.4,a,i3.3,a,i2.2)') &
                      &'restart_', iyear,'-',jday,'-',ihour
          open(1,file=trim(pointer_filename),form='formatted', &
                 status='unknown')
          write(1,'(a)')trim(flnmra)
          close(1)
          flnmrb = trim(flnmra)
        else
          flnmra = flnmrso  !.a extension added by restart_out
          flnmrb = flnmrso  !.b extension added by restart_out
        endif
!
        call restart_out(nstep,dtime, flnmra,flnmrb, nstep.ge.nstep2, &
                         restart_cpl )
!
! ---   set layer thickness (incl.bottom pressure) at u,v points
! ---   needed because restart_out may have modified dp
!
        call dpthuv
!
        call xctilr(    dp(1-nbdy,1-nbdy,1,1),1,2*kk, nbdy,nbdy, halo_ps)
        call xctilr(dpmixl(1-nbdy,1-nbdy,1  ),1,2   , nbdy,nbdy, halo_ps)
!
        margin = nbdy
!
        nstep = nstep+1  !for pipe_compare
        do nm=1,2
!$OMP     PARALLEL DO PRIVATE(j,k,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            if     (nm.eq.mod(nstep+1,2)+1) then
              do i=1-margin,ii+margin
                if (SEA_P) then
                  dpbl( i,j)=dpmixl(i,j,nm)
                endif !ip
              enddo !i
            endif !nm
            do i=1-margin,ii+margin
              if (SEA_P) then
                p(i,j,1)=0.0
                do k=1,kk
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,nm)
                enddo !k
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
          call dpudpv(dpu(1-nbdy,1-nbdy,1,nm), &
                      dpv(1-nbdy,1-nbdy,1,nm), &
                      p,depthu,depthv, margin,max(0,margin-1))
!            
        enddo  !nm=1,2
        nstep = nstep-1  !restore
!
! ---   update montg to preserve bit for bit reproducability on restart
        m=mod(nstep  ,2)+1
        n=mod(nstep+1,2)+1
        if     (meanfq.ne.0.0) then
! ---     assume meanfq not changed to zero when rerun from restart
          call momtum_hs(n,m)
          surflx(:,:) = 0.0
          salflx(:,:) = 0.0
          wtrflx(:,:) = 0.0
          if     (histmn) then
            call mean_zero(dtime)
            call mean_add(n, 0.5)
          endif  ! histmn
        endif !meanfq
        call momtum_hs(n,m)
!
        call xctmr1(51)
      endif  ! restrt
!
      if     (histry .or. hiprof .or. hitile .or. hisurf) then
        if (mnproc.eq.1) then
        write (lp,'(a,i9,a,f13.5,a)') &
         &' step',nstep,' day',dtime,' -- archiving completed --'
        call flush(lp)
        endif !1st tile
      endif
!
! --- read next incremental update.
!
      if (incflg.ne.0) then
        call incupd_rd(dtime)
      endif
#if defined(STOKES)
!
! --- calculate Stokes Drift for the next time step
      if     (nsdzi.gt.0) then
        call stokes_forfun(dtime,n)
      endif !nsdzi
#endif
!
      delt1=baclin+baclin
!
#if defined (ESPC_COUPLE)
      if(nstep.ge.nstep2_cpl) then
        end_of_run_cpl = .true.
        nstep1_cpl=nstep2_cpl
      else
        end_of_run_cpl = .false.
      endif
#endif

      end_of_run = nstep.ge.nstep2
!
! --- at end: output float restart file 
!
      if     (synflt .and. end_of_run) then
        call floats_restart
      endif !synflt+end_of_run
#if defined(USE_ESMF4) || defined (ESPC_COUPLE)
      call xctmr1(80)   !time  inside HYCOM_Run
      call xctmr0(79)   !time outside HYCOM_Run
#endif
      end subroutine HYCOM_Run

#if defined(USE_ESMF4)
      subroutine HYCOM_Final &
                 (gridComp, impState, expState, extClock, rc)
!
! --- Calling parameters
      type(ESMF_GridComp)  :: gridComp
      type(ESMF_State)     :: impState
      type(ESMF_State)     :: expState
      type(ESMF_Clock)     :: extClock
      integer, intent(out) :: rc
!
! --- Report
      call ESMF_LogWrite("HYCOM finalize routine called", &
           ESMF_LOG_INFO, rc=rc)
      call ESMF_LogFlush(rc=rc)
!
! --- Destroy internal ocean clock
!     call ESMF_ClockDestroy(intClock, rc=rc)
!
! --- print active timers.
      call xctmr1(79)   !time outside HYCOM_Run
      call xctmrp
!
      if (mnproc .eq. 1) then
        write(nod,'(a)') 'normal stop'
        call flush(nod)
      endif
!
      end subroutine HYCOM_Final
#else
      subroutine HYCOM_Final

# if defined (ESPC_COUPLE)

!
! --- print active timers.
      call xctmr1(79)   !time outside HYCOM_Run
      call xctmrp
#endif
!
! --- end of the run.
      if (mnproc.eq.1) then
      write(nod,'(a)') 'normal stop'
      call flush(nod)
      endif !1st tile
# if ! defined (ESPC_COUPLE)
      call xcstop('(normal)')  !calls xctmrp
             stop '(normal)'   !won't get here
#endif
      end subroutine HYCOM_Final
#endif /* USE_ESMF4:else */

      end module mod_hycom
!
!
!> Revision history:
!>
!> May  1997 - removed statement "theta(1)=-thbase" after loop 14
!> June 1997 - added loop 60 to fix bug in vertical summation of -diaflx-
!> Oct. 1999 - option for krt or kpp mixed layer model - convec and diapfl
!>             not called for kpp mixing model
!> Oct. 1999 - dpbl (boundary layer thickness) is output in addition to
!>             dpmixl when the kpp mixing model is selected
!> May  2000 - conversion to SI units
!> Aug. 2000 - added isopycnic (MICOM) vertical coordinate option
!> Oct. 2000 - added option for high frequency atmospheric forcing
!> Nov. 2000 - archive time stamp is either time step or YYYY_DDD_HH
!> Aug. 2002 - added support for multiple tracers
!> Nov. 2002 - more basin-wide surface flux statistics
!> Dec. 2003 - more basin-wide mean statistics
!> Mar. 2005 - more accurate ice statistics
!> Jan. 2006 - mod_hycom with HYCOM_Init, HYCOM_Run, HYCOM_Final
!> Nov. 2006 - version 2.2
!> Nov. 2006 - added incremental update (for data assimilation)
!> Mar. 2007 - added srfhgt
!> Mar. 2010 - removed  DETIDE
!> Apr. 2010 - put back DETIDE
!> Apr. 2010 - added proffq
!> Apr. 2010 - added archiv_init
!> Mar. 2011 - mean archive now the same with and without tides
!> Apr. 2011 - surface archives separate from 3-D archives
!> Nov. 2011 - can now start from climatology when yrflag=3
!> Jan. 2012 - smooth imported ice drift and exported ocean currents
!> Jan. 2012 - added thkcdw
!> Aug. 2012 - use CICE fields for ice statistics when available
!> Aug. 2012 - call pipe_init after blkdat
!> Dec. 2012 - initialize l0-3 and w0-3 before first call to momtum_hs
!> Jan. 2013 - replace dragrh with drgten
!> Apr. 2013 - added incupd_si_c
!> Aug. 2013 - optionally added mod_stokes
!> Nov. 2013 - added jerlv0=-1 and forfunc
!> Jan. 2014 - added mslprf and forfunhp
!> Apr. 2014 - added ice shelf logic (ishlf)
!> Apr. 2014 - replace ip with ipa for mass sums, reduces need for jja
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Feb. 2015 - removed dtime0 from calculation of dtime
!> Feb. 2015 - fixed issues with bit for bit reproducability on restart
!> Sep. 2015 - NaN test on mean surf thk (SST and SSS)
!> Dec. 2015 - time all instances of HY_Ini, HY_Out, HY_Run
!> Aug. 2017 - added restrt to call to incupd
!> Sep. 2017 - added call to stfupd
!> Aug. 2018 - always use mass-conserving diagnostics, i.e. include oneta
!> Aug. 2018 - added mod_asselen, call asselin_filter before hybgen
!> Aug. 2018 - tsadvc now called after barotp
!> Aug. 2018 - added mod_barotp and call to barotp_init
!> Aug. 2018 - added wtrflx
!> Nov  2018 - added yrflag=4 for 365 days no-leap calendar (CESM)
!> Dec. 2018 - added /* USE_NUOPC_CESMBETA */ macro for CESMBETA coupled simulation
!> Dec. 2018 - added /* ESPC_COUPLE */ macro for coupling with NAVYESPC
!> Feb. 2019 - replaced onetai by 1.0
!> Sep. 2019 - added oneta0
