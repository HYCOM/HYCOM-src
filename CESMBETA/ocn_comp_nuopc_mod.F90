module ocn_comp_nuopc_mod

  !-----------------------------------------------------------------------------
  ! OCN Component for CESM-BETA; use ESMF and  NUOPC 
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS      => SetServices, &
    model_label_Advance   => label_Advance
  
  use MOD_HYCOM, only : HYCOM_Init, HYCOM_Run, HYCOM_Final, &
    end_of_run, end_of_run_cpl
    
  use mod_hycom_nuopc_glue
  use mod_cb_arrays_nuopc_glue
  use ESMF_IOScripMod
  use esmfshr_util_mod, only : esmfshr_util_ArrayGetIndex
  use esmfshr_util_mod, only : esmfshr_util_ArrayGetSize
  use esmf2mct_mod    , only : esmf2mct_init

  use mct_mod,          only : mct_gsMap, mct_gGrid, mct_gsMap_gsize
  use shr_string_mod,   only : shr_string_listGetNum
  use seq_flds_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use esmfshr_nuopc_mod

  implicit none
  
  private
  
  ! private internal state to keep instance data
  type InternalStateStruct
    type(hycom_nuopc_glue_type)   :: glue
    integer                       :: slice
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type
  
  public SetServices


  type HYCOM_CESM_FIELD_TABLE_ENTRY
    character(len=128)              :: hycom_stdname
    character(len=128)              :: cesm_stdname
    character(len=32)               :: unit
    logical                         :: connected
  end type

  public loadHycomDictionary

  type(mct_gsMap), public, pointer  :: gsmap_o
  type(mct_gGrid), public, pointer  :: dom_o
  integer                           :: OCNID      
  type(ESMF_Routehandle), save      :: HYCOM2CESM_RHR8, CESM2HYCOM_RHR8
  type(ESMF_Routehandle), save      :: HYCOM2CESM_RHI4, CESM2HYCOM_RHI4
  integer, parameter, public        :: number_import_fields = 30
  integer, parameter, public        :: number_export_fields = 8

  type(HYCOM_CESM_FIELD_TABLE_ENTRY):: cesm2hycom_table(number_import_fields)
  type(HYCOM_CESM_FIELD_TABLE_ENTRY):: hycom2cesm_table(number_export_fields)

  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Provide InitializeP0 to switch from default IPDv00 to IPDv02
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
!    ! overwrite Finalize
!    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
!      userRoutine=Finalize, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Run routine to execute run functionality for tight coupling
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      userRoutine=routine_Run2, phaseLabelList=(/"OCN_TIGHT_RUN_PHASE"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Run routine to execute run functionality for loose coupling
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      userRoutine=routine_Run3, phaseLabelList=(/"OCN_LOOSE_RUN_PHASE"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS

!#ifndef HYCOM_IN_CESM
!    ! importable fields:
!    call NUOPC_Advertise(importState, &
!      StandardNames=(/ &
!      "surface_downward_eastward_stress       ",    & ! from ATM
!      "surface_downward_northward_stress      ",    & ! from ATM
!      "wind_speed_height10m                   ",    & ! from ATM
!      "friction_speed                         ",    & ! from ATM
!      "mean_net_sw_flx                        ",    & ! from ATM
!      "mean_down_lw_flx                       ",    & ! from ATM
!      "mean_up_lw_flx                         ",    & ! from ATM
!      "mean_lat_flx                           ",    & ! from ATM
!      "mean_sens_flx                          ",    & ! from ATM
!      "inst_temp_height2m                     ",    & ! from ATM
!      "mean_prec_rate                         ",    & ! from ATM
!      "inst_spec_humid_height2m               ",    & ! from ATM
!      "sea_surface_temperature                ",    & ! from ATM
!      "water_flux_into_sea_water              ",    & ! from ATM
!      "frozen_water_flux_into_sea_water       ",    & ! from ATM
!      "sea_ice_area_fraction                  ",    & ! from SEA-ICE
!      "downward_x_stress_at_sea_ice_base      ",    & ! from SEA-ICE
!      "downward_y_stress_at_sea_ice_base      ",    & ! from SEA-ICE
!      "downward_sea_ice_basal_solar_heat_flux ",    & ! from SEA-ICE
!      "upward_sea_ice_basal_heat_flux         ",    & ! from SEA-ICE
!      "downward_sea_ice_basal_salt_flux       ",    & ! from SEA-ICE
!      "downward_sea_ice_basal_water_flux      ",    & ! from SEA-ICE
!      "sea_ice_temperature                    ",    & ! from SEA-ICE
!      "sea_ice_thickness                      ",    & ! from SEA-ICE
!      "sea_ice_x_velocity                     ",    & ! from SEA-ICE
!      "sea_ice_y_velocity                     "/),  & ! from SEA-ICE
!      rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    
!    ! exportable fields:
!    call NUOPC_Advertise(exportState, &
!      StandardNames=(/ &
!      "sea_surface_temperature                  ",    &
!      "upward_sea_ice_basal_available_heat_flux ",    &
!      "sea_lev                                  ",    &
!      "mixed_layer_depth                        ",    &
!      "s_surf                                   ",    &
!      "eastward_sea_surface_slope               ",    &
!      "northward_sea_surface_slope              ",    &
!      "ocn_current_zonal                        ",    &
!      "ocn_current_merid                        "/),  &
!      rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!#endif

    ! set Component name so it becomes identifiable
    call ESMF_GridCompSet(gcomp, name="HYCOM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Load HYCOM Field dictionary into CESM NUOPC framework, should be done at driver level
    call loadHycomDictionary(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! Initialize 2-way connection tables between HYCOM 2D and CESM 1D layers
    call initialize_HYCOM_CESM_tables(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Setup import-able CESM 1D fields according to names in ocn_import_fields
    ! Same names are used in HYCOM_CESM_tables
    call esmfshr_nuopc_advertise_fields( &
      ocn_import_fields, importState, tag='HYCOM import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Setup export-able CESM 1D fields according to names in ocn_export_fields
    call esmfshr_nuopc_advertise_fields( &
      ocn_export_fields, exportState, tag='HYCOM export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
      
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! This method assumes that the incoming clock indicates the full simulation
    ! interval between startTime and stopTime.
    
    ! local variables    
    type(ESMF_Field)            :: field
    type(ESMF_Grid)             :: gridIn
    type(ESMF_Grid)             :: gridOut
    type(ESMF_array)            :: array
    type(ESMF_VM)               :: vm
    integer                     :: mpiComm
    TYPE(ESMF_Time)             :: startTime, stopTime, hycomRefTime, currTime
    TYPE(ESMF_TimeInterval)     :: interval,timeStep
    real(ESMF_KIND_R8)          :: startTime_r8, stopTime_r8
    type(InternalState)         :: is
    integer                     :: stat
    type(ESMF_CALKIND_FLAG)     :: calkind
    logical                     :: restFlag = .false.        ! initial/restart run (F/T)
    character(len=80)           :: pointer_filename          ! restart pointer file !!Alex
    logical                     :: restart_write = .false.   ! write restart
    integer                     :: cplfrq
    type(ESMF_State)            :: import_state, export_state
    type(ESMF_Mesh)             :: mesh, meshIn, meshOut
    type(ESMF_DistGrid)         :: distgrid, distgrid2D
    type(ESMF_Array)            :: o2x, x2o, dom
    type(ESMF_ArraySpec)        :: arrayspec
    type(ESMF_Delayout)         :: delayout
    integer                     :: ldeCount, eleCount, lsize, mpicom_ocn, nfields, lde, i, n_elem
    integer                     :: maxIndex(2, 1)
    integer, pointer            :: fptrSeqIndex(:)
    character(len=32)           :: starttype                 ! infodata start type
    real(ESMF_KIND_R8)          :: l_startTime_r8
    
    rc = ESMF_SUCCESS
    
    ! Allocate memory for the internal state and set it in the Component.
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of the internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Prepare to call into HYCOM_Init
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpiComm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Translate currTime and stopTime into HYCOM format
    call ESMF_ClockGet(clock, currTime=currTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Define HYCOM Ref Time
!!Alex    call ESMF_TimeSet(hycomRefTime, yy=1901, mm=01, dd=01, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    call ESMF_TimeSet(hycomRefTime, yy=0001, mm=01, dd=01, calkindflag=ESMF_CALKIND_NOLEAP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    interval = currTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=startTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    interval = stopTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=stopTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    print *, " HYCOM_INIT -->> startTime_r8=", startTime_r8, "stopTime_r8=", stopTime_r8

    ! get coupling frequency from ocean clock for 1st export
     call ESMF_ClockGet(ccsm_Eclock_o, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_TimeIntervalGet(timeStep, d_r8=ocn_cpl_frq, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return  ! bail out
       
    ! Get start type run : start-up or continuous run
    call ESMF_AttributeGet(exportState, name="start_type", value=starttype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       l_startTime_r8=-startTime_r8
       restFlag = .false.
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       l_startTime_r8=startTime_r8
       restFlag = .true.
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       l_startTime_r8=startTime_r8
       restFlag = .true.
    else
       call ESMF_LogWrite('hycom_nuopc ERROR: unknown starttype',  &
          ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_OBJ_BAD
       return
    end if

    ! Get the pointer restart name !!Alex
    pointer_filename = 'rpointer.ocn' 
    restart_write = seq_timemgr_RestartAlarmIsOn(ccsm_EClock_o)

    call ESMF_LOGWRITE("BEFORE HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)
    
    ! Call into the HYCOM initialization  
    call HYCOM_Init(mpiComm, & ! -->> call into HYCOM <<--
       hycom_start_dtg=l_startTime_r8, hycom_end_dtg=stopTime_r8, &
       pointer_filename=pointer_filename,  restart_write=restart_write)

    call ESMF_LOGWRITE("AFTER HYCOM_INIT", ESMF_LOGMSG_INFO, rc=rc)

    ! Test if coupling frequency is consistent between HYCOM and framework
    cplfrq = nint( ocn_cpl_frq*(86400.0/baclin) )
    if (icefrq .ne. cplfrq) then 
       call ESMF_LogWrite('hycom_nuopc ERROR: coupling frq between HYCOM and &
   &   NUOPC not consistent ',  &
          ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_OBJ_BAD
       return
    end if
    
    ! Fill in the glue structure.
    call HYCOM_GlueInitialize(is%wrap%glue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Write some HYCOM distribution info into the Log.
    call HYCOM_TileInfo(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! use the HYCOM Grid that was setup inside of the glue structure
    gridIn  = is%wrap%glue%grid ! for imported Fields
    gridOut = is%wrap%glue%grid ! for exported Fields
    
    ! conditionally realize or remove Fields in import and export States
    ! also keep track of these Fields on the glue layer
    
    ! importable fields:
    call HYCOM_GlueFieldsRealize(is%wrap%glue, importState, &
      StandardNames=(/ &
      "surface_downward_eastward_stress       ",    & ! from ATM
      "surface_downward_northward_stress      ",    & ! from ATM
      "wind_speed_height10m                   ",    & ! from ATM
      "friction_speed                         ",    & ! from ATM
      "mean_net_sw_flx                        ",    & ! from ATM
      "mean_down_lw_flx                       ",    & ! from ATM
      "mean_up_lw_flx                         ",    & ! from ATM
      "mean_lat_flx                           ",    & ! from ATM
      "mean_sens_flx                          ",    & ! from ATM
      "inst_temp_height2m                     ",    & ! from ATM
      "mean_prec_rate                         ",    & ! from ATM
      "inst_spec_humid_height2m               ",    & ! from ATM
      "sea_surface_temperature                ",    & ! from ATM
      "water_flux_into_sea_water              ",    & ! from ATM
      "frozen_water_flux_into_sea_water       ",    & ! from ATM
      "sea_ice_area_fraction                  ",    & ! from SEA-ICE
      "downward_x_stress_at_sea_ice_base      ",    & ! from SEA-ICE
      "downward_y_stress_at_sea_ice_base      ",    & ! from SEA-ICE
      "downward_sea_ice_basal_solar_heat_flux ",    & ! from SEA-ICE
      "upward_sea_ice_basal_heat_flux         ",    & ! from SEA-ICE
      "downward_sea_ice_basal_salt_flux       ",    & ! from SEA-ICE
      "downward_sea_ice_basal_water_flux      ",    & ! from SEA-ICE
      "sea_ice_temperature                    ",    & ! from SEA-ICE
      "sea_ice_thickness                      ",    & ! from SEA-ICE
      "sea_ice_x_velocity                     ",    & ! from SEA-ICE
      "sea_ice_y_velocity                     "/),  & ! from SEA-ICE
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !! exportable fields:
    call HYCOM_GlueFieldsRealize(is%wrap%glue, exportState, &
      StandardNames=(/ &
      "sea_surface_temperature                  ",    &
      "upward_sea_ice_basal_available_heat_flux ",    &
      "sea_lev                                  ",    &
      "mixed_layer_depth                        ",    &
      "s_surf                                   ",    &
      "eastward_sea_surface_slope               ",    &
      "northward_sea_surface_slope              ",    &
      "ocn_current_zonal                        ",    &
      "ocn_current_merid                        "/),  &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Import data to HYCOM native structures through glue fields.
   call HYCOM_GlueFieldsDataImport(is%wrap%glue, .not. restFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Connect 2D glue layer with 1D cap so internal HYCOM data will be copied to 2D Fields in glue layer
    do i  = 1, number_export_fields
      if(hycom2cesm_table(i)%connected) then
        call HYCOM_RedistHYCOM2CESM(is%wrap%glue%exportFields, hycom2cesm_table(i)%hycom_stdname, &
          exportState, hycom2cesm_table(i)%cesm_stdname, connectOnly=.true., rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo
    ! Export HYCOM native structures to data through glue fields.
    CALL HYCOM_GlueFieldsDataExport(is%wrap%glue, .not. restFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Reset the slice counter
    is%wrap%slice = 1


    ! Export precip_fact for precipitation adjustment !!Alex
!    PRINT*,'Alex precip_fact=',pcp_fact
    call ESMF_AttributeSet(exportState, name="precip_fact", value=pcp_fact, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call seq_infodata_PutData(infodata, precip_fact=pcp_fact)

    ! Prepare for CESM SPECIFIC DATA STRUCTURES
    import_state = importState
    export_state = exportState

    ! duplicate the mpi communicator from the current VM 
    call MPI_Comm_dup(mpicomm, mpicom_ocn, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! satisfy some external needs
    call ESMF_AttributeGet(export_state, name="ID", value=OCNID, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_GridGet(is%wrap%glue%grid, distgrid=distgrid2D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_DistGridGet(distgrid2D, 0, elementCount=n_elem, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    allocate(fptrSeqIndex(n_elem))
    call ESMF_DistGridGet(distgrid2D, 0, seqIndexList=fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    distgrid = ESMF_DistGridCreate(fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    
    deallocate(fptrSeqIndex)

    meshIn = ESMF_MeshCreate(distgrid, nodalDistgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    meshOut = meshIn

    !-----------------------------------------
    !  Set arrayspec for dom, o2x and x2o
    !-----------------------------------------
    
    call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_DistGridGet(distgrid, delayout=delayout, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_DelayoutGet(delayout, localDeCount=ldeCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !print *, 'HYCOM DG DELAYOUT localDECount: ', ldeCount

    lsize = 0
    do lde = 0, ldeCount-1
      call ESMF_DistGridGet(distgrid, lde, elementCount=eleCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      lsize = lsize + eleCount
    enddo

    !print *, 'HYCOM DG DELAYOUT lsize: ', lsize

    !-----------------------------------------
    ! Create dom 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_dom_fields))

    dom = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="domain", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    ! Set values of dom (needs ocn initialization info)

    call ocn_domain_esmf(dom, is%wrap%glue%grid, rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
   
    !----------------------------------------- 
    !  Create o2x 
    !-----------------------------------------

    ! 1d undistributed index of fields, 2d is packed data

    nfields = shr_string_listGetNum(trim(seq_flds_o2x_fields))

    o2x = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="d2x", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(o2x, name="mct_names", value=trim(seq_flds_o2x_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    !  Create x2o 
    !-----------------------------------------

    nfields = shr_string_listGetNum(trim(seq_flds_x2o_fields))

    x2o = ESMF_ArrayCreate(distgrid=distgrid, arrayspec=arrayspec, distgridToArrayMap=(/2/), &
         undistLBound=(/1/), undistUBound=(/nfields/), name="x2d", rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2o, name="mct_names", value=trim(seq_flds_x2o_fields), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !----------------------------------------- 
    ! Add esmf arrays to import and export state 
    !-----------------------------------------

    call ESMF_StateAdd(export_state, (/dom/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(export_state, (/o2x/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    call ESMF_StateAdd(import_state, (/x2o/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate(gsmap_o)
    allocate(dom_o)
   
    call esmf2mct_init(distgrid, OCNID, gsmap_o, mpicom_ocn, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call esmf2mct_init(dom, dom_o, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !print *, 'HYCOM GSMAP gsize: ', mct_gsMap_gsize(gsmap_o)

    call ESMF_GridGet(is%wrap%glue%grid, distgrid=distgrid2D, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_DistGridGet(distgrid2D, maxIndexPTile=maxIndex, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="gsize", value=mct_gsMap_gsize(gsmap_o), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocnrof_prognostic", value=.true., rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_nx", value=maxIndex(1,1), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_AttributeSet(export_state, name="ocn_ny", value=maxIndex(2,1), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !! Create and Realize Importable fields
    call esmfshr_nuopc_create_fields( &
      ocn_import_fields, meshIn, importState, tag='HYCOM import', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Create and Realize Exportable fields
    call esmfshr_nuopc_create_fields( &
      ocn_export_fields, meshOut, exportState, tag='HYCOM export', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! Redist Export Fields from internal HYCOM to CESM 1D Fields, but zero dhdx, dhdy, Fioo_q, u, v, bldepth
    do i  = 1, number_export_fields
      if(hycom2cesm_table(i)%connected) then
        call HYCOM_RedistHYCOM2CESM(is%wrap%glue%exportFields, hycom2cesm_table(i)%hycom_stdname, &
          exportState, hycom2cesm_table(i)%cesm_stdname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo
    ! Copy data in individual 1D Fields to o2x Array. Need this during Initialization as CESM Driver accesses this directly.
    call esmfshr_nuopc_copy(ocn_export_fields, exportState, 'd2x', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !! Redist 1D Field to 2D Field and write result in .nc Files
    !call RedistAndWriteField(is%wrap%glue%grid, exportState, filePrefix="field_ocn_init_export_", &
    !  timeslice=1, relaxedFlag=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    !call HYCOM_WriteFieldBundle(is%wrap%glue%grid, is%wrap%glue%exportFields, filePrefix="fieldbundle_ocn_init_export_", &
    !  timeslice=1, relaxedFlag=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    
  end subroutine

  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
!#ifndef HYCOM_IN_CESM
!    call HYCOM_ModelAdvance(gcomp, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!    return ! bail out
!#endif

  end subroutine

  SUBROUTINE HYCOM_ModelAdvance(gcomp, rc)
    use seq_infodata_mod
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep_O
    type(ESMF_Time)             :: hycomRefTime
    type(ESMF_TimeInterval)     :: interval
    real(ESMF_KIND_R8)          :: endTime_r8    ! end of coupling sequence
    type(InternalState)         :: is
    logical                     :: initFlag
    type(ESMF_CALKIND_FLAG)     :: calkind

    integer                     :: fieldCount, i
    character(len=128), allocatable :: fieldNameList(:)
    type(ESMF_Field)            :: field
    character(len=80)           :: pointer_filename     ! restart pointer file !!Alex
    logical                     :: restart_write = .false.
    character(len=32)           :: starttype            ! infodata start type
    logical                     :: restFlag = .false.
    character(len=128)          :: msg
    character*80 :: filenc !!Alex

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    ! on the internal Clock object. The NUOPC Layer will update the Clock 
    ! automatically each time before entering ModelAdvance(), but the HYCOM
    ! model must be stepped forward within this method.
    
    CALL ESMF_ClockGet(clock, currTime=currTime,  rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Translate currTime + timeStep into HYCOM format
!!Alex    call ESMF_TimeSet(hycomRefTime, yy=1901, mm=01, dd=01, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    call ESMF_TimeSet(hycomRefTime, yy=0001, mm=01, dd=01, calkindflag=ESMF_CALKIND_NOLEAP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get endtime from driver coupling frequency ocn_cpl
    ! There is an issue with ocn integration time, fix it by using global ocn
    ! clock. The ocn attempts to run on most frequent coupling frequency and 
    ! will only run when alarm is on. Integration should still use ocn coupling
    ! timestep.
    call ESMF_ClockGet(ccsm_Eclock_o, timeStep=timeStep_O, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
!!Alex incompatibility with ESMF 7.0.0    call ESMF_TimePrint(currTime, &
!      "--------------> HYCOM_Run() advancing from: ", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out
!    call ESMF_TimePrint(currTime + timeStep_O, &
!      "--------------> HYCOM_Run() advancing to: ", rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

    ! get endtime 
    interval = currTime - hycomRefTime
    call ESMF_TimeIntervalGet(interval, d_r8=endTime_r8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    endTime_r8 = endTime_r8 + ocn_cpl_frq
    
    call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
!!Alex    print *, 'FIELDBUNDLEPRINT: -> number of field in HYCOM import FieldBundle: ', fieldcount
    allocate(fieldNameList(fieldCount))
    call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    do i = 1, fieldCount
!!Alex      print *, 'FIELDBUNDLEPRINT: -> ', fieldNameList(i)
      call ESMF_FieldBundleGet(is%wrap%glue%importFields, fieldName=fieldNameList(i), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call ESMF_AttributeSet(field, name="HYCOM_IN_CESM_Connected", value=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    end do
    deallocate(fieldNameList)

    ! Redistribute 1D CESM import Fields to HYCOM Fields stored in glue import fieldbundle
    do i  = 1, number_import_fields
      if(cesm2hycom_table(i)%connected) then
        call HYCOM_RedistCESM2HYCOM(importState, cesm2hycom_table(i)%cesm_stdname, &
          is%wrap%glue%importFields, cesm2hycom_table(i)%hycom_stdname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo

    !! Redistribute 1D CESM import Fields to 2D Fields and write to .nc files
    !call RedistAndWriteField(is%wrap%glue%grid, importState, filePrefix="field_ocn_import_", &
    !  timeslice=is%wrap%slice, relaxedFlag=.true., overwrite=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    
    ! Get start type run : start-up or continuous run !!Alex add restFlag
    call ESMF_AttributeGet(exportState, name="start_type", value=starttype, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if (     trim(starttype) == trim(seq_infodata_start_type_start)) then
       restFlag = .false.
    else if (trim(starttype) == trim(seq_infodata_start_type_cont) ) then
       restFlag = .true.
    else if (trim(starttype) == trim(seq_infodata_start_type_brnch)) then
       restFlag = .true.
    else
       call ESMF_LogWrite('hycom_nuopc ERROR: unknown starttype',  &
          ESMF_LOGMSG_ERROR, rc=rc)
       rc = ESMF_RC_OBJ_BAD
       return
    end if

    !TODO: don't need the additional initialization step once data-dependency
    !TODO: is taken care of during initialize.
    initFlag = .false.
    if (is%wrap%slice==1 .and. (.not. restFlag)) initFlag = .true.
    
    
    ! Import data to HYCOM native structures through glue fields.
    call HYCOM_GlueFieldsDataImport(is%wrap%glue, initFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    
    !call ESMF_VMLogMemInfo('MEMORY Usage BEFORE HYCOM_RUN', rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
  
    ! Get the pointer restart name 
    pointer_filename = 'rpointer.ocn' 
    restart_write = seq_timemgr_RestartAlarmIsOn(ccsm_EClock_o)

!!Alex output some fields
!       write(filenc,'(a,a,a)') 'frzmlti_hycom.nc'
!       call ncput_2d(frzh(:,:),'frzmlt',TRIM(filenc))


    ! Enter the advancing loop over HYCOM_run...
    do
      ! ...on return the end-of-run flags indicate whether HYCOM has advanced
      ! far enough...
      CALL HYCOM_Run(endtime=endTime_r8,pointer_filename=pointer_filename, restart_write=restart_write) ! -->> call into HYCOM <<--
      !print *, "HYCOM_Run returned with end_of_run, end_of_run_cpl:", &
      !  end_of_run, end_of_run_cpl
      if (end_of_run .or. end_of_run_cpl) exit
    enddo

    !call ESMF_VMLogMemInfo('MEMORY Usage AFTER HYCOM_RUN', rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    ! Export HYCOM native data through the glue fields.
    call HYCOM_GlueFieldsDataExport(is%wrap%glue, initFlag, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Export precip_fact for precipitation adjustment !!Alex 
    call ESMF_AttributeSet(exportState, name="precip_fact", value=pcp_fact, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call seq_infodata_PutData(infodata, precip_fact=pcp_fact)

    !! Copy o2x Array to CESM 1D Fields if forcing is used.
    !call esmfshr_nuopc_copy(ocn_export_fields, 'd2x', exportState, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &

    ! Redistribute HYCOM field stored in glue export fieldbundle to CESM 1D Fields
    do i  = 1, number_export_fields
      if(hycom2cesm_table(i)%connected) then
        call HYCOM_RedistHYCOM2CESM(is%wrap%glue%exportFields, hycom2cesm_table(i)%hycom_stdname, &
          exportState, hycom2cesm_table(i)%cesm_stdname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    enddo
    ! Don't need this during Run Phase, no secret data flow during run phase
    !call esmfshr_nuopc_copy(ocn_export_fields, exportState, 'd2x', rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &

    !call HYCOM_WriteFieldBundle(is%wrap%glue%grid, is%wrap%glue%exportFields, filePrefix="fieldbundle_ocn_export_", &
    !  timeslice=is%wrap%slice, relaxedFlag=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out

    !call RedistAndWriteField(is%wrap%glue%grid, exportState, filePrefix="field_ocn_export_", &
    !  timeslice=is%wrap%slice, relaxedFlag=.true., rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    
    ! advance the time slice counter
    is%wrap%slice = is%wrap%slice + 1

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    type(InternalState)  :: is
    integer              :: stat

    rc = ESMF_SUCCESS
  
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! TODO: Destroy objects inside of internal state.

    ! Deallocate the internal state memory.
    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_RouteHandleRelease(CESM2HYCOM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_RouteHandleRelease(CESM2HYCOM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    deallocate(gsmap_o)
    deallocate(dom_o)
      
  end subroutine

  !-----------------------------------------------------------------------------
  subroutine routine_Run2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: compclock
    logical                       :: ocn_present
    logical                       :: ocnrun_alarm
    logical                       :: tight_coupling
    type(ESMF_Array)              :: d2x

    rc = ESMF_SUCCESS
    tight_coupling = .false.

    call ESMF_GridCompGet(gcomp, clock=compclock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(importState, name="ocean_tight_coupling", &
      value=tight_coupling, defaultvalue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocn_present", value=ocn_present, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocnrun_alarm", value=ocnrun_alarm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (ocn_present  .and.  ocnrun_alarm  .and.  tight_coupling) then
      call ESMF_LogWrite(trim('HYCOM RUN2 --->'), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call HYCOM_ModelAdvance(gcomp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine routine_Run3(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp

    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    logical                       :: ocn_present
    logical                       :: ocnrun_alarm
    logical                       :: tight_coupling
    type(ESMF_Array)              :: d2x
    real(ESMF_KIND_R8), pointer   :: rawdata(:,:)

    rc = ESMF_SUCCESS
    tight_coupling = .false.

    call ESMF_AttributeGet(importState, name="ocean_tight_coupling", &
      value=tight_coupling, defaultvalue=.false., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocn_present", value=ocn_present, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(exportState, &
      name="ocnrun_alarm", value=ocnrun_alarm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (ocn_present  .and.  ocnrun_alarm  .and.  (.not. tight_coupling)) then
      call ESMF_LogWrite(trim('HYCOM RUN3 --->'), ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      call HYCOM_ModelAdvance(gcomp, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    endif

  end subroutine
  !-----------------------------------------------------------------------------

  subroutine ocn_forcing(state, d2x, filename, rc)
    use shr_string_mod

    type(ESMF_State)                :: state
    type(ESMF_Array)                :: d2x
    character(len=*), intent(in)    :: filename
    integer, intent(out)            :: rc

    real(ESMF_KIND_R8), allocatable :: rawdata(:,:)
    integer                         :: gsize, nfields, elb(2,1), eub(2,1), lpet, rec_len
    type(ESMF_VM)                   :: vm

    rc = ESMF_SUCCESS

    call ESMF_AttributeGet(state, name="gsize", value=gsize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_ArrayGet(d2x, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    nfields = shr_string_listGetNum(trim(seq_flds_o2x_fields))
    !print *, 'ocn_forcing nfields: ', nfields, ' gsize: ', gsize

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_VMGet(vm, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) then
      allocate(rawdata(nfields, gsize))
      inquire (IOLENGTH=rec_len) rawdata
      open(1901,file=filename,status = 'unknown', form='unformatted', access='direct',recl=rec_len)
      read(1901,rec=1) rawdata
      close(1901)
    endif

    call ESMF_ArrayScatter(d2x, rawdata, 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) deallocate(rawdata)

  end subroutine

  subroutine HYCOM_WriteFieldBundle(grid, fieldbundle, fieldNameList, filePrefix, overwrite, &
    status, timeslice, relaxedflag, rc)
    type(ESMF_Grid),            intent(in)            :: grid
    type(ESMF_FieldBundle),     intent(in)            :: fieldbundle
    character(len=*),           intent(in),  optional :: fieldNameList(:)
    character(len=*),           intent(in),  optional :: filePrefix
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
  !-----------------------------------------------------------------------------
    ! local variables
    integer                         :: i, itemCount
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    character(len=80)               :: fileName
    character(len=80), allocatable  :: fieldNameList_loc(:)
    type(ESMF_Field)                :: dst2DField
    real(ESMF_KIND_R8), pointer     :: ptr1D(:), ptr2D(:,:)

    if (present(rc)) rc = ESMF_SUCCESS

    if (present(fieldNameList)) then
      allocate(fieldNameList_loc(size(fieldNameList)))
      do i=1, size(fieldNameList)
        fieldNameList_loc(i) = trim(fieldNameList(i))
      enddo
    else
      call ESMF_FieldBundleGet(fieldbundle, fieldCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList_loc(itemCount))
      call ESMF_FieldBundleGet(fieldbundle, fieldNameList=fieldNameList_loc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    do i=1, size(fieldNameList_loc)
      call ESMF_FieldBundleGet(fieldbundle, fieldName=fieldNameList_loc(i), field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out
      ! -> output to file
      if (present(filePrefix)) then
        write (fileName,"(A)") filePrefix//trim(fieldNameList_loc(i))//".nc"
      else
        write (fileName,"(A)") trim(fieldNameList_loc(i))//".nc"
      endif

      call ESMF_FieldGet(field, farrayPtr=ptr2D, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      
      !call ESMF_LogWrite(trim('HYCOM_WriteFieldBundle: '//fieldNameList_loc(i)), ESMF_LOGMSG_INFO, rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__, &
      !  file=__FILE__)) &
      !return ! bail out

!!Alex incompatibility ESMF 7.0.0      call ESMF_FieldWrite(field, file=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
!        overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
      call ESMF_FieldWrite(field, fileName=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
        overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out

    enddo

    deallocate(fieldNameList_loc)

  end subroutine

  subroutine RedistAndWriteField(grid, state, fieldNameList, filePrefix, overwrite, &
    status, timeslice, relaxedflag, rc)
    type(ESMF_Grid),            intent(in)            :: grid
    type(ESMF_State),           intent(in)            :: state
    character(len=*),           intent(in),  optional :: fieldNameList(:)
    character(len=*),           intent(in),  optional :: filePrefix
    logical,                    intent(in),  optional :: overwrite
    type(ESMF_FileStatus_Flag), intent(in),  optional :: status
    integer,                    intent(in),  optional :: timeslice
    logical,                    intent(in),  optional :: relaxedflag
    integer,                    intent(out), optional :: rc
  !-----------------------------------------------------------------------------
    ! local variables
    integer                         :: i, itemCount, elb(2), eub(2), j, k
    type(ESMF_Field)                :: field
    type(ESMF_StateItem_Flag)       :: itemType
    character(len=80)               :: fileName, msg
    character(len=80), allocatable  :: fieldNameList_loc(:)
    type(ESMF_Field)                :: dst2DField
    real(ESMF_KIND_R8), pointer     :: ptr1D(:), ptr2D(:,:)

    if (present(rc)) rc = ESMF_SUCCESS

    dst2DField = ESMF_FieldCreate(grid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if (present(fieldNameList)) then
      allocate(fieldNameList_loc(size(fieldNameList)))
      do i=1, size(fieldNameList)
        fieldNameList_loc(i) = trim(fieldNameList(i))
      enddo
    else
      call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(fieldNameList_loc(itemCount))
      call ESMF_StateGet(state, itemNameList=fieldNameList_loc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    do i=1, size(fieldNameList_loc)
      call ESMF_StateGet(state, itemName=fieldNameList_loc(i), &
        itemType=itemType, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=FILENAME)) &
        return  ! bail out
      if (itemType == ESMF_STATEITEM_FIELD) then
        ! field is available in the state
        call ESMF_StateGet(state, itemName=fieldNameList_loc(i), field=field, &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out
        ! -> output to file
        if (present(filePrefix)) then
          write (fileName,"(A)") filePrefix//trim(fieldNameList_loc(i))//".nc"
        else
          write (fileName,"(A)") trim(fieldNameList_loc(i))//".nc"
        endif

        call ESMF_FieldGet(field, farrayPtr=ptr1D, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
        return ! bail out

        call ESMF_FieldGet(dst2Dfield, farrayPtr=ptr2D, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
        return ! bail out
        ptr2D = 0.0_ESMF_KIND_R8

        write(msg, *) elb, eub, lbound(ptr1D), ubound(ptr1D)
        !call ESMF_LogWrite(trim('RedistAndWriteField: bounds: ')// trim(msg), ESMF_LOGMSG_INFO, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        !  line=__LINE__, &
        !  file=__FILE__)) &
        !return ! bail out

        !call ESMF_FieldRedist(field, dst2DField, routehandle=CESM2HYCOM_RHR8, rc=rc)
        call copy_1D_to_2D(field, dst2DField, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
        return ! bail out
        
        !call ESMF_LogWrite(trim('RedistAndWriteField: '//fieldNameList_loc(i)), ESMF_LOGMSG_INFO, rc=rc)
        !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        !  line=__LINE__, &
        !  file=__FILE__)) &
        !return ! bail out

!!Alex incompatibility ESMF 7.0.0        call ESMF_FieldWrite(dst2DField, file=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
!          overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)
        call ESMF_FieldWrite(dst2DField, fileName=trim(fileName), variableName=trim(fieldNameList_loc(i)), &
          overwrite=overwrite, status=status, timeslice=timeslice, rc=rc)

        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=FILENAME)) &
          return  ! bail out

      endif
    enddo

    deallocate(fieldNameList_loc)

    call ESMF_FieldDestroy(dst2DField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

  end subroutine

  subroutine copy_1D_to_2D(field1D, field2D, zeroDst, rc)
    use iso_c_binding

    type(ESMF_Field), intent(inout)    :: field1D, field2D
    logical, intent(in), optional      :: zeroDst
    integer, intent(out)               :: rc

    real(ESMF_KIND_R8), pointer        :: fptr1D(:), fptr2D(:,:), fptr2D_new(:,:)
    integer                            :: i,j,k
    logical                            :: l_zeroDst

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field1D, farrayPtr=fptr1D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_FieldGet(field2D, farrayPtr=fptr2D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    l_zeroDst = .false.
    if(present(zeroDst)) l_zeroDst = zeroDst
    if(l_zeroDst) then
      fptr2D = 0.
    else 
      k = 1
      do j = 1, ubound(fptr2D, 2)
        do i = 1, ubound(fptr2D, 1)
          fptr2D(i, j) = fptr1D(k)
          k = k + 1
        enddo
      enddo
    endif

    return
  end subroutine

  subroutine copy_2D_to_1D(field2D, field1D, rc)
    use iso_c_binding

    type(ESMF_Field), intent(inout)    :: field1D, field2D
    integer, intent(out)               :: rc

    real(ESMF_KIND_R8), pointer        :: fptr1D(:), fptr2D(:,:), fptr1D_new(:)
    integer                            :: i,j,k

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field1D, farrayPtr=fptr1D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_FieldGet(field2D, farrayPtr=fptr2D, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    k = 1
    do j = 1, ubound(fptr2D, 2)
      do i = 1, ubound(fptr2D, 1)
        fptr1D(k) = fptr2D(i,j)
        k = k + 1
      enddo
    enddo

    return
  end subroutine
    

  ! make the proper connections
  ! redist hycom 2d data to cesm 1d data
  subroutine HYCOM_RedistHYCOM2CESM(exportFields, hycom_field_stdname, &
    exportState, cesm_field_stdname, hycom_field_shortname, &
    initFlag, fptr, twolevel, connectOnly, zeroDst, rc)

    type(ESMF_FieldBundle)     , intent(in)                :: exportFields 
    character(len=*)           , intent(in)                :: hycom_field_stdname

    type(ESMF_State)           , intent(in)                :: exportState
    character(len=*)           , intent(in)                :: cesm_field_stdname

    character(len=*)           , intent(in),  optional     :: hycom_field_shortname

    logical                    , intent(in) , optional     :: initFlag
    real(ESMF_KIND_R8), pointer, intent(in) , optional     :: fptr(:,:)
    logical                    , intent(in) , optional     :: twolevel

    logical                    , intent(in) , optional     :: connectOnly
    logical                    , intent(in) , optional     :: zeroDst

    integer                    , intent(out), optional     :: rc

    !local
    type(ESMF_Field)                                       :: cesm_field  ! cesm field
    type(ESMF_Field)                                       :: hycom_field ! hycom field
    character(len=256)                                     :: cesm_field_shortname
    character(len=256)                                     :: l_hycom_field_shortname

    logical                                                :: l_connectOnly, l_zeroDst
    real(ESMF_KIND_R8), pointer                            :: fptr1D(:)

    rc = ESMF_SUCCESS

    !call ESMF_LogWrite(trim('HYCOM_RedistHYCOM2CESM: '// trim(hycom_field_stdname) // ' ---> ' //trim(cesm_field_stdname)), &
    !  ESMF_LOGMSG_INFO, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    ! retrieve 1D field from cesm export State
    call ESMF_StateGet(exportState, itemName=cesm_field_stdname, field =cesm_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(.not. present(hycom_field_shortname)) then
      call esmfshr_FieldDictionaryGetEntry(standardName=trim(hycom_field_stdname), &
        defaultShortName=l_hycom_field_shortname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else 
      l_hycom_field_shortname = hycom_field_shortname
    endif
    call ESMF_FieldBundleGet(exportFields, fieldname=l_hycom_field_shortname, field=hycom_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    l_zeroDst = .false.
    if(present(zeroDst)) l_zeroDst = zeroDst
    if(l_zeroDst) then
      call ESMF_FieldGet(cesm_field, farrayPtr=fptr1D, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      fptr1D = 0.
      return
    endif

    l_connectOnly = .false.
    if(present(connectOnly)) l_connectOnly = connectOnly

    call ESMF_AttributeSet(hycom_field, name="HYCOM_IN_CESM_Connected", value=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(.not. l_connectOnly) then
      !call ESMF_FieldRedist(hycom_field, cesm_field, routehandle=HYCOM2CESM_RHR8, rc=rc)
      call copy_2D_to_1D(hycom_field, cesm_field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
    endif

    !call ESMF_LogWrite(trim('HYCOM_RedistHYCOM2CESM: '// trim(hycom_field_stdname) // &
    !  ' : ' // trim(l_hycom_field_shortname) //' ---> ' //trim(cesm_field_stdname)//' : '//trim(cesm_field_shortname)), &
    !  ESMF_LOGMSG_INFO, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
    
  end subroutine

  ! redist cesm 1d data to hycom 2d data
  subroutine HYCOM_RedistCESM2HYCOM(importState, cesm_field_stdname, &
    importFields, hycom_field_stdname, hycom_field_shortname, &
    initFlag, fptr, twolevel, rc)

    type(ESMF_State)           , intent(in)                :: importState
    character(len=*)           , intent(in)                :: cesm_field_stdname

    type(ESMF_FieldBundle)     , intent(in)                :: importFields 
    character(len=*)           , intent(in)                :: hycom_field_stdname
    character(len=*)           , intent(in),  optional     :: hycom_field_shortname

    logical                    , intent(in),  optional     :: initFlag
    real(ESMF_KIND_R8), pointer, intent(in),  optional     :: fptr(:,:)
    logical                    , intent(in),  optional     :: twolevel

    integer                    , intent(out), optional     :: rc

    !local
    type(ESMF_Field)                                       :: cesm_field  ! cesm field
    type(ESMF_Field)                                       :: hycom_field ! hycom field
    character(len=256)                                     :: cesm_field_shortname
    character(len=256)                                     :: l_hycom_field_shortname

    rc = ESMF_SUCCESS

    !call ESMF_LogWrite(trim('HYCOM_RedistCESM2HYCOM: '// trim(cesm_field_stdname) // ' ---> ' //trim(hycom_field_stdname)), &
    !  ESMF_LOGMSG_INFO, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    ! retrieve 1D field from cesm import State
    call ESMF_StateGet(importState, itemName=cesm_field_stdname, field =cesm_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(.not. present(hycom_field_shortname)) then
      call esmfshr_FieldDictionaryGetEntry(standardName=trim(hycom_field_stdname), &
        defaultShortName=l_hycom_field_shortname, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      l_hycom_field_shortname = hycom_field_shortname
    endif
    call ESMF_FieldBundleGet(importFields, fieldname=l_hycom_field_shortname, field=hycom_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !call ESMF_FieldRedist(cesm_field, hycom_field, routehandle=CESM2HYCOM_RHR8, rc=rc)
    call copy_1D_to_2D(cesm_field, hycom_field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_AttributeSet(hycom_field, name="HYCOM_IN_CESM_Connected", value=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    !call ESMF_LogWrite(trim('HYCOM_RedistCESM2HYCOM: '// trim(cesm_field_stdname) // &
    !  ' : ' // trim(cesm_field_shortname) //' ---> ' //trim(hycom_field_stdname)//' : '//trim(l_hycom_field_shortname)), &
    !  ESMF_LOGMSG_INFO, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

  end subroutine

!***********************************************************************
!BOP
! !IROUTINE: ocn_domain_esmf
! !INTERFACE:

 subroutine ocn_domain_esmf( dom, grid, rc)

    use iso_c_binding

! !DESCRIPTION:
!  This routine creates the ocean domain and necessary communication routehandles
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT/OUTPUT PARAMETERS:

    implicit none
    type(ESMF_Array), intent(inout)     :: dom      ! CESM DOMAIN INFO
    type(ESMF_Grid),  intent(in)        :: grid     ! Native HYCOM 2D Grid
    integer, intent(out)                :: rc

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


    integer ::   &
      i,j, n, iblock

    integer ::   &
      klon,klat,karea,kmask,kfrac ! domain fields

    real(ESMF_KIND_R8),    pointer ::  &
      fptr (:,:)          ! data pointer into ESMF array

    real(ESMF_KIND_R8)  :: &
      frac                ! temporary var to compute frac/mask from KMT

    type(ESMF_DistGrid)  :: distgrid, distgrid2d
    type(ESMF_VM)        :: vm
    character(len=256)   :: msg
    integer              :: n_elem, n_pet, lpet, k
    integer, pointer     :: indexlist(:)
    logical              :: arbIndexFlag
    type(ESMF_Array)     :: lon1d, lat1d, area1d, mask1d
    type(ESMF_Array)     :: plon, plat, area, mask, area2d, mask2d
    integer              :: elb(2,1), eub(2,1), elb1(1,1), eub1(1,1)
    real(ESMF_KIND_R8), pointer  :: tlon(:), tlat(:), tarea(:), fptrLon(:,:), fptrLat(:,:)
    integer(ESMF_KIND_I4), pointer :: tmask(:), fptrSeqIndex(:)
    real(ESMF_KIND_R8)   :: radian, radius, pi
    type(ESMF_TYPEKIND_FLAG)  :: tkf

    type(ESMF_Array)     :: dummy1D, dummy2D
    real(ESMF_KIND_R8), pointer     :: fptr2D(:,:), fptr1D(:), fptr2D_new(:,:)

    type(ESMF_RouteHandle)          :: redist_padding_rh

!-----------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Retrieve domain data pointer
    call ESMF_ArrayGet(dom, localDe=0, farrayPtr=fptr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Retrieve the HYCOM Grid coordinate and mask arrays for reference      
    call ESMF_GridGetCoord(grid, coordDim=1, array=plon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, coordDim=2, array=plat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, array=mask, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(mask, typekind=tkf, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetItem(grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_AREA, array=area, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGet(grid, distgrid=distgrid2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    area2d = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    mask2d = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_I4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_VMGet(vm, petCount=n_pet, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Use the mesh based 1D distgrid to create DOM elements
    call ESMF_ArrayGet(dom, distgrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lon1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(lon1D, farrayPtr = tlon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    lat1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(lat1D, farrayPtr = tlat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    area1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(area1D, farrayPtr = tarea, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    mask1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_I4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(mask1D, farrayPtr = tmask, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! CESM uses 1 DE per PET
    call ESMF_DistGridGet(distgrid, 0, arbSeqIndexFlag=arbIndexFlag, elementCount=n_elem, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    allocate(indexList(n_elem))
    call ESMF_DistGridGet(distgrid, 0, seqIndexList=indexlist, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

!-------------------------------------------------------------------
!
!  initialize domain type, lat/lon in degrees,
!  area in radians^2, mask is 1 (ocean), 0 (non-ocean)
!  Fill in correct values for domain components
!
!-------------------------------------------------------------------

    klon  = esmfshr_util_ArrayGetIndex(dom,'lon ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    klat  = esmfshr_util_ArrayGetIndex(dom,'lat ',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    karea = esmfshr_util_ArrayGetIndex(dom,'area',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kmask = esmfshr_util_ArrayGetIndex(dom,'mask',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    kfrac = esmfshr_util_ArrayGetIndex(dom,'frac',rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    fptr(:,:) = -9999.0_ESMF_KIND_R8
    n=0

    write(msg, *) 'DUMPING HYCOM INDICES BEGINS:'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'total number of ocean pet', n_pet, ' local pet number', lpet
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'number of elements on this pet:', n_elem, ' arbflag', arbIndexFlag
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'lpet ', 'n ', 'Index ', 'lon ', 'lat ', &
      'area ', 'frac ', 'mask'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(plon, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(msg, *) 'src shape: ', elb, eub, ' dst shape: ', lbound(fptr), ubound(fptr)
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    write(msg, *) 'plon: ', elb, eub
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(lon1d, exclusiveLBound=elb1, exclusiveUBound=eub1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write(msg, *) 'lon1d: ', elb1, eub1
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayRedistStore(plon, lon1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(plon, lon1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(plat, lat1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(mask, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedist(area, area1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    pi = 3.14159265358
    radian = 180.0_ESMF_KIND_R8/pi
    radius    = 6370.0e3_ESMF_KIND_R8
    do n = elb1(1,1), eub1(1,1)
      fptr(klon , n)          = TLON(n)
      fptr(klat , n)          = TLAT(n)
      fptr(karea, n)          = TAREA(n)/radius/radius
      frac                    = TMASK(n)
      if (frac > 1.0_ESMF_KIND_R8) frac = 1.0_ESMF_KIND_R8
      fptr(kfrac, n)          = frac
      fptr(kmask, n)          = frac
      !write(msg, '(I4,A1,I8,A7,I8,2F10.3,E15.7,2F10.3)') lpet, ' ', n, ' INDEX=',indexlist(n), fptr(klon, n), &
      !  fptr(klat , n), fptr(karea, n), fptr(kfrac, n), fptr(kmask, n)
      !call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
      !if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    enddo

    write(msg, *) 'DUMPING HYCOM INDICES ENDS:'
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_ArrayGet(plon, farrayPtr=fptrLon, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_ArrayGet(plat, farrayPtr=fptrLat, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_DistGridGet(distgrid2D, 0, arbSeqIndexFlag=arbIndexFlag, elementCount=n_elem, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    allocate(fptrSeqIndex(n_elem))
    call ESMF_DistGridGet(distgrid2D, 0, seqIndexList=fptrSeqIndex, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    write(msg, *) 'HYCOM 2D distribution is arbitrary? ', arbIndexFlag, n_elem
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    !do i = elb(1,1), eub(1,1)
    !  do j = elb(2,1), eub(2,1) 
    !    write(msg, '(2I6,2F10.3)') i, j, fptrLon(i,j), fptrLat(i,j)
    !    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    !    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    !  enddo
    !enddo

    !do n = 1, n_elem
    !  write(msg, '(I6,I8)') n, fptrSeqIndex(n)
    !  call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    !  if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    !enddo

    ! Release these routehandles because they have paddings(halo) in hycom memory
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_RouteHandleRelease(HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    ! Recompute routehandles that has no paddings. 

    call ESMF_ArrayRedistStore(area1d, area2d, routehandle=CESM2HYCOM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask1D, mask2d, routehandle=CESM2HYCOM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(area2d, area1d, routehandle=HYCOM2CESM_RHR8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayRedistStore(mask2D, mask1d, routehandle=HYCOM2CESM_RHI4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !call ESMF_ArrayRedist(area1d, area2d, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !dummy1D = ESMF_ArrayCreate(distgrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out
    !dummy2D = ESMF_ArrayCreate(distgrid2d, typekind=ESMF_TYPEKIND_R8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayGet(dummy2D, farrayPtr=fptr2D, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !fptr2D(:,:) = lpet

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=1, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy2D, dummy1D, routehandle=HYCOM2CESM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy1D, dummy2D, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=2, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !k = 1
    !do j = 1, ubound(fptr2D, 2)
    !  do i = 1, ubound(fptr2D, 1)
    !    fptr2D(i, j) = fptrSeqIndex(k)
    !    k = k + 1
    !  enddo
    !enddo

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=3, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &

    !call ESMF_ArrayRedist(dummy2D, dummy1D, routehandle=HYCOM2CESM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayRedist(dummy1D, dummy2D, routehandle=CESM2HYCOM_RHR8, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    !call ESMF_ArrayWrite(dummy2D, file='scramble.nc', variableName='dummy2D', &
    !  overwrite=.true., timeslice=4, rc=rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !return ! bail out

    deallocate(indexlist)
    deallocate(fptrSeqIndex)

    return

!-----------------------------------------------------------------------
!EOC

  end subroutine ocn_domain_esmf

  subroutine initialize_HYCOM_CESM_tables(rc)

    integer, intent(out) :: rc

    integer              :: n_entries, i

    rc = ESMF_SUCCESS

    do i = 1, 30
      cesm2hycom_table(i)%connected = .true.
      cesm2hycom_table(i)%unit = ''
    enddo
    do i = 1, 8
      hycom2cesm_table(i)%connected = .true.
      hycom2cesm_table(i)%unit = ''
    enddo

    cesm2hycom_table(1)%hycom_stdname  = "surface_downward_eastward_stress"
    cesm2hycom_table(1)%cesm_stdname   = "Foxx_taux"
    cesm2hycom_table(1)%unit           = "Pa"

    cesm2hycom_table(2)%hycom_stdname  = "surface_downward_northward_stress"
    cesm2hycom_table(2)%cesm_stdname   = "Foxx_tauy"
    cesm2hycom_table(2)%unit           = "Pa"

    cesm2hycom_table(3)%hycom_stdname  = "wind_speed_height10m"
    cesm2hycom_table(3)%cesm_stdname   = "Sx_u10"
    cesm2hycom_table(3)%unit           = "m/s"
    cesm2hycom_table(3)%connected      = .false.

    cesm2hycom_table(4)%hycom_stdname  = "friction_speed"
    cesm2hycom_table(4)%cesm_stdname   = "So_ustar"
    cesm2hycom_table(4)%unit           = "m/s"
    cesm2hycom_table(4)%connected      = .false.

    cesm2hycom_table(5)%hycom_stdname  = "mean_down_sw_flx"                      
    cesm2hycom_table(5)%cesm_stdname   = ""
    cesm2hycom_table(5)%unit           = "W/m^2"
    cesm2hycom_table(5)%connected      = .false.

    cesm2hycom_table(6)%hycom_stdname  = "mean_net_sw_flx"
    cesm2hycom_table(6)%cesm_stdname   = "Foxx_swnet"
    cesm2hycom_table(6)%unit           = "W/m^2"

    cesm2hycom_table(7)%hycom_stdname  = "mean_net_lw_flx"
    cesm2hycom_table(7)%cesm_stdname   = ""
    cesm2hycom_table(7)%unit           = "W/m^2"
    cesm2hycom_table(7)%connected      = .false.

    cesm2hycom_table(8)%hycom_stdname  = "mean_down_lw_flx"
    cesm2hycom_table(8)%cesm_stdname   = "Faxa_lwdn"
    cesm2hycom_table(8)%unit           = "W/m^2"

    cesm2hycom_table(9)%hycom_stdname  = "mean_up_lw_flx"
    cesm2hycom_table(9)%cesm_stdname   = "Foxx_lwup"
    cesm2hycom_table(9)%unit           = "W/m^2"

    cesm2hycom_table(10)%hycom_stdname  = "inst_temp_height2m"
    cesm2hycom_table(10)%cesm_stdname   = "Sa_tbot"
    cesm2hycom_table(10)%unit           = "K"
    cesm2hycom_table(10)%connected      = .false.

    cesm2hycom_table(11)%hycom_stdname = "mean_prec_rate"
    cesm2hycom_table(11)%cesm_stdname  = "Faxa_prec"
    cesm2hycom_table(11)%unit          = "Kg/m^2/s"

    cesm2hycom_table(12)%hycom_stdname = "inst_spec_humid_height2m"
    cesm2hycom_table(12)%cesm_stdname  = "Sa_shum"
    cesm2hycom_table(12)%unit          = "Kg/Kg"
    cesm2hycom_table(12)%connected     = .false.

    cesm2hycom_table(13)%hycom_stdname = "sea_surface_temperature"
    cesm2hycom_table(13)%cesm_stdname  = "Sa_ptem"
    cesm2hycom_table(13)%unit          = "K"
    cesm2hycom_table(13)%connected     = .false.

    cesm2hycom_table(14)%hycom_stdname = "sea_ice_area_fraction"
    cesm2hycom_table(14)%cesm_stdname  = "Si_ifrac"
    cesm2hycom_table(14)%unit          = "1"

    cesm2hycom_table(15)%hycom_stdname = "downward_x_stress_at_sea_ice_base"
    cesm2hycom_table(15)%cesm_stdname  = "Fioi_taux"
    cesm2hycom_table(15)%unit          = "Pa"
    cesm2hycom_table(15)%connected     = .false.

    cesm2hycom_table(16)%hycom_stdname = "downward_y_stress_at_sea_ice_base"
    cesm2hycom_table(16)%cesm_stdname  = "Fioi_tauy"
    cesm2hycom_table(16)%unit          = "Pa"
    cesm2hycom_table(16)%connected     = .false.

    cesm2hycom_table(17)%hycom_stdname = "downward_sea_ice_basal_solar_heat_flux"
    cesm2hycom_table(17)%cesm_stdname  = ""
    cesm2hycom_table(17)%unit          = "W/m^2"
    cesm2hycom_table(17)%connected     = .false.

    cesm2hycom_table(18)%hycom_stdname = "upward_sea_ice_basal_heat_flux"
    cesm2hycom_table(18)%cesm_stdname  = "Fioi_melth"
    cesm2hycom_table(18)%unit          = "W/^2"

    cesm2hycom_table(19)%hycom_stdname = "downward_sea_ice_basal_salt_flux"
    cesm2hycom_table(19)%cesm_stdname  = "Fioi_salt"
    cesm2hycom_table(19)%unit          = "Kg/m^2/s"

    cesm2hycom_table(20)%hycom_stdname = "downward_sea_ice_basal_water_flux"
    cesm2hycom_table(20)%cesm_stdname  = "Fioi_meltw"
    cesm2hycom_table(20)%unit          = "Kg/m^2/s"

    cesm2hycom_table(21)%hycom_stdname = "sea_ice_temperature"
    cesm2hycom_table(21)%cesm_stdname  = ""
    cesm2hycom_table(21)%unit          = "K"
    cesm2hycom_table(21)%connected     = .false.

    cesm2hycom_table(22)%hycom_stdname = "sea_ice_thickness"
    cesm2hycom_table(22)%cesm_stdname  = ""
    cesm2hycom_table(22)%unit          = "m"
    cesm2hycom_table(22)%connected     = .false.

    cesm2hycom_table(23)%hycom_stdname = "sea_ice_x_velocity"
    cesm2hycom_table(23)%cesm_stdname  = ""
    cesm2hycom_table(23)%unit          = "m/s"
    cesm2hycom_table(23)%connected     = .false.

    cesm2hycom_table(24)%hycom_stdname = "sea_ice_y_velocity"
    cesm2hycom_table(24)%cesm_stdname  = ""
    cesm2hycom_table(24)%unit          = "m/s"
    cesm2hycom_table(24)%connected     = .false.

    cesm2hycom_table(25)%hycom_stdname = "downward_x_stress_ocean"
    cesm2hycom_table(25)%cesm_stdname  = "Faox_taux"
    cesm2hycom_table(25)%unit          = "Pa"
    cesm2hycom_table(25)%connected     = .false.

    cesm2hycom_table(26)%hycom_stdname = "downward_y_stress_ocean"
    cesm2hycom_table(26)%cesm_stdname  = "Faox_tauy"
    cesm2hycom_table(26)%unit          = "Pa"
    cesm2hycom_table(26)%connected     = .false.

    cesm2hycom_table(27)%hycom_stdname = "mean_lat_flx"
    cesm2hycom_table(27)%cesm_stdname  = "Foxx_lat"
    cesm2hycom_table(27)%unit          = "W/m^2"

    cesm2hycom_table(28)%hycom_stdname = "mean_sens_flx"
    cesm2hycom_table(28)%cesm_stdname  = "Foxx_sen"
    cesm2hycom_table(28)%unit          = "W/m^2"

    cesm2hycom_table(29)%hycom_stdname = "water_flux_into_sea_water"
    cesm2hycom_table(29)%cesm_stdname  = "Foxx_rofl"
    cesm2hycom_table(29)%unit          = "Kg/m^2/s"

    cesm2hycom_table(30)%hycom_stdname = "frozen_water_flux_into_sea_water"
    cesm2hycom_table(30)%cesm_stdname  = "Foxx_rofi"
    cesm2hycom_table(30)%unit          = "Kg/m^2/s"

    ! --------------------------------------------------------------------------
    hycom2cesm_table(1)%hycom_stdname = "sea_surface_temperature"
    hycom2cesm_table(1)%cesm_stdname  = "So_t"
    hycom2cesm_table(1)%unit          = "K"

    hycom2cesm_table(2)%hycom_stdname = "s_surf"
    hycom2cesm_table(2)%cesm_stdname  = "So_s"
    hycom2cesm_table(2)%unit          = "g/Kg"

    hycom2cesm_table(3)%hycom_stdname = "ocn_current_zonal"
    hycom2cesm_table(3)%cesm_stdname  = "So_u"
    hycom2cesm_table(3)%unit          = "m/s"

    hycom2cesm_table(4)%hycom_stdname = "ocn_current_merid"
    hycom2cesm_table(4)%cesm_stdname  = "So_v"
    hycom2cesm_table(4)%unit          = "m/s"

    hycom2cesm_table(5)%hycom_stdname = "eastward_sea_surface_slope"
    hycom2cesm_table(5)%cesm_stdname  = "So_dhdx"
    hycom2cesm_table(5)%unit          = ""

    hycom2cesm_table(6)%hycom_stdname = "northward_sea_surface_slope"
    hycom2cesm_table(6)%cesm_stdname  = "So_dhdy"
    hycom2cesm_table(6)%unit          = ""

    hycom2cesm_table(7)%hycom_stdname = "upward_sea_ice_basal_available_heat_flux"
    hycom2cesm_table(7)%cesm_stdname  = "Fioo_q"
    hycom2cesm_table(7)%unit          = "W/m^2"

    hycom2cesm_table(8)%hycom_stdname = "mixed_layer_depth"
    hycom2cesm_table(8)%cesm_stdname  = "So_bldepth"
    hycom2cesm_table(8)%unit          = "m"

  end subroutine

  subroutine dumpRawData(state, filename, distArrayName, fieldlist, rc)

    use shr_string_mod

    type(ESMF_State), intent(in)    :: state
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: distArrayName
    character(len=*), intent(in)    :: fieldlist
    integer, intent(out)            :: rc

    type(ESMF_Array)                :: d2x
    real(ESMF_KIND_R8), allocatable :: rawdata(:,:)
    integer                         :: gsize, nfields, elb(2,1), eub(2,1), lpet, rec_len
    type(ESMF_VM)                   :: vm

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemName=distArrayName, array=d2x, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
    call ESMF_AttributeGet(state, name="gsize", value=gsize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_ArrayGet(d2x, exclusiveLBound=elb, exclusiveUBound=eub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    nfields = shr_string_listGetNum(trim(fieldlist))
    !nfields = eub(2,1)
!!Alex    print *, 'dumpRawData: nfields: ', nfields

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    call ESMF_VMGet(vm, localPet=lpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) allocate(rawdata(nfields, gsize))

    call ESMF_ArrayGather(d2x, rawdata, 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out

    if(lpet == 0) then
      inquire (IOLENGTH=rec_len) rawdata
      open(1901,file=filename,status = 'unknown', form='unformatted', access='direct',recl=rec_len)
      write(1901,rec=1) rawdata
      close(1901)
    endif

    if(lpet == 0) deallocate(rawdata)

  end subroutine

  ! compute regridding weights for all regridding pairs
  subroutine HYCOM_CESM_REGRID_TOOCN(srcGridFile, dstGrid, weightFile, regridMethod, rc)
    character(len=*), intent(in)             :: srcGridFile
    type(ESMF_Grid), intent(in)              :: dstGrid
    character(len=*), intent(in)             :: weightFile
    type(ESMF_REGRIDMETHOD_FLAG), intent(in) :: regridMethod
    integer, intent(out)                     :: rc

    ! local
    character(len=125)                       :: path="/glade/p/work/feiliu/weights/T62/"

    type(ESMF_Grid)                          :: srcGrid
    type(ESMF_Field)                         :: srcField, dstField
    type(ESMF_Field)                         :: srcFracField, dstFracField
    integer, pointer                         :: factorIndexList(:,:)
    real(ESMF_KIND_R8), pointer              :: factorList(:)

    rc = ESMF_SUCCESS

    srcGrid = ESMF_GridCreate(srcGridFile, ESMF_FILEFORMAT_SCRIP, (/5,8/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    srcField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dstField = ESMF_FieldCreate(dstGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    srcFracField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dstFracField = ESMF_FieldCreate(dstGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldRegridStore(srcField, dstField, regridmethod=regridMethod, &
      !srcFracField=srcFracField, dstFracField=dstFracField, &
      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
      factorList=factorList, factorIndexList=factorIndexList, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_OutputScripWeightFile(trim(path)//trim(weightFile), factorList, factorIndexList, &
      method=ESMF_REGRIDMETHOD_BILINEAR, &
      srcFile=srcGridFile, dstFile="/glade/p/cesm/cseg/mapping/grids/gx1v6_090205.nc", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  ! compute regridding weights for all regridding pairs
  subroutine HYCOM_CESM_REGRID_FROMOCN(srcGrid, dstGridFile, weightFile, regridMethod, rc)
    type(ESMF_Grid), intent(in)              :: srcGrid
    character(len=*), intent(in)             :: dstGridFile
    character(len=*), intent(in)             :: weightFile
    type(ESMF_REGRIDMETHOD_FLAG), intent(in) :: regridMethod
    integer, intent(out)                     :: rc

    ! local
    character(len=125)                       :: path="/glade/p/work/feiliu/weights/T62/"

    type(ESMF_Grid)                          :: dstGrid
    type(ESMF_Field)                         :: srcField, dstField
    type(ESMF_Field)                         :: srcFracField, dstFracField
    integer, pointer                         :: factorIndexList(:,:)
    real(ESMF_KIND_R8), pointer              :: factorList(:)

    rc = ESMF_SUCCESS

    dstGrid = ESMF_GridCreate(dstGridFile, ESMF_FILEFORMAT_SCRIP, (/5,8/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    srcField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dstField = ESMF_FieldCreate(dstGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    srcFracField = ESMF_FieldCreate(srcGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    dstFracField = ESMF_FieldCreate(dstGrid, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldRegridStore(srcField, dstField, regridmethod=regridMethod, &
      ! enable for conservative regridding
      !srcFracField=srcFracField, dstFracField=dstFracField, &
      unmappedaction=ESMF_UNMAPPEDACTION_IGNORE, &
      factorList=factorList, factorIndexList=factorIndexList, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_OutputScripWeightFile(trim(path)//trim(weightFile), factorList, factorIndexList, &
      method=ESMF_REGRIDMETHOD_BILINEAR, &
      srcFile="/glade/p/cesm/cseg/mapping/grids/gx1v6_090205.nc", dstFile=dstGridFile, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  ! All standard names here are defined by HYCOM, not CESM dependent
  ! A translation happens with hycom2cesm and cesm2hycom tables
  subroutine loadHycomDictionary(rc)
    integer, intent(out)                     :: rc
    rc = ESMF_SUCCESS

    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "surface_downward_eastward_stress")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="surface_downward_eastward_stress", &
        canonicalUnits="Pa", &
        defaultLongName="Surface Downward Eastward Stress", &
        defaultShortName="taux", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "surface_downward_northward_stress")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="surface_downward_northward_stress", &
        canonicalUnits="Pa", &
        defaultLongName="Surface Downward Northward Stress", &
        defaultShortName="tauy", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "wind_speed_height10m")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="wind_speed_height10m", &
        canonicalUnits="m s-1", &
        defaultLongName="Wind Speed at 10m", &
        defaultShortName="wndspd", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "friction_speed")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="friction_speed", &
        canonicalUnits="m s-1", &
        defaultLongName="Friction Speed", &
        defaultShortName="ustara", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "mean_down_sw_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_down_sw_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Downward Short Wave Radiation Flux", &
        defaultShortName="mdswfx", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "mean_net_sw_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_net_sw_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Net Short Wave Radiation Flux", &
        defaultShortName="mnswfx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "mean_down_lw_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_down_lw_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Downward Long Wave Radiation Flux", &
        defaultShortName="mdlwfx", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "mean_up_lw_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_up_lw_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Upward Long Wave Radiation Flux", &
        defaultShortName="mulwfx", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "mean_net_lw_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_net_lw_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Net Long Wave Radiation Flux", &
        defaultShortName="mnlwfx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "mean_lat_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_lat_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Latent Heat Flux", &
        defaultShortName="mnlatfx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "mean_sens_flx")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_sens_flx", &
        canonicalUnits="W m-2", &
        defaultLongName="Mean Sensible Heat Flux", &
        defaultShortName="mnsenfx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "inst_temp_height2m")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="inst_temp_height2m", &
        canonicalUnits="K", &
        defaultLongName="Instantaneous Temperature 2m Above Ground", &
        defaultShortName="ith2m", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "inst_spec_humid_height2m")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="inst_spec_humid_height2m", &
        canonicalUnits="kg kg-1", &
        defaultLongName="Instantaneous Specific Humidity 2m Above Ground", &
        defaultShortName="ishh2m", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not.NUOPC_FieldDictionaryHasEntry( &
      "mean_prec_rate")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mean_prec_rate", &
        canonicalUnits="kg s m-2", &
        defaultLongName="Mean Liquid Precipitation Rate", &
        defaultShortName="lprec", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "air_surface_temperature")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="air_surface_temperature", &
        canonicalUnits="K", &
        defaultLongName="Air surface temerature", &
        defaultShortName="surtmp", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "upward_sea_ice_basal_available_heat_flux")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="upward_sea_ice_basal_available_heat_flux", &
        canonicalUnits="W m-2", &
        defaultLongName="Oceanic Heat Flux Available to Sea Ice", &
        defaultShortName="ssfi", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_lev")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_lev", &
        canonicalUnits="m", &
        defaultLongName="sea level", &
        defaultShortName="sea_lev", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "mixed_layer_depth")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="mixed_layer_depth", &
        canonicalUnits="m", &
        defaultLongName="Mixed Layer Depth in Ocean", &
        defaultShortName="mld", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "s_surf")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="s_surf", &
        canonicalUnits="psu", &
        defaultLongName="sea surface salinity on t-cell", &
        defaultShortName="sss", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_ice_area_fraction")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_ice_area_fraction", &
        canonicalUnits="1", &
        defaultLongName="Sea Ice Concentration", &
        defaultShortName="sic", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_x_stress_at_sea_ice_base")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_x_stress_at_sea_ice_base", &
        canonicalUnits="Pa", &
        defaultLongName="Sea Ice X-Stress", &
        defaultShortName="sitx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_y_stress_at_sea_ice_base")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_y_stress_at_sea_ice_base", &
        canonicalUnits="Pa", &
        defaultLongName="Sea Ice Y-Stress", &
        defaultShortName="sity", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_sea_ice_basal_solar_heat_flux")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_sea_ice_basal_solar_heat_flux", &
        canonicalUnits="W m-2", &
        defaultLongName="Solar Heat Flux thru Ice to Ocean", &
        defaultShortName="siqs", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "upward_sea_ice_basal_heat_flux")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="upward_sea_ice_basal_heat_flux", &
        canonicalUnits="W m-2", &
        defaultLongName="Ice Freezing/Melting Heat Flux", &
        defaultShortName="sifh", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_sea_ice_basal_salt_flux")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_sea_ice_basal_salt_flux", &
        canonicalUnits="kg m-2 s-1", &
        defaultLongName="Ice Freezing/Melting Salt Flux", &
        defaultShortName="sifs", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_sea_ice_basal_water_flux")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_sea_ice_basal_water_flux", &
        canonicalUnits="kg m-2 s-1", &
        defaultLongName="Ice Net Water Flux", &
        defaultShortName="sifw", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_ice_temperature")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_ice_temperature", &
        canonicalUnits="K", &
        defaultLongName="Sea Ice Temperature", &
        defaultShortName="sit", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_ice_thickness")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_ice_thickness", &
        canonicalUnits="m", &
        defaultLongName="Sea Ice Thickness", &
        defaultShortName="sih", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_ice_x_velocity")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_ice_x_velocity", &
        canonicalUnits="m s-1", &
        defaultLongName="Sea Ice X-Velocity", &
        defaultShortName="siu", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "sea_ice_y_velocity")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="sea_ice_y_velocity", &
        canonicalUnits="m s-1", &
        defaultLongName="Sea Ice Y-Velocity", &
        defaultShortName="siv", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "dummyfield")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="dummyfield", &
        canonicalUnits="1", &
        defaultLongName="Dummy Test Field", &
        defaultShortName="dummyfield", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "ocn_current_zonal")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="ocn_current_zonal", &
        canonicalUnits="m s-1", &
        defaultLongName="ocean current zonal component", &
        defaultShortName="ocncz", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "ocn_current_merid")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="ocn_current_merid", &
        canonicalUnits="m s-1", &
        defaultLongName="ocean current meridional component", &
        defaultShortName="ocncm", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_x_stress_ocean")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_x_stress_ocean", &
        canonicalUnits="Pa", &
        defaultLongName="ocean downward eastward stress", &
        defaultShortName="sotx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "downward_y_stress_ocean")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="downward_y_stress_ocean", &
        canonicalUnits="Pa", &
        defaultLongName="ocean downward northward stress", &
        defaultShortName="soty", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "eastward_sea_surface_slope")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="eastward_sea_surface_slope", &
        canonicalUnits="", &
        defaultLongName="eastward sea surface slope", &
        defaultShortName="dhdx", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "northward_sea_surface_slope")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="northward_sea_surface_slope", &
        canonicalUnits="", &
        defaultLongName="northward sea surface slope", &
        defaultShortName="dhdy", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 

    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "water_flux_into_sea_water")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="water_flux_into_sea_water", &
        canonicalUnits="kg m-2 s-1", &
        defaultLongName="Water flux due to runoff (liquid)", &
        defaultShortName="rofl", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 
    if (.not. NUOPC_FieldDictionaryHasEntry( &
      "frozen_water_flux_into_sea_water")) then
      call esmfshr_FieldDictionaryAddEntry( &
        standardName="frozen_water_flux_into_sea_water", &
        canonicalUnits="kg m-2 s-1", &
        defaultLongName="Water flux due to runoff (frozen)", &
        defaultShortName="rofi", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif 

  end subroutine

  
end module
