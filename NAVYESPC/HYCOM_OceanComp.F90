! ESMF macros for logging
#define FILENAME "HYCOM_Ocean_Comp.F90"
!#define CONTEXT  __LINE__,FILENAME,METHOD
#define CONTEXT  line=__LINE__,file=__FILE__

! Define ESMF real kind to match HYCOM single/double precision
#if defined(ESPC_IMPEXP_SINGLE)
#define ESMF_KIND_RX ESMF_KIND_R4
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R4
#else
#define ESMF_KIND_RX ESMF_KIND_R8
#define ESMF_TYPEKIND_RX ESMF_TYPEKIND_R8
#endif


!=========================================================================================
! HYCOM OCEAN ESMF Module
!=========================================================================================
#if defined ESPC_OCN
  MODULE OCEAN_Mod
#else
  MODULE HYCOM_Mod   
#endif
 
! ESMF Framework module
  USE ESMF

  use NUOPC

  use NUOPC_Model, only: &
    model_routine_SS    => SetServices, &
    model_label_Advance => label_Advance, &
    model_label_Finalize => label_Finalize

! HYCOM OCEAN forecast module
  use mod_hycom, only : end_of_run_cpl, &
                        end_of_run, &
                        HYCOM_Init, &
                        HYCOM_Run,  &
                        HYCOM_Final
  use mod_archiv, only: archiv_exchange
  use hycom_couple, only: idim_size,jdim_size,deBList, &
    lon_e,lat_e,mask_e,lon_q,lat_q, &
    hycom_couple_init,hycom_couple_final,set_hycom_import_flag, &
    import_to_hycom_deb, export_from_hycom_deb,ocn_import_forcing

# ifdef ESPC_COUPLE
  use read_impexp_config_mod
  use impexpField_cdf_mod
# endif

  IMPLICIT NONE
  PRIVATE
  SAVE
    
! Public member functions
  PUBLIC OCEAN_SetServices

! Miscellaneous
  TYPE(ESMF_VM) :: vm
!  INTEGER :: mpiCommunicator
  INTEGER :: nPets, lPet
  INTEGER, PARAMETER :: localDE=0
  TYPE(ESMF_Clock) :: intClock
  TYPE(ESMF_Config) :: config
  real ocean_start_dtg, ocean_end_dtg

# ifdef ESPC_COUPLE

  INTEGER :: numImpFields, numExpFields
  type(ESMF_Field), dimension(:), allocatable  :: impField
  type(ESMF_Field), dimension(:), allocatable  :: expField

  character(len=30), pointer :: expFieldName(:),impFieldName(:) => NULL()
  character(len=60), pointer :: expStandName(:),impStandName(:) => NULL()
  character(len=30), pointer :: expFieldUnit(:),impFieldUnit(:) => NULL()
  logical, pointer :: expFieldEnable(:),impFieldEnable(:) => NULL()


  integer cdf_impexp_freq
  real cpl_time_step
  logical ocn_esmf_exp_output, ocn_esmf_imp_output, hycom_arche_output
!!Alex  character (len=10) :: base_dtg
  character (len=15) :: base_dtg

# endif

  real endtime

  integer itdmx,jtdmx

  logical show_minmax

# ifdef ESPC_TIMER
  real(kind=ESMF_KIND_R8) :: timer_beg,timer_end

! espc_timer(1): Init Phase
! espc_timer(2): Run Phase
! espc_timer(3): Final Phase
! espc_timer(4): Run Phase import
! espc_timer(5): Run Phase Core
! espc_timer(6): Run Phase export
  real(kind=ESMF_KIND_R8) :: espc_timer(6)
# endif


!=========================================================================================
! HYCOM OCEAN ESMF Module Subroutines
!=========================================================================================
CONTAINS


#undef METHOD
#define METHOD "OCEAN_SetServices"
SUBROUTINE OCEAN_SetServices(model, rc)

! Calling parameters
  TYPE(ESMF_GridComp) :: model  
  INTEGER,INTENT(OUT) :: rc

! Assume failure
  rc = ESMF_FAILURE


    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN


    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN


    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
       specRoutine=OCEAN_Final, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="set final entry point failed", CONTEXT)) RETURN

!   Return success
    rc = ESMF_SUCCESS

END SUBROUTINE OCEAN_SetServices


!======================================================================
#undef METHOD
#define METHOD "InitializeP1"
SUBROUTINE InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model 
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer i,irc

# ifdef ESPC_COUPLE
    character(len=30) :: ocn_impexp_file
# endif

    INTEGER :: cpl_hour, cpl_min, cpl_sec

    rc = ESMF_SUCCESS


!   Get VM and Config info
    CALL ESMF_GridCompGet(model,config=config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="GridCompGet failed", CONTEXT)) RETURN

    CALL ESMF_VMGet(vm, petCount=nPets,localPet=lPet,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="VMGet failed", CONTEXT)) RETURN

  if(lPet.eq.0) print *,"hycom, InitializeP1 called,nPets=",nPets

# ifdef ESPC_TIMER
    do i=1,6
      espc_timer(i)=0.
    enddo

    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif



# ifdef ESPC_COUPLE


    CALL ESMF_ConfigGetAttribute(config, cdf_impexp_freq, &
       label=TRIM("cdf_impexp_freq="), default=9999, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get cdf_impexp_freq failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, cpl_hour, &
       label=TRIM("cpl_hour="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get cpl_hour failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, cpl_min, &
       label=TRIM("cpl_min="), default=0, rc=rc)

    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get cpl_min failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, cpl_sec, &
       label=TRIM("cpl_sec="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get cpl_sec failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    cpl_time_step=(cpl_hour+cpl_min/60.+cpl_sec/3600.)

    CALL ESMF_ConfigGetAttribute(config, base_dtg, &
       label=TRIM("base_dtg="), default='9999999999', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get base_dtg failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
    
    CALL ESMF_ConfigGetAttribute(config, hycom_arche_output, &
       label=TRIM("hycom_arche_output="), default=.true., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get hycom_arche_output failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
   
#if defined ESPC_OCN 
    CALL ESMF_ConfigGetAttribute(config, ocn_esmf_exp_output, &
       label=TRIM("ocn_esmf_exp_output="), default=.true., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get ocn_esmf_output failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
#else
    CALL ESMF_ConfigGetAttribute(config, ocn_esmf_exp_output, &
       label=TRIM("hyc_esmf_exp_output="), default=.true., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get hyc_esmf_output failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
#endif   

#if defined ESPC_OCN 
    CALL ESMF_ConfigGetAttribute(config, ocn_esmf_imp_output, &
       label=TRIM("ocn_esmf_imp_output="), default=.true., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get ocn_esmf_output failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
#else
    CALL ESMF_ConfigGetAttribute(config, ocn_esmf_imp_output, &
       label=TRIM("hyc_esmf_imp_output="), default=.true., rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get hyc_esmf_output failed", &
      CONTEXT,  rcToReturn=rc)) RETURN
#endif

#if defined ESPC_OCN
    CALL ESMF_ConfigGetAttribute(config, ocn_impexp_file, &
       label=TRIM("ocn_impexp_file="), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get ocn_impexp_file failed", &
     CONTEXT, rcToReturn=rc)) RETURN
#else
    CALL ESMF_ConfigGetAttribute(config, ocn_impexp_file, &
       label=TRIM("hyc_impexp_file="), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get hyc_impexp_file failed", &
     CONTEXT, rcToReturn=rc)) RETURN
#endif
    ocn_impexp_file=trim(ocn_impexp_file)
    if(lPet.eq.0) print *,"ocn_impexp_file=",ocn_impexp_file

    call read_impexp_config(ocn_impexp_file,numExpFields,numImpFields,expFieldName,impFieldName, &
    expStandName,impStandName,expFieldUnit,impFieldUnit, &
    expFieldEnable,impFieldEnable,irc)


    allocate(impField(numImpFields))
    allocate(expField(numExpFields))

    if(lPet.eq.0) print *,"hycom,expFieldName=",(expFieldName(i),i=1,numExpFields)
    if(lPet.eq.0) print *,"hycom,expFieldEnable=",(expFieldEnable(i),i=1,numExpFields)
    if(lPet.eq.0) print *,"hycom,impFieldName=",(impFieldName(i),i=1,numImpFields)
    if(lPet.eq.0) print *,"hycom,impFieldEnable=",(impFieldEnable(i),i=1,numImpFields)


    do i=1,numImpFields

      if(impFieldEnable(i)) then
        if(lPet.eq.0) print *,"hycom,import field advertised, name=",impFieldName(i),impStandName(i)

        call NUOPC_Advertise(importState, name=impFieldName(i), &
        StandardName=impStandName(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN
      endif
    enddo


    do i=1,numExpFields
      if(expFieldEnable(i)) then
        if(lPet.eq.0) print *,"hycom,export field advertised, name=",expFieldName(i),expStandName(i)

        call NUOPC_Advertise(exportState, name=expFieldName(i), &
        StandardName=expStandName(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT)) RETURN

      endif
    enddo
#endif

   CALL ESMF_ConfigGetAttribute(config, show_minmax, &
      label=TRIM("espc_show_impexp_minmax="), default=.true., rc=rc)

!#ifdef ESPC_SHOW_IMPEXP_MINMAX
!  show_minmax=.true.
!#else
!  show_minmax=.false.
!#endif

# ifdef ESPC_TIMER
  call ESMF_VMBarrier(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_VMWtime(timer_end, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  espc_timer(1)=timer_end-timer_beg
!  if(lPet.eq.0) print *,"hycom, InitializeP1, timer=",espc_timer(1)
  call print_timer_stat('hycom, Init1:',timer_end-timer_beg,lPet,nPets,vm,rc)

# endif

  if(lPet.eq.0) print *,"hycom, InitializeP1 end called..."

END SUBROUTINE InitializeP1

!======================================================================
#undef METHOD
#define METHOD "InitializeP2"
  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Grid)         :: gridIn
    type(ESMF_Grid)         :: gridOut

# ifdef ESPC_COUPLE
    INTEGER :: tlb(2), tub(2)
    INTEGER(ESMF_KIND_I4), POINTER :: iptr(:,:)
    REAL(ESMF_KIND_R8), POINTER :: fptr(:,:)
# endif

    integer :: i,j
    integer decomp(2)

    INTEGER :: end_hour,end_min,end_sec
    INTEGER :: start_hour,start_min,start_sec
    real*8 h_start_dtg,h_end_dtg
    integer mpiCommunicator

    type(ESMF_DistGrid) :: ocnDistGrid
    type(ESMF_DistGridConnection), allocatable, dimension(:) :: connectionList

# ifdef ESPC_COUPLE

    type(ESMF_Field)        :: lon_field, lat_field,mask_field
    real(ESMF_KIND_R8), dimension(:,:), pointer :: lon_data,lat_data,mask_data
    type(ESMF_ArraySpec) :: arraySpec2Dr
    real(kind=ESMF_KIND_R8), allocatable, dimension(:,:) :: tmp_e


    integer status

# endif

    rc = ESMF_SUCCESS



    if(lPet.eq.0) print *,"hycom, InitializeP2 called"

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_beg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif



!   HYCOM Start DTG
    CALL ESMF_ConfigGetAttribute(config, ocean_start_dtg, &
       label=TRIM("ocean_start_dtg="), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get start_dtg failed", &
     CONTEXT, rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, start_hour, &
       label=TRIM("start_hour="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get start_hour failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, start_min, &
       label=TRIM("start_min="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get start_min failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, start_sec, &
       label=TRIM("start_sec="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get start_sec failed", &
      CONTEXT,  rcToReturn=rc)) RETURN


! End time
    CALL ESMF_ConfigGetAttribute(config, end_hour, &
       label=TRIM("end_hour="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get end_hour failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, end_min, &
       label=TRIM("end_min="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get end_min failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    CALL ESMF_ConfigGetAttribute(config, end_sec, &
       label=TRIM("end_sec="), default=0, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="config get end_sec failed", &
      CONTEXT,  rcToReturn=rc)) RETURN

    if (ocean_start_dtg.lt.0.0) then
!     start from rest
      ocean_start_dtg=-ocean_start_dtg
      h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg
      h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
      h_start_dtg=-h_start_dtg
    else
!     normal start
      h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg
      h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
    endif

!    h_end_dtg=(end_hour+end_min/60.+end_sec/3600.)/24.+ocean_start_dtg
!    h_start_dtg=(start_hour+start_min/60.+start_sec/3600.)/24.+ocean_start_dtg

    if(lPet.eq.0) then
      print *,"HYCOM_OceanComp, HYCOM_Init called..."
      print *,"end_hour,end_min,end_sec=",end_hour,end_min,end_sec
      print *,"h_start_dtg,h_end_dtg=",h_start_dtg,h_end_dtg
    endif

    CALL ESMF_VMGet(vm, mpiCommunicator=mpiCommunicator,rc=rc)


!   Call into ocean init
    CALL HYCOM_Init(mpiCommunicator,h_start_dtg,h_end_dtg)


# ifdef ESPC_COUPLE
!!Alex    call ocn_import_init()
    do i=1,numImpFields
      if(impFieldEnable(i)) then
        call set_hycom_import_flag(i,impFieldName(i))
      endif
    enddo


# endif

    call hycom_couple_init(nPets,rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="hycom_couple_init failed", CONTEXT)) RETURN

    if(lPet.eq.0) print *,"hycom, InitializeP2 called 2,idim_size,jdim_size", &
       idim_size,jdim_size

#if defined(ARCTIC)
    itdmx=idim_size; jtdmx=jdim_size-1
#else
    itdmx=idim_size; jtdmx=jdim_size
#endif

    if(lPet.eq.0) print *,"hycom, itdmx,jtdmx..=",itdmx,jtdmx


# ifdef ESPC_COUPLE
    if(lPet.eq.0) then
      allocate(tmp_e(itdmx,jtdmx))
    else
      allocate(tmp_e(1,1))
    endif

#endif


# ifdef ESPC_COUPLE
    allocate(connectionList(1)) ! one connection
    call ESMF_DistGridConnectionSet(connection=connectionList(1), &
    tileIndexA=1, tileIndexB=1, positionVector=(/itdmx, 0/), rc=rc)

    ocnDistGrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/itdmx,jtdmx/), &
      indexflag=ESMF_INDEX_GLOBAL, &
      deBlockList=deBList,connectionList=connectionList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", CONTEXT)) RETURN 


    gridIn = ESMF_GridCreate(distGrid=ocnDistGrid, &
            indexflag=ESMF_INDEX_GLOBAL, &
            coordSys=ESMF_COORDSYS_SPH_DEG, &
            name="OCEAN:grid", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", CONTEXT)) RETURN 

    if(allocated(connectionList))deallocate(connectionList)

#else

   do i=int(sqrt(real(nPets))),2,-1
      if(mod(nPets,i).eq.0) exit
    enddo
    decomp=(/nPets/i,i/)
    if(lPet.eq.0) print *,"nPets,decomp=",nPets,decomp

    gridIn = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), maxIndex=(/idim_size,jdim_size/), &
            indexflag=ESMF_INDEX_GLOBAL, &
            regDecomp=decomp, name="OCEAN:grid", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid create failed", CONTEXT)) RETURN

#endif


# ifdef ESPC_COUPLE


! Add ESMF grid coordinate arrays
  CALL ESMF_GridAddCoord(gridIn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid add coord failed", CONTEXT)) RETURN

! Add ESMF mask array
  CALL ESMF_GridAddItem(gridIn, itemflag=ESMF_GRIDITEM_MASK, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid add mask failed", CONTEXT)) RETURN



    if((ocn_esmf_imp_output.or.ocn_esmf_exp_output).and.lPet.eq.0) then

      call impexp_cdf_put_latlonmsk('hycom',itdmx,jtdmx,lat_e, &
        lon_e,mask_e,status,rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_latlonmsk failed", CONTEXT)) RETURN

      call impexp_cdf_put_latlonmsk('hycom-orig',idim_size,jdim_size,lat_e, &
        lon_e,mask_e,status,rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_latlonmsk failed", CONTEXT)) RETURN

      call impexp_cdf_put_latlonmsk_corner('hycom',idim_size,jdim_size,lat_e,lon_e,int(mask_e),lat_q, &
        lon_q,status,rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_latlonmsk_corner failed", CONTEXT)) RETURN

    endif

  CALL ESMF_ArraySpecSet(arraySpec2Dr, rank=2, typekind=ESMF_TYPEKIND_R8, rc=rc)

  lon_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
                  indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                  name=TRIM("lon_field"), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  if(lPet.eq.0) then
    do i=1,itdmx
    do j=1,jtdmx
      tmp_e(i,j)=lon_e(i,j)
    enddo
    enddo
  endif

  call ESMF_FieldScatter(lon_field,tmp_e,rootPet=0,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  lat_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
                  indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                  name=TRIM("lat_field"), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  if(lPet.eq.0) then
    do i=1,itdmx
    do j=1,jtdmx
      tmp_e(i,j)=lat_e(i,j)
    enddo
    enddo
  endif

  call ESMF_FieldScatter(lat_field,tmp_e,rootPet=0,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  mask_field = ESMF_FieldCreate(gridIn, arraySpec2Dr, &
                  indexFlag=ESMF_INDEX_GLOBAL, staggerLoc=ESMF_STAGGERLOC_CENTER, &
                  name=TRIM("mask_field"), rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  if(lPet.eq.0) then
    do i=1,itdmx
    do j=1,jtdmx
      tmp_e(i,j)=mask_e(i,j)
    enddo
    enddo
  endif

  call ESMF_FieldScatter(mask_field,tmp_e,rootPet=0,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_FieldGet(lon_field,localDE,lon_data,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_FieldGet(lat_field,localDE,lat_data,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_FieldGet(mask_field,localDE,mask_data,rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN


! Copy in coordinate data
  CALL ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=1, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, &
       totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid get coord 1 failed", CONTEXT)) RETURN
!  fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_e(tlb(1):tub(1),tlb(2):tub(2))
  fptr(tlb(1):tub(1),tlb(2):tub(2)) = lon_data(tlb(1):tub(1),tlb(2):tub(2))


  CALL ESMF_GridGetCoord(gridIn, localDE=localDE, coordDim=2, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, &
       totalLBound=tlb, totalUBound=tub, farrayPtr=fptr, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid get coord 2 failed", CONTEXT)) RETURN
!  fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_e(tlb(1):tub(1),tlb(2):tub(2))
  fptr(tlb(1):tub(1),tlb(2):tub(2)) = lat_data(tlb(1):tub(1),tlb(2):tub(2))

! Copy in land/sea mask (integer)
  CALL ESMF_GridGetItem(gridIn, localDE=localDE, itemflag=ESMF_GRIDITEM_MASK, &
       staggerLoc=ESMF_STAGGERLOC_CENTER, &
       totalLBound=tlb, totalUBound=tub, farrayPtr=iptr, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="grid get mask failed", CONTEXT)) RETURN
!  iptr(tlb(1):tub(1),tlb(2):tub(2)) = NINT(mask_e(tlb(1):tub(1),tlb(2):tub(2)))
  iptr(tlb(1):tub(1),tlb(2):tub(2)) = NINT(mask_data(tlb(1):tub(1),tlb(2):tub(2)))

  call ESMF_FieldDestroy(lon_field,rc=rc)
  call ESMF_FieldDestroy(lat_field,rc=rc)
  call ESMF_FieldDestroy(mask_field,rc=rc)

# endif

    gridOut = gridIn ! for now out same as in

# ifdef ESPC_COUPLE

    do i=1,numImpFields

      if(impFieldEnable(i)) then
        if(lPet.eq.0) print *,"hycom, import field created, name=",impFieldName(i)

        impField(i) = ESMF_FieldCreate(name=impFieldName(i), grid=gridIn, &
        typekind=ESMF_TYPEKIND_RX, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

        call NUOPC_Realize(importState, field=impField(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

      endif
    enddo


    do i=1,numExpFields

      if(expFieldEnable(i)) then
        if(lPet.eq.0) print *,"hycom, export field created, name=",expFieldName(i)

      ! exportable field:
        expField(i) = ESMF_FieldCreate(name=expFieldName(i), grid=gridOut, &
        typekind=ESMF_TYPEKIND_RX, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

        call NUOPC_Realize(exportState, field=expField(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN


        if(lPet.eq.0) print *,"hycom, export field done creating, name=",expFieldName(i)

      endif
    enddo

    do i=1,numExpFields
      if(expFieldEnable(i)) then
        call do_export(i,expField(i),rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
      endif
    enddo


    endtime=0
    if(ocn_esmf_exp_output.and.(mod(endtime,float(cdf_impexp_freq)).eq.0)) then

      call impexp_cdf_put_flds('hycom', base_dtg, endtime, &
        itdmx,jtdmx, &
        numExpFields,expFieldEnable,expFieldName,expStandName,expFieldUnit,expField, &
        status,lPet,rc,'exp')
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", CONTEXT)) RETURN

    endif

    deallocate(tmp_e)

#endif

    call hycom_couple_final()

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_end, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    espc_timer(1)=espc_timer(1)+timer_end-timer_beg
!    if(lPet.eq.0) print *,"hycom, InitializeP2,timer=",espc_timer(1),timer_end-timer_beg
    call print_timer_stat('hycom, Init2:',timer_end-timer_beg,lPet,nPets,vm,rc)

# endif

    if(lPet.eq.0) print *,"hycom, InitializeP2 end called..."

END SUBROUTINE InitializeP2

!======================================================================
#undef METHOD
#define METHOD "ModelAdvance"

SUBROUTINE ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model 
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState

    integer i,status

!!!!
  TYPE(ESMF_Time) :: extCurrTime
  TYPE(ESMF_Time) :: extRefTime
  Type(ESMF_TimeInterval) :: extTimeStep
  character(ESMF_MAXSTR) :: currtimeString,reftimeString
  TYPE(ESMF_TimeInterval) :: extTimeSinceStart
  real(kind=ESMF_KIND_R8) :: extSecSinceStarti
  real(kind=ESMF_KIND_R8) :: extSecTimeStep
  real  begtime
  real*8 endtime8
  real endtimex

# ifdef ESPC_TIMER
  real(kind=ESMF_KIND_R8) :: timer_tmp_beg,timer_tmp_end
# endif

    rc = ESMF_SUCCESS

  if(lPet.eq.0) print *,"hycom, ModelAdvance called"
!  return


# ifdef ESPC_TIMER
  call ESMF_VMBarrier(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_VMWtime(timer_beg, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

!!!!!!
! Compute the step count from the external clock
  CALL ESMF_ClockGet(clock, currTime=extCurrTime, refTime=extRefTime, &
      timeStep=extTimeStep,rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="Get extCurrTime failed", CONTEXT)) RETURN


  call ESMF_TimeGet(extCurrTime,timeString=currtimeString)
  call ESMF_TimeGet(extRefTime,timeString=reftimeString)

!  if(lPet.eq.0) print *,"hycom,extCurrTime, extRefTime=",currtimeString,reftimeString

! Elapsed time in seconds for internal clock
  extTimeSinceStart = extCurrTime - extRefTime
  CALL ESMF_TimeIntervalGet(extTimeSinceStart, s_r8=extSecSinceStarti, rc=rc)
  CALL ESMF_TimeIntervalGet(extTimeStep, s_r8=extSecTimeStep, rc=rc)
  IF (ESMF_LogFoundError(rcToCheck=rc, msg="Get time interval failed", CONTEXT)) RETURN

  if(lPet.eq.0) print *,"hycom,extSecTimeStep=",extSecTimeStep

  begtime=extSecSinceStarti/3600.
  endtime=(extSecSinceStarti+extSecTimeStep)/3600.

  if(lPet.eq.0) print *,"hycom,begtime,endtime=",begtime,endtime

  endtime=endtime+ocean_start_dtg*24


! Run atmos forward

  if(lPet.eq.0) print *,"HYCOM_OceanCom, ModelAdvance, Run ocean forward...,endtime=",endtime/24


# ifdef ESPC_COUPLE

    endtimex=endtime-ocean_start_dtg*24


    if(ocn_esmf_imp_output.and.( (mod(endtimex,float(cdf_impexp_freq)).eq.0 .or. &
      endtimex.eq.0.5)  )) then

      call impexp_cdf_put_flds('hycom', base_dtg, endtimex-cpl_time_step, &
        itdmx,jtdmx, &
        numImpFields,impFieldEnable,impFieldName,impStandName,impFieldUnit,impField,status, &
        lPet,rc,'imp')
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", CONTEXT)) RETURN

    endif

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif


    do i=1,numImpFields
      if(impFieldEnable(i)) then
        call do_import(i,impField(i),.false.,rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
      endif
    enddo

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    espc_timer(4)=espc_timer(4)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Import):',timer_tmp_end-timer_tmp_beg,lPet,nPets,vm,rc)

# endif

!    call redef_radflx()

!ajw
!move to mod_hycom.F
#  ifdef ESPC_COUPLE
    call ocn_import_forcing()
    if(hycom_arche_output .and. begtime.eq.0) &
      call archiv_exchange  !arche file of fields exchanged with ice component
#endif

# endif

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

    ! Because of the way that the internal Clock was set in InitializeP2()
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    if(lPet.eq.0) then
      call ESMF_ClockPrint(clock, options="currTime",&
        preString="------>Advancing OCN from: ", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

      call ESMF_ClockPrint(clock, options="stopTime",&
        preString="--------------------------------> to: ", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
    endif


    endtime8=endtime

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif

  do !until end of run
    CALL HYCOM_Run(endtime8/24)
    if     (end_of_run .or. end_of_run_cpl ) then
      exit
    endif

  enddo

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    espc_timer(5)=espc_timer(5)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Core):',timer_tmp_end-timer_tmp_beg,lPet,nPets,vm,rc)

# endif

# ifdef ESPC_COUPLE

# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_beg, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif

    do i=1,numExpFields
      if(expFieldEnable(i)) then
        call do_export(i,expField(i),rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

      endif
    enddo



# ifdef ESPC_TIMER
    call ESMF_VMBarrier(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_VMWtime(timer_tmp_end, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    espc_timer(6)=espc_timer(6)+timer_tmp_end-timer_tmp_beg
    call print_timer_stat('hycom, Run(Export):',timer_tmp_end-timer_tmp_beg,lPet,nPets,vm,rc)

# endif

!move to mod_hycom.F
# ifdef ESPC_COUPLE
    if(hycom_arche_output) call archiv_exchange  !arche file of fields exchanged with ice component
# endif

    if(ocn_esmf_exp_output.and.( (mod(endtimex,float(cdf_impexp_freq)).eq.0 .or. &
      endtimex.eq.0.5)  )) then

      call impexp_cdf_put_flds('hycom', base_dtg,endtimex, &
        itdmx,jtdmx, &
        numExpFields,expFieldEnable,expFieldName,expStandName,expFieldUnit,expField, &
        status,lPet,rc,'exp')
      IF (ESMF_LogFoundError(rcToCheck=rc, msg="impexp_cdf_put_flds failed", CONTEXT)) RETURN

    endif

# endif

# ifdef ESPC_TIMER
  call ESMF_VMBarrier(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_VMWtime(timer_end, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  espc_timer(2)=espc_timer(2)+timer_end-timer_beg
!  if(lPet.eq.0) print *,"hycom, ModelAdvance,timer=",espc_timer(2), &
!     timer_end-timer_beg
  call print_timer_stat('hycom, Run:',timer_end-timer_beg,lPet,nPets,vm,rc)

# endif

  if(lPet.eq.0) print *,"hycom, ModelAdvance end..."
  if(lPet.eq.0) print *,"hycom, ModelAdvance end...",begtime,endtime
END SUBROUTINE ModelAdvance



!======================================================================
#undef METHOD
#define METHOD "OCEAN_Final"
SUBROUTINE OCEAN_Final(model, rc)

  TYPE(ESMF_GridComp) :: model
  INTEGER,INTENT(OUT) :: rc

! Local variables
  INTEGER :: lrc, i
# ifdef ESPC_TIMER
  integer j,ij
  real(kind=ESMF_KIND_R8), allocatable:: espc_all_timer(:)
  real, allocatable:: all_timer(:,:)
  real timer_min(6),timer_max(6),timer_mean(6),timer_stdev(6)
# endif
! Assume failure
  rc = ESMF_FAILURE

! Report
  CALL ESMF_LogWrite("HYCOM finalize routine called", ESMF_LOGMSG_INFO)

# ifdef ESPC_TIMER
  call ESMF_VMBarrier(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_VMWtime(timer_beg, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN
# endif



! Finalize ocean 
   if(lPet.eq.0) print *,"HYCOM Final called.."
   CALL HYCOM_Final

# ifdef ESPC_COUPLE

    if (associated(expFieldName)) deallocate(expFieldName)
    if (associated(impFieldName)) deallocate(impFieldName)
    if (associated(expStandName)) deallocate(expStandName)
    if (associated(impStandName)) deallocate(impStandName)
    if (associated(expFieldUnit)) deallocate(expFieldUnit)
    if (associated(impFieldUnit)) deallocate(impFieldUnit)
    !if (associated(expFieldAddOffset)) deallocate(expFieldAddOffset)
    !if (associated(impFieldAddOffset)) deallocate(impFieldAddOffset)
    !if (associated(expFieldScaleFac)) deallocate(expFieldScaleFac)
    !if (associated(impFieldScaleFac)) deallocate(impFieldScaleFac)
    if (associated(expFieldEnable)) deallocate(expFieldEnable)
    if (associated(impFieldEnable)) deallocate(impFieldEnable)


    deallocate(impField)
    deallocate(expField)

# endif

# ifdef ESPC_TIMER
  call ESMF_VMBarrier(vm, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  call ESMF_VMWtime(timer_end, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

  espc_timer(3)=timer_end-timer_beg

!  if(lPet.eq.0) print *,"hycom, Final,timer=",espc_timer(3)
  call print_timer_stat('hycom, Final:',timer_end-timer_beg,lPet,nPets,vm,rc)

  call print_timer_stat('       HYCOM_Init Phase:',espc_timer(1),lPet,nPets,vm,rc)
  call print_timer_stat('        HYCOM_Run Phase:',espc_timer(2),lPet,nPets,vm,rc)
  call print_timer_stat('      HYCOM_Final Phase:',espc_timer(3),lPet,nPets,vm,rc)
  call print_timer_stat('HYCOM_Run Phase(Import):',espc_timer(4),lPet,nPets,vm,rc)
  call print_timer_stat('  HYCOM_Run Phase(Core):',espc_timer(5),lPet,nPets,vm,rc)
  call print_timer_stat('HYCOM_Run Phase(Export):',espc_timer(6),lPet,nPets,vm,rc)


# endif

  if(lPet.eq.0) print *,"hycom, OCEAN_Final end..."

! Return success
  rc = ESMF_SUCCESS

!END SUBROUTINE HYCOM_OCEAN_Final
END SUBROUTINE OCEAN_Final


# ifdef ESPC_COUPLE

!=======================================================================================
#undef METHOD
#define METHOD "do_export"
SUBROUTINE do_export(k,field,rc)

    integer k,i,j
    type(ESMF_Field)        :: field

    real*8, allocatable, dimension(:,:) :: expData
    real(ESMF_KIND_RX), dimension(:,:), pointer :: field_data
    integer tlb(2),tub(2)

    integer :: rc
    character(len=30) fieldName
    integer status

    fieldName=expFieldName(k)


    call ESMF_FieldGet(field,localDE,field_data,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_FieldGetBounds(field, localDE=localDE,  &
      totalLBound=tlb, totalUBound=tub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) RETURN


!    write(*,990) lPet,tlb(1),tub(1),1+i0,ii+i0,tlb(2),tub(2),1+j0,jj+j0
!990 format('lPet,tlb...=',I4, 4I6,6x,4I6)
!    write(*,991) lPet, i0,ii,j0,jj
!991 format('lPet, i0...=',I4,4I6)

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

    field_data(:,:)=0.

    allocate(expData(tlb(1):tub(1),tlb(2):tub(2) ))
    call export_from_hycom_deb(tlb,tub,expData,fieldName,show_minmax)

    do j   = tlb(2),tub(2)
    do i   = tlb(1),tub(1)
        field_data(i,j)=expData(i,j)
    enddo
    enddo

    if(allocated(expData)) deallocate(expData)


!  Return success
   rc = ESMF_SUCCESS

END SUBROUTINE do_export

!=======================================================================================
#undef METHOD
#define METHOD "do_import"
SUBROUTINE do_import(k,field,data_init_flag,rc)

    integer k,i,j
    type(ESMF_Field)        :: field

    real*8, allocatable, dimension(:,:) :: impData
    real(ESMF_KIND_RX), dimension(:,:), pointer :: field_data
    integer tlb(2),tub(2)

    integer :: rc,status
    character(len=30) fieldName
    logical data_init_flag

    fieldName=impFieldName(k)


    call ESMF_FieldGet(field,localDE,field_data,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU,CONTEXT)) RETURN

    call ESMF_FieldGetBounds(field, localDE=localDE,  &
      totalLBound=tlb, totalUBound=tub, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, &
      rcToReturn=rc)) RETURN

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

    allocate(impData(tlb(1):tub(1),tlb(2):tub(2) ))

    impData(:,:)=0.

    do j   = tlb(2),tub(2)
    do i   = tlb(1),tub(1)
        impData(i,j)=field_data(i,j)
    enddo
    enddo


    call import_to_hycom_deb(tlb,tub,impData,fieldName,show_minmax,data_init_flag)

    if(allocated(impData)) deallocate(impData)


!  Return success
   rc = ESMF_SUCCESS

END SUBROUTINE do_import

# endif


!=========================================================================================
!END MODULE HYCOM_OCEAN_Mod
#if defined ESPC_OCN
END MODULE OCEAN_Mod
#else
END MODULE HYCOM_Mod
#endif


