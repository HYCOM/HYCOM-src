module mod_hycom_nuopc_glue

  !-----------------------------------------------------------------------------
  ! NUOPC glue code for HYCOM
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  use mod_xc
  use mod_cb_arrays_nuopc_glue

  implicit none

  integer :: jja, jtdma 
  
  private
  
  type hycom_nuopc_glue_type
    ! native grid:
    type(ESMF_Grid)         :: grid
    ! native mesh:
    type(ESMF_Mesh)         :: mesh           ! from Grid center stagger 
    ! import fields:
    type(ESMF_FieldBundle)  :: importFields
    ! export fields:
    type(ESMF_FieldBundle)  :: exportFields
#ifdef WORKAROUND_HOLES
    ! --------------------------------------------------------------------------
    ! Below members are only temporarily needed to work-around a Regrid 
    ! limitation wrt holes in the src/dst DistGrid. The work-around ensures
    ! that the Fields listed above this comment are defined on DistGrids
    ! without index space holes. The Fields on the native HYCOM decomposition
    ! (which has holes in index space) are held in the _shadow Fields below.
    ! Extra _shadow RouteHandles are necessary to deal with the extra level of
    ! indirection.
    ! --------------------------------------------------------------------------
    ! shadow grid:
    type(ESMF_Grid)         :: shadow_grid
    ! import shadow fields:
    type(ESMF_FieldBundle)  :: shadow_importFields
    ! export shadow fields:
    type(ESMF_FieldBundle)  :: shadow_exportFields
    ! other members
    logical                 :: rh_import2shadow_ready, rh_shadow2export_ready
    type(ESMF_RouteHandle)  :: rh_import2shadow
    type(ESMF_RouteHandle)  :: rh_shadow2export
#endif
  endtype
  
  public hycom_nuopc_glue_type
  public HYCOM_TileInfo, HYCOM_GlueInitialize
  public HYCOM_GlueFieldRealize, HYCOM_GlueFieldsRealize
  public HYCOM_GlueFieldsDataImport, HYCOM_GlueFieldsDataExport
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine HYCOM_TileInfo(rc)
    integer, intent(out), optional :: rc
    
    character(len=1600)  :: msg
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    call ESMF_LogWrite("<HYCOM Tile Info>", ESMF_LOGMSG_INFO)
    
    write (msg, *) "ipr=", ipr, "   jpr=", jpr, "   ijpr=", ijpr
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
    write (msg, *) "mproc=", mproc, "   nproc=", nproc, "   mnproc=", mnproc
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
    write (msg, *) "i0=", i0, "   ii=", ii, "   j0=", j0, "   jj=", jja
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)
!!Alex    write (msg, *) "margin=", margin, "   nreg=", nreg, "   vland=", vland
    write (msg, *) "nbdy=", nbdy, "   nreg=", nreg, "   vland=", vland
    call ESMF_LogWrite(msg, ESMF_LOGMSG_INFO)

    call ESMF_LogWrite("</HYCOM Tile Info>", ESMF_LOGMSG_INFO)
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine HYCOM_GlueInitialize(glue, rc)
    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    integer, intent(out), optional              :: rc
    
    integer, allocatable, target  :: deBlockList(:,:,:)
    integer, allocatable, target  :: deBlockList_buf(:,:,:)
    integer, pointer      :: sendData(:), recvData(:)
    type(ESMF_VM)         :: vm
    integer               :: localPet, petCount
    type(ESMF_DistGridConnection), allocatable :: connectionList(:)
    type(ESMF_DistGrid)   :: dg
    type(ESMF_Array)      :: array_plon, array_plat, array_msk, array_area
    type(ESMF_Array)      :: array_qlon, array_qlat
    
    real(kind=ESMF_KIND_R8) :: dump_lat(1500,1100)
    real(kind=ESMF_KIND_R8) :: dump_lon(1500,1100)
    integer                 :: dump_msk(1500,1100)
    integer                 :: i, j, dumpUnit
    
    real(kind=ESMF_KIND_R8) :: tst_lat(1500,1100)
    real(kind=ESMF_KIND_R8) :: tst_lon(1500,1100)
    real(kind=ESMF_KIND_R8) :: tst_msk(1500,1100)
    real(kind=ESMF_KIND_R8), allocatable :: tst_lmsk(:,:)
    
    real(kind=ESMF_KIND_R8), pointer  :: farrayPtr(:,:)
    integer, pointer                  :: farrayPtrI(:,:)

    character(len=80)       :: dumpFile
    
#ifdef WORKAROUND_HOLES
    type(ESMF_DistGrid) :: dg_worka
    type(ESMF_Array)    :: array_plon_worka, array_plat_worka, array_msk_worka
    type(ESMF_RouteHandle)  :: rh
#endif
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (localPet /= mnproc-1) then
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="localPet and mnproc must correspond!", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return
    endif

    ! change the domain if HYCOM Arctic patch   
    if (ltripolar) then
       jtdma = jtdm-1
       jja   = min(jj,jtdma-j0)
    else
       jtdma = jtdm
       jja   = jj
    endif
   
    ! prepare the deBlockList needed for DistGrid creation
    
    ! first step: set the local piece of information
    
    allocate(deBlockList(2, 2, ijpr)) ! dimCount, 2, deCount
    allocate(deBlockList_buf(2, 2, ijpr)) ! dimCount, 2, deCount
    
    deBlockList(1, 1, mnproc) = i0+1   ! minIndex 1st dim
    deBlockList(2, 1, mnproc) = j0+1   ! minIndex 2nd dim
    
    deBlockList(1, 2, mnproc) = i0+ii  ! maxIndex 1st dim
    deBlockList(2, 2, mnproc) = j0+jja ! maxIndex 2nd dim

    deBlockList_buf = deBlockList

    ! second step: all-to-all the information so every PET has full deBlockList

    sendData => deBlockList(:,1,mnproc)
    recvData => deBlockList_buf(:,1,1)
    
    call ESMF_VMAllGather(vm, sendData, recvData, 4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    deBlockList = deBlockList_buf
    deallocate(deBlockList_buf)

#if 0
    write(dumpFile,"('deBlockList',I3.3,'.dat')") localPet
    call ESMF_UtilIOUnitGet(dumpUnit)
    open(unit=dumpUnit, file=dumpFile, status="replace", action="write")
    do i=1,petCount
      write (dumpUnit, "(5I5)") i, deBlockList(1,1,i), deBlockList(2,1,i), &
        deBlockList(1,2,i), deBlockList(2,2,i)
    enddo
    close(unit=dumpUnit)
#endif

    ! prepare a single connection for periodicity along i-axis (longitude)
    allocate(connectionList(1))
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/itdm, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! ready to create the HYCOM DistGrid from deBlockList with periodic connect.
    dg = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/itdm,jtdma/), &
      deBlockList=deBlockList, connectionList=connectionList, &
      indexflag=ESMF_INDEX_GLOBAL, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    deallocate(connectionList)
    deallocate(deBlockList)
    
    ! dress up "plon" array, considering HYCOM memory layout with halo + padding
    array_plon = ESMF_ArrayCreate(dg, farray=plon, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! dress up "plat" array, considering HYCOM memory layout with halo + padding
    array_plat = ESMF_ArrayCreate(dg, farray=plat, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
     
    ! dress up "qlon" array, considering HYCOM memory layout with halo + padding
    array_qlon = ESMF_ArrayCreate(dg, farray=qlon, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! dress up "qlat" array, considering HYCOM memory layout with halo + padding
    array_qlat = ESMF_ArrayCreate(dg, farray=qlat, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
 
    ! dress up "mask" array, considering HYCOM memory layout with halo + padding
    array_msk = ESMF_ArrayCreate(dg, farray=ip, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! dress up "area" array, considering HYCOM memory layout with halo + padding
    array_area = ESMF_ArrayCreate(dg, farray=scp2, &
      indexflag=ESMF_INDEX_DELOCAL, &
      computationalLWidth=(/nbdy,nbdy/), computationalUWidth=(/nbdy,nbdy/), &
      totalLWidth=(/nbdy,nbdy/), & ! lower corner pinned to memory alloc, float upper
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! ready to create the HYCOM Grid from DistGrid and coordinate Arrays
    glue%grid = ESMF_GridCreate(dg, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger longitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=1, array=array_plon, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger latitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=2, array=array_plat, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the corner stagger latitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      coordDim=1, array=array_qlon, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger latitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CORNER, &
      coordDim=2, array=array_qlat, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger mask Array
    call ESMF_GridSetItem(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, array=array_msk, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger area Array
    call ESMF_GridSetItem(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_AREA, array=array_area, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! attach the center staggered mesh too for later reference
    glue%mesh = ESMF_GridToMesh(glue%grid, ESMF_STAGGERLOC_CENTER, 1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    return ! bail out
     
    ! create the import and export FieldBundles
    glue%importFields = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    glue%exportFields = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
#ifdef WORKAROUND_HOLES
    ! Push the native field_sst into the "shadow" to allow the work-around 
    ! Field (without holes in the DistGrid) play the part of field_sst.
    glue%shadow_grid = glue%grid
    ! Create work-around DistGrid.
    allocate(connectionList(1))
    call ESMF_DistGridConnectionSet(connectionList(1), tileIndexA=1, &
      tileIndexB=1, positionVector=(/itdm, 0/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    dg_worka = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/itdm,jtdma/), &
      connectionList=connectionList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    deallocate(connectionList)
    ! Create work-around coordinate and mask Arrays.
    array_plon_worka = ESMF_ArrayCreate(dg_worka, ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    array_plat_worka = ESMF_ArrayCreate(dg_worka, ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    array_msk_worka = ESMF_ArrayCreate(dg_worka, ESMF_TYPEKIND_I4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Redistribute coordinates from HYCOM-native distribution to work-around
    call ESMF_ArrayRedistStore(array_plon, array_plon_worka, routehandle=rh, &
      rc=rc) ! this rh will work for plon & plat
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(array_plon_worka, farrayPtr=farrayPtr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    farrayPtr = 0.D0  ! zero out before Redist
    call ESMF_ArrayRedist(array_plon, array_plon_worka, routehandle=rh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(array_plat_worka, farrayPtr=farrayPtr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    farrayPtr = 0.D0  ! zero out before Redist
    call ESMF_ArrayRedist(array_plat, array_plat_worka, routehandle=rh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ArrayRedistStore(array_msk, array_msk_worka, routehandle=rh, &
      rc=rc) ! mask needs its own rh because it is of different typekind
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayGet(array_msk_worka, farrayPtr=farrayPtrI, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    farrayPtrI = 0  ! zero out before Redist
    call ESMF_ArrayRedist(array_msk, array_msk_worka, routehandle=rh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! ready to create the HYCOM Grid from DistGrid and coordinate Arrays
    glue%grid = ESMF_GridCreate(dg_worka, coordSys=ESMF_COORDSYS_SPH_DEG, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger longitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=1, array=array_plon_worka, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger latitude coordinate Array
    call ESMF_GridSetCoord(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      coordDim=2, array=array_plat_worka, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set the center stagger mask Array
    call ESMF_GridSetItem(glue%grid, staggerLoc=ESMF_STAGGERLOC_CENTER, &
      itemflag=ESMF_GRIDITEM_MASK, array=array_msk_worka, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! create the shadow import and export FieldBundles
    glue%shadow_importFields = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    glue%shadow_exportFields = ESMF_FieldBundleCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! reset the RH guard flags
    glue%rh_import2shadow_ready = .false.
    glue%rh_shadow2export_ready = .false.
#endif

    ! initialize coupling flags
    cpl_taux      =.false.
    cpl_tauy      =.false.
    cpl_wndspd    =.false.
    cpl_ustara    =.false.
    cpl_airtmp    =.false.
    cpl_vapmix    =.false.
    cpl_swflx     =.false.
    cpl_lwmdnflx  =.false.
    cpl_lwmupflx  =.false.
    cpl_latflx    =.false.
    cpl_sensflx   =.false.
    cpl_precip    =.false.
    cpl_surtmp    =.false.
    cpl_seatmp    =.false.
    cpl_implicit  =.false.
    cpl_orivers   =.false.
    cpl_irivers   =.false.
   
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine HYCOM_GlueFieldsRealize(glue, state, standardNames, rc)
    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    type(ESMF_State)                            :: state
    character(len=*)                            :: standardNames(:)
    integer, intent(out), optional              :: rc
    
    ! local variables
    integer                           :: i
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    do i=1, size(standardNames)

      call ESMF_LogWrite(trim('Realize glue field: '//standardNames(i)), ESMF_LOGMSG_INFO, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

      call HYCOM_GlueFieldRealize(glue, state, &
        standardName=standardNames(i), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine HYCOM_GlueFieldRealize(glue, state, standardName, rc)
    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    type(ESMF_State)                            :: state
    character(len=*)                            :: standardName
    integer, intent(out), optional              :: rc
    
    ! local variables
    logical                           :: connected
    character(len=80)                 :: fieldName
    type(ESMF_Field)                  :: field, shadow
    type(ESMF_StateIntent_Flag)       :: stateIntent
    real(kind=ESMF_KIND_R8), pointer  :: farrayPtr(:,:)
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    ! determine shortName from field dictionary, to be used as fieldName
    fieldName=standardName

    call ESMF_StateGet(state, stateintent=stateIntent, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (stateIntent == ESMF_STATEINTENT_UNSPECIFIED) then
      ! not a valid intent
      call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
        msg="The incoming state must have import or export intent!", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif

    connected = .true.
    
    if (connected) then
      call ESMF_LogWrite(trim('HYCOM: Create connected glue field: '//standardName), ESMF_LOGMSG_INFO, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      ! create the Field object
      field = ESMF_FieldCreate(glue%grid, ESMF_TYPEKIND_R8, name=fieldName, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      farrayPtr = -999999999999999999.9999999999999 
      call ESMF_AttributeSet(field, name="StandardName", value=standardName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      ! add the Field to the correct glue FieldBundle
      if (stateIntent == ESMF_STATEINTENT_IMPORT) then
        call ESMF_LogWrite(trim('HYCOM: Add connected import glue field: '//standardName), ESMF_LOGMSG_INFO, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
        ! import
        call ESMF_FieldBundleAdd(glue%importFields, fieldList=(/field/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        call ESMF_LogWrite(trim('HYCOM: Add connected export glue field: '//standardName), ESMF_LOGMSG_INFO, rc=rc)
        if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
        ! export
        call ESMF_FieldBundleAdd(glue%exportFields, fieldList=(/field/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
      
#ifdef WORKAROUND_HOLES
      
      ! deal with shadow objects
      
      shadow = ESMF_FieldCreate(glue%shadow_grid, ESMF_TYPEKIND_R8, &
        name=fieldName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! add the Field to the correct shadow glue FieldBundle
      if (stateIntent == ESMF_STATEINTENT_IMPORT) then
        ! import
        call ESMF_FieldBundleAdd(glue%shadow_importFields, &
          fieldList=(/shadow/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! potentially precompute the Redist RH
        if (.not.glue%rh_import2shadow_ready) then
          ! zero out src side
          call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          farrayPtr = 0.D0  ! zero out before Redist
          ! precompute RH import2shadow
          ! Requires the use of the "ignoreUnmatchedIndices" argument b/c dst 
          ! DistGrid has holes in the index space and not every src DistGrid 
          ! elements has an associated dst DistGrid element.
          call ESMF_FieldRedistStore(field, shadow, &
            routehandle=glue%rh_import2shadow, ignoreUnmatchedIndices=.true., &
            rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! set flag
          glue%rh_import2shadow_ready = .true. 
        endif
      else
        ! export
        call ESMF_FieldBundleAdd(glue%shadow_exportFields, &
          fieldList=(/shadow/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
#define WHILE_REDIST_HAS_NO_ZEROREGIONFLAG
#ifdef WHILE_REDIST_HAS_NO_ZEROREGIONFLAG
        ! zero out the entire field was added to the exportState
        ! This is strictly not necessary, because the Connector that connects
        ! to this Field is supposed to take into consideration the mask, and
        ! therefore uninitialized elements (where holes are on the native grid)
        ! are not used during Regrid in the Connector. However, when looking at
        ! the diagnostic output field files there will be noise showing up over
        ! the holes, which is is disturbing, and can be eliminated by zero'ing
        ! out the Fields when they are created here.
        call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        farrayPtr = 0.D0  ! zero out so that areas with holes hold zero
#endif
        ! potentially precompute the Redist RH
        if (.not.glue%rh_shadow2export_ready) then
          ! zero out src side
          call ESMF_FieldGet(shadow, farrayPtr=farrayPtr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          farrayPtr = 0.D0  ! zero out before Redist
          ! precompute RH shadow2export
          call ESMF_FieldRedistStore(shadow, field, &
            routehandle=glue%rh_shadow2export, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! set flag
          glue%rh_shadow2export_ready = .true. 
        endif
      endif
      
#endif

    else
      ! remove a not connected Field from State
      !call ESMF_StateRemove(state, (/fieldName/), rc=rc)
      !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      !  line=__LINE__, &
      !  file=__FILE__)) &
      !  return  ! bail out
    endif
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine HYCOM_GlueFieldsDataImport(glue, initFlag, rc)
    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    logical                                     :: initFlag
    integer, intent(out), optional              :: rc
    
    integer                           :: fieldCount, iField, stat
    type(ESMF_Field), allocatable     :: fieldList(:)
    character(len=80), allocatable    :: fieldNameList(:)
    character(len=80)                 :: fieldName, fieldStdName
    type(ESMF_Field)                  :: field, shadow
    type(ESMF_StateIntent_Flag)       :: stateIntent
    real(kind=ESMF_KIND_R8), pointer  :: farrayPtr(:,:)
    real(kind=ESMF_KIND_R8), pointer  :: impPtr(:,:), impPtr2(:,:,:)
    integer                           :: i,j
    logical                           :: twoLevel
    real, parameter                   :: sstmin = -1.8 ! Celsius
    real, parameter                   :: sstmax = 35.0 ! Celsius
    logical                           :: isConnected
    
    if (present(rc)) rc = ESMF_SUCCESS
    
    ! access the members inside of the importFields FieldBundle
    call ESMF_FieldBundleGet(glue%importFields, fieldCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(fieldList(fieldCount), fieldNameList(fieldCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocations failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldBundleGet(glue%importFields, fieldList=fieldList, &
      fieldNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! loop over all of the import Fields and pull the data in
    do iField=1, fieldCount
   
      field = fieldList(iField)
      fieldName = fieldNameList(iField)
    
      call ESMF_AttributeGet(field, name="StandardName", value=fieldStdName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      isConnected = .false.
      call ESMF_AttributeGet(field, name="HYCOM_IN_CESM_Connected", value=isConnected, defaultValue=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      if(.not. isConnected) cycle
      !call ESMF_LogWrite("HYCOM_GlueFieldsDataImport(): "// &
      !  trim(fieldStdName)//" - "//trim(fieldName), ESMF_LOGMSG_INFO)
        
#ifndef WORKAROUND_HOLES
      ! Access the HYCOM distributed field.
      call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#else
      ! look for this field in the shadow_importFields
      call ESMF_FieldBundleGet(glue%shadow_importFields, fieldName=fieldName, &
        field=shadow, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! Redist from import to shadow
      call ESMF_FieldRedist(field, shadow, &
        routehandle=glue%rh_import2shadow, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! Access the HYCOM distributed field.
      call ESMF_FieldGet(shadow, farrayPtr=farrayPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif
 
      ! identify the exact import field by standard name and copy the data
      twoLevel = .false. ! reset
      if (fieldStdName == "surface_downward_eastward_stress") then
        cpl_taux = .not.initFlag
        impPtr2 => imp_taue
        twoLevel = .true.
      elseif (fieldStdName == "surface_downward_northward_stress") then
        cpl_tauy = .not.initFlag
        impPtr2 => imp_taun
        twoLevel = .true.
      elseif (fieldStdName == "wind_speed_height10m") then
        cpl_wndspd = .not.initFlag
        impPtr2 => imp_wndspd
        twoLevel = .true.
      elseif (fieldStdName == "friction_speed") then
        cpl_ustara = .not.initFlag
        impPtr2 => imp_ustara
        twoLevel = .true.
      elseif (fieldStdName == "mean_net_sw_flx") then
        cpl_swflx = .not.initFlag
        impPtr2 => imp_swflx
        twoLevel = .true. 
      elseif (fieldStdName == "mean_down_lw_flx") then
        cpl_lwmdnflx = .not.initFlag
        impPtr2 => imp_lwdflx
        twoLevel = .true.
      elseif (fieldStdName == "mean_up_lw_flx") then
        cpl_lwmupflx = .not.initFlag
        impPtr2 => imp_lwuflx
        twoLevel = .true.
      elseif (fieldStdName == "mean_lat_flx") then
        cpl_latflx = .not.initFlag
        impPtr2 => imp_latflx
        twoLevel = .true.
      elseif (fieldStdName == "mean_sens_flx") then
        cpl_sensflx = .not.initFlag
        impPtr2 => imp_sensflx
        twoLevel = .true.
      elseif (fieldStdName == "inst_temp_height2m") then
        cpl_airtmp = .not.initFlag
        impPtr2 => imp_airtmp
        twoLevel = .true.
      elseif (fieldStdName == "mean_prec_rate") then
        cpl_precip = .not.initFlag
        impPtr2 => imp_precip
        twoLevel = .true.
      elseif (fieldStdName == "inst_spec_humid_height2m") then
        cpl_vapmix = .not.initFlag
        impPtr2 => imp_vapmix
        twoLevel = .true.
      elseif (fieldStdName == "sea_surface_temperature") then
        cpl_surtmp = .not.initFlag
        impPtr2 => imp_surtmp
        twoLevel = .true.
      elseif (fieldStdName == "sea_ice_area_fraction") then
        cpl_sic = .true.
        impPtr => sic_import
        twoLevel = .false.
      elseif (fieldStdName == "downward_x_stress_at_sea_ice_base") then
        cpl_sitx = .true.
        impPtr => sitx_import
        twoLevel = .false.
      elseif (fieldStdName == "downward_y_stress_at_sea_ice_base") then
        cpl_sity = .true.
        impPtr => sity_import
        twoLevel = .false.
      elseif (fieldStdName == "downward_sea_ice_basal_solar_heat_flux") then
        cpl_siqs = .true.
        impPtr => siqs_import
        twoLevel = .false.
      elseif (fieldStdName == "upward_sea_ice_basal_heat_flux") then
        cpl_sifh = .true.
        impPtr => sifh_import
        twoLevel = .false.
      elseif (fieldStdName == "downward_sea_ice_basal_salt_flux") then
        cpl_sifs = .true.
        impPtr => sifs_import
        twoLevel = .false.
      elseif (fieldStdName == "downward_sea_ice_basal_water_flux") then
        cpl_sifw = .true.
        impPtr => sifw_import
        twoLevel = .false.
      elseif (fieldStdName == "sea_ice_temperature") then
        cpl_sit = .true.
        impPtr => sit_import
        twoLevel = .false.
      elseif (fieldStdName == "sea_ice_thickness") then
        cpl_sih = .true.
        impPtr => sih_import
        twoLevel = .false.
      elseif (fieldStdName == "sea_ice_x_velocity") then
        cpl_siu = .true.
        impPtr => siu_import
        twoLevel = .false.
      elseif (fieldStdName == "sea_ice_y_velocity") then
        cpl_siv = .true.
        impPtr => siv_import
        twoLevel = .false.
!To Alex, connect the water flux from river here then 
! convert Kg/m^2/s -> m/s as needed below for these two fields
      elseif (fieldStdName == "water_flux_into_sea_water") then
        cpl_orivers = .true.
        impPtr2 => imp_orivers
        twoLevel = .true.
      elseif (fieldStdName == "frozen_water_flux_into_sea_water") then
        cpl_irivers = .true.
        impPtr2 => imp_irivers
        twoLevel = .true.
      endif
      
      
      ! copy the data into the right import location
      if (twoLevel) then
        if (initFlag) then
          ! initial condition set #2 to initial
          do j=1,jja
          do i=1,ii
            impPtr2(i,j,2) = farrayPtr(i,j)
          enddo
          enddo
        else
          ! shift #1 -> #2
          do j=lbound(impPtr2,2), ubound(impPtr2,2)
           do i=lbound(impPtr2,1), ubound(impPtr2,1)
            impPtr2(i,j,2) = impPtr2(i,j,1)
           enddo
          enddo
        endif
        ! fill #1
        do j=1,jja
        do i=1,ii
          impPtr2(i,j,1) = farrayPtr(i,j)
        enddo
        enddo
        
        ! special treatment for setting the internal sst
        if (fieldStdName=="sea_surface_temperature") then
          if (initFlag) then
            ! unit change: K -> C
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) - 273.15
              impPtr2(i,j,2) = impPtr2(i,j,2) - 273.15
            enddo
            enddo
          else
            ! unit change: K -> C
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) - 273.15
            enddo
            enddo
          endif
          if (sstflg/=3) then
            ! consider atmos air surface temperature
            cpl_seatmp = .true.
            do j=1,jja
            do i=1,ii
              imp_seatmp(i,j,1) = max(sstmin, min(imp_surtmp(i,j,1), sstmax))
              imp_seatmp(i,j,2) = max(sstmin, min(imp_surtmp(i,j,2), sstmax))
            enddo
            enddo
          endif
        
        ! special treatment for temperature unit change
        else if (fieldStdName=="inst_temp_height2m") then
          if (initFlag) then
            ! unit change: K -> C
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) - 273.15
              impPtr2(i,j,2) = impPtr2(i,j,2) - 273.15
            enddo
            enddo
          else
            ! unit change: K -> C
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) - 273.15
            enddo
            enddo
          endif
          
        ! special treatment for precipitation rate
        else if (fieldStdName=="mean_prec_rate") then
          if (initFlag) then
            ! unit change: mm/s (same as kg s-1 m-2) -> m/s
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.D-3
              impPtr2(i,j,2) = impPtr2(i,j,2) * 1.D-3
            enddo
            enddo
          else
            ! unit change: mm/s (same as kg s-1 m-2) -> m/s
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.D-3
            enddo
            enddo
          endif

        ! special treatment for latent heat flux
        else if (fieldStdName=="mean_lat_flx") then
          if (initFlag) then
            ! change the sign to comply with HYCOM convention 
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = -impPtr2(i,j,1) 
              impPtr2(i,j,2) = -impPtr2(i,j,2) 
            enddo
            enddo
          else
            ! change the sign to comply with HYCOM convention 
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = -impPtr2(i,j,1) 
            enddo
            enddo
          endif
       
        ! special treatment for sensible heat flux
        else if (fieldStdName=="mean_sens_flx") then
          if (initFlag) then
            ! change the sign to comply with HYCOM convention 
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = -impPtr2(i,j,1) 
              impPtr2(i,j,2) = -impPtr2(i,j,2) 
            enddo
            enddo
          else
           ! change the sign to comply with HYCOM convention 
           do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = -impPtr2(i,j,1) 
            enddo
            enddo
          endif
        
        ! special treatment for rivers water flux
        else if (fieldStdName=="water_flux_into_sea_water") then
          if (initFlag) then
            ! convert (kg m-2 s-1) to (m s-1) 
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.e-3
              impPtr2(i,j,2) = impPtr2(i,j,2) * 1.e-3
            enddo
            enddo
          else
           !  convert (kg m-2 s-1) to (m s-1) 
           do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.e-3
            enddo
            enddo
          endif
        ! special treatment for rivers  flux from ice
        else if (fieldStdName=="frozen_water_flux_into_sea_water") then
          if (initFlag) then
            ! convert (kg m-2 s-1) to (m s-1)
            do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.e-3
              impPtr2(i,j,2) = impPtr2(i,j,2) * 1.e-3
            enddo
            enddo
          else
           !  convert (kg m-2 s-1) to (m s-1) 
           do j=1,jja
            do i=1,ii
              impPtr2(i,j,1) = impPtr2(i,j,1) * 1.e-3
            enddo
            enddo
          endif
       endif
       
      else ! single-level imports
        do j=1,jja
        do i=1,ii
          impPtr(i,j) = farrayPtr(i,j)
        enddo
        enddo
      endif
      
    enddo

    ! Rotate the wind stress to local grid    
    if (cpl_taux .and. cpl_tauy) then
       do j=1,jja
          do i=1,ii
             imp_taux(i,j,1)= imp_taue(i,j,1) * cos(pang(i,j)) + imp_taun(i,j,1) * sin(pang(i,j))
             imp_taux(i,j,2)= imp_taue(i,j,2) * cos(pang(i,j)) + imp_taun(i,j,2) * sin(pang(i,j))
             
             imp_tauy(i,j,1)= imp_taun(i,j,1) * cos(pang(i,j)) - imp_taue(i,j,1) * sin(pang(i,j))
             imp_tauy(i,j,2)= imp_taun(i,j,2) * cos(pang(i,j)) - imp_taue(i,j,2) * sin(pang(i,j))

          enddo
       enddo
    endif


    ! transfer SEA-ICE imports into native HYCOM variables
    ! NOTE: For accurate conservation of salt, update_ocn_f should be set to true in CICE
    do j=1,jja
    do i=1,ii
      if (iceflg.ge.2 .and. icmflg.ne.3) then
        covice(i,j) = sic_import(i,j) !Sea Ice Concentration
        si_c  (i,j) = sic_import(i,j) !Sea Ice Concentration
        if (covice(i,j).gt.0.0) then
          flxice(i,j) = sifh_import(i,j) !Sea Ice Heat Flux Melting potential
          si_tx (i,j) =  sitx_import(i,j) !Sea Ice X-Stress into ocean
          si_ty (i,j) =  sity_import(i,j) !Sea Ice Y-Stress into ocean
          fswice(i,j) =  siqs_import(i,j) !Solar Heat Flux thru Ice to Ocean already in swflx
          sflice(i,j) =  sifs_import(i,j)*1.e3 !Ice Freezing/Melting Salt Flux
          wflice(i,j) =  sifw_import(i,j) !Ice Water Flux
          temice(i,j) =   sit_import(i,j) !Sea Ice Temperature
          si_t  (i,j) =   sit_import(i,j) !Sea Ice Temperature
          thkice(i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_h  (i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_u  (i,j) =   siu_import(i,j) !Sea Ice X-Velocity
          si_v  (i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
        else
          si_tx (i,j) = 0.0
          si_ty (i,j) = 0.0
          fswice(i,j) = 0.0
          flxice(i,j) = 0.0
          sflice(i,j) = 0.0
          wflice(i,j) = 0.0
          temice(i,j) = 0.0
          si_t  (i,j) = 0.0
          thkice(i,j) = 0.0
          si_h  (i,j) = 0.0
          si_u  (i,j) = 0.0
          si_v  (i,j) = 0.0
        endif !covice
      elseif (iceflg.ge.2 .and. icmflg.eq.3) then
        si_c(i,j) =  sic_import(i,j) !Sea Ice Concentration
        if (si_c(i,j).gt.0.0) then
          si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
          si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
          si_h (i,j) =   sih_import(i,j) !Sea Ice Thickness
          si_t (i,j) =   sit_import(i,j) !Sea Ice Temperature
          si_u (i,j) =   siu_import(i,j) !Sea Ice X-Velocity
          si_v (i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
        else
          si_tx(i,j) = 0.0
          si_ty(i,j) = 0.0
          si_h (i,j) = 0.0
          si_t (i,j) = 0.0
          si_u (i,j) = 0.0
          si_v (i,j) = 0.0
        endif !covice
      endif !iceflg>=2 (icmflg)
    enddo
    enddo

    ! clean-up
    deallocate(fieldList, fieldNameList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocations failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------
  subroutine HYCOM_GlueFieldsDataExport(glue, initFlag, rc)
    type(hycom_nuopc_glue_type), intent(inout)  :: glue
    logical                                     :: initFlag
    integer, intent(out), optional              :: rc
    
    integer                           :: fieldCount, iField, stat
    type(ESMF_Field), allocatable     :: fieldList(:)
    character(len=80), allocatable    :: fieldNameList(:)
    character(len=80)                 :: fieldName, fieldStdName, msg
    type(ESMF_Field)                  :: field, shadow
    type(ESMF_StateIntent_Flag)       :: stateIntent
    real(kind=ESMF_KIND_R8), pointer  :: farrayPtr(:,:)
    integer                           :: i,j
    logical                           :: isConnected
    
    real(kind=ESMF_KIND_R8) :: hfrz, t2f, tfrz, smxl, tmxl, ssfi
    real(kind=ESMF_KIND_R8) :: usur1, usur2, vsur1, vsur2, utot, vtot
    real(kind=ESMF_KIND_R8) :: ssh_n,ssh_s,ssh_e,ssh_w,dhdx,dhdy
    integer                 :: cplfrq
    
    if (present(rc)) rc = ESMF_SUCCESS

    ! access the members inside of the exportFields FieldBundle
    call ESMF_FieldBundleGet(glue%exportFields, fieldCount=fieldCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(fieldList(fieldCount), fieldNameList(fieldCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocations failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_FieldBundleGet(glue%exportFields, fieldList=fieldList, &
      fieldNameList=fieldNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !write(msg, *) 'Number of fields in export Fields', fieldCount
    !call ESMF_LogWrite("HYCOM_GlueFieldsDataExport(): "// &
    !  trim(msg), ESMF_LOGMSG_INFO)
    
    ! loop over all of the export Fields and fill in the data
    do iField=1, fieldCount
    
      field = fieldList(iField)
      fieldName = fieldNameList(iField)
    
      call ESMF_AttributeGet(field, name="StandardName", value=fieldStdName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      isConnected = .false.
      call ESMF_AttributeGet(field, name="HYCOM_IN_CESM_Connected", value=isConnected, defaultValue=.false., rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
      return ! bail out
      !write(msg, *) trim(fieldStdName), ' connected: ', isConnected
      !call ESMF_LogWrite("HYCOM_GlueFieldsDataExport(): "// &
      !  trim(msg), ESMF_LOGMSG_INFO)
      if(.not. isConnected) cycle
      
      !call ESMF_LogWrite("HYCOM_GlueFieldsDataExport(): "// &
      !  trim(fieldStdName)//" - "//trim(fieldName), ESMF_LOGMSG_INFO)
        
#ifdef WORKAROUND_HOLES
      ! look for this field in the shadow_exportFields
      call ESMF_FieldBundleGet(glue%shadow_exportFields, fieldName=fieldName, &
        field=shadow, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      ! pointer to the export field data
      call ESMF_FieldGet(shadow, farrayPtr=farrayPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#else
      ! pointer to the export field data
      call ESMF_FieldGet(field, farrayPtr=farrayPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif

      ! identify the exact export field by standard name and copy the data
      if (fieldStdName == "sea_surface_temperature") then
        do j=1,jja
        do i=1,ii
          farrayPtr(i,j) = 0.5*(temp(i,j,1,2)+temp(i,j,1,1))  ! construct SST [C]
         farrayPtr(i,j) = farrayPtr(i,j) + 273.15d0            ! [C] -> [K]
        enddo
        enddo
      elseif (fieldStdName == "upward_sea_ice_basal_available_heat_flux") then
      ! get coupling frequency in time steps
        cplfrq = nint( ocn_cpl_frq*(86400.0/baclin) )
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then
! ---     quantities for available freeze/melt heat flux
! ---     relax to tfrz with e-folding time of cplfrq time steps
! ---     assuming the effective surface layer thickness is hfrz
! ---     multiply by dpbl(i,j)/hfrz to get the actual e-folding time
              hfrz = min( thkfrz*onem, dpbl(i,j) )
              t2f  = (spcifh*hfrz)/(baclin*cplfrq*g)
              ! ---     average both available time steps, to avoid time splitting.
              !smxl = 0.5*(saln(i,j,1,2)+saln(i,j,1,1)) !!Alex calculated in mod_hycom.F
              !tmxl = 0.5*(temp(i,j,1,2)+temp(i,j,1,1)) !!Alex calculated in mod_hycom.F 
              smxl = sml(i,j)
              tmxl = tml(i,j)
              tfrz = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
              ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
        
              farrayPtr(i,j) = max(-1000.0,min(1000.0,ssfi))
              frzh     (i,j) = max(-1000.0,min(1000.0,ssfi)) ! > 0. freezing potential of flxice
           else
              farrayPtr(i,j) = 0.
              frzh     (i,j) = 0.                            ! > 0. freezing potential of flxice
           endif
        enddo
        enddo        
      elseif (fieldStdName == "sea_lev") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
!!Alex              farrayPtr(i,j) = srfhgt(i,j)/g   ! convert ssh to m
              farrayPtr(i,j) = sshm(i,j)/g   ! convert ssh to m
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "mixed_layer_depth") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
              farrayPtr(i,j) = dpbl(i,j) * qonem    ! convert mld to m
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      elseif (fieldStdName == "s_surf") then
        do j=1,jja
        do i=1,ii
!!Alex          farrayPtr(i,j) = 0.5d0 * (saln(i,j,1,2)+saln(i,j,1,1))
          farrayPtr(i,j) = sml(i,j)

        enddo
        enddo
      elseif (fieldStdName == "ocn_current_zonal") then
        do j=1,jja
           do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = umxl(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "ocn_current_merid") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = vmxl(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
        enddo
        enddo
      elseif (fieldStdName == "eastward_sea_surface_slope") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = dhde(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      elseif (fieldStdName == "northward_sea_surface_slope") then
        do j=1,jja
        do i=1,ii
           if (.not. initFlag) then 
           ! Now calculated in mod_hycom.F with accurate loop boundaries to get a bit for bit restart
           ! Results now differs from when the calculation was done here
              farrayPtr(i,j) = dhdn(i,j)
           else
              farrayPtr(i,j) = 0.
           endif
         enddo
        enddo
      endif
      
#ifdef WORKAROUND_HOLES
      ! Redist from shadow to export
      call ESMF_FieldRedist(shadow, field, &
        routehandle=glue%rh_shadow2export, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
#endif
      
    enddo
    
    ! clean-up
    deallocate(fieldList, fieldNameList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocations failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

end module
