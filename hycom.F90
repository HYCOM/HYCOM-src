#if defined(USE_ESMF4)
      program hycom
!
! --- ESMF driver for stand-alone HYCOM ocean model
!
      use ESMF_Mod
!     use mod_hycom, only : OCEAN_SetServices => HYCOM_SetServices
      use mod_hycom, only : end_of_run, &
                            OCEAN_SetServices => HYCOM_SetServices
!
      implicit none
!
! --- Local variables
!
! --- Gridded Components
      type(ESMF_GridComp) :: oceanGridComp
!
! --- States, Virtual Machines, and Layouts
      type(ESMF_VM) :: worldVM
      type(ESMF_State) :: oceanImpState, oceanExpState
      integer :: petCount, localPet, split
!
! --- Calendars and clocks
      type(ESMF_Clock) :: worldClock
      type(ESMF_Clock) :: oceanClock
!
! --- Return codes for error checks
      integer :: rc
!
! --- Miscellaneous
      integer :: i
!
!-------------------------------------------------------------------------------
!  Initialize the ESMF Framework
!-------------------------------------------------------------------------------
!
! --- Set default calendar and log type; get world VM
      call ESMF_Initialize(defaultCalendar=ESMF_CAL_GREGORIAN, &
                           defaultLogType=ESMF_LOG_SINGLE, &
                           vm=worldVM, rc=rc)
      if (rc .ne. ESMF_SUCCESS) stop 99
!
! --- Get VM info
      call ESMF_VMGet(worldVM, petCount=petCount, localPET=localPet, &
                      rc=rc)
      if (ESMF_LogMsgFoundError(rc, "ESMF_VMGet failed", rc))  &
         goto 10
!
!-------------------------------------------------------------------------------
! --- Create section
!-------------------------------------------------------------------------------
!
! --- Create the OCEAN gridded component
      oceanGridComp = ESMF_GridCompCreate(vm=worldVM, &
                      name="OCEAN Gridded Component", &
                      gridCompType=ESMF_OCEAN, &
                      rc=rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN GridCompCreate failed", rc)) &
         goto 10
!
! --- Create empty OCEAN import/export states
      oceanImpState = ESMF_StateCreate(stateName="OCEAN Import State", &
                                       stateType=ESMF_STATE_IMPORT, &
                                       rc=rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN ImpState Create failed", rc)) &
         goto 10
      oceanExpState = ESMF_StateCreate(stateName="OCEAN Export State", &
                                       stateType=ESMF_STATE_EXPORT, &
                                       rc=rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN ExpState Create failed", rc)) &
         goto 10
!
!-------------------------------------------------------------------------------
! --- Register section
!-------------------------------------------------------------------------------
!
! --- Register the OCEAN gridded component
      call ESMF_GridCompSetServices(oceanGridComp, &
                                    OCEAN_SetServices, rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN Registration failed", rc)) &
         goto 10
!
!-------------------------------------------------------------------------------
! --- Initalize Section
!-------------------------------------------------------------------------------
!
! --- Initialize OCEAN gridded component
      call ESMF_GridCompInitialize(gridComp=oceanGridComp, &
                                   importState=oceanImpState, &
                                   exportState=oceanExpState, &
                                   clock=worldClock, &
                                   phase=ESMF_SINGLEPHASE, &
                                   blockingflag=ESMF_NONBLOCKING, &
                                   rc=rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN Initialize failed", rc)) &
         goto 10
!c
!c --- Get copy of OCEAN clock
!      call ESMF_GridCompGet(oceanGridComp, clock=oceanClock, rc=rc)
!c
!c --- Initialize WORLD clock using OCEAN clock
!      worldClock = ESMF_ClockCreate(clock=oceanClock, rc=rc)
!
!-------------------------------------------------------------------------------
! --- Run Section
!-------------------------------------------------------------------------------
!
      do !until end of run
        call ESMF_GridCompRun(gridComp=oceanGridComp, &
                              importState=oceanImpState, &
                              exportState=oceanExpState, &
                              clock=worldClock, &
                              phase=ESMF_SINGLEPHASE, &
                              blockingflag=ESMF_NONBLOCKING, &
                              rc=rc)
!
! ---   use end_of_run, rather than a ESMF Clock
        if     (end_of_run) then
          exit
        endif
      enddo
!
!-------------------------------------------------------------------------------
!  Finalize Section
!-------------------------------------------------------------------------------
!
! --- Finalize OCEAN gridded component
      call ESMF_GridCompFinalize(gridComp=oceanGridComp, &
                                 importState=oceanImpState, &
                                 exportState=oceanExpState, &
                                 clock=worldClock, &
                                 phase=ESMF_SINGLEPHASE, &
                                 blockingflag=ESMF_NONBLOCKING, &
                                 rc=rc)
      if (ESMF_LogMsgFoundError(rc, "OCEAN Finalize failed", rc))  &
         goto 10
!
10    continue
      call ESMF_VMBarrier(worldVM)
      call ESMF_Finalize(rc=rc)
!
      stop
      end program hycom
#else
      program hycom
!
! --- Non-ESMF driver for stand-alone HYCOM ocean model
!
      use mod_hycom, only : end_of_run, &
                            HYCOM_Init, &
                            HYCOM_Run, &
                            HYCOM_Final
!
      implicit none
!
! --- Initialize HYCOM.
      call HYCOM_Init

! --- Run HYCOM.
      do !until end of run
        call HYCOM_Run
        if     (end_of_run) then
          exit
        endif
      enddo
!
! --- Finalize HYCOM.
      call HYCOM_Final
!
      stop
      end program hycom
#endif
