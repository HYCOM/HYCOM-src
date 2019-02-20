module ocn_comp_mct

! !USES:

  use shr_sys_mod
  use shr_kind_mod     , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS
  use shr_file_mod     , only: shr_file_getunit, shr_file_getlogunit, shr_file_getloglevel, &
                               shr_file_setlogunit, shr_file_setloglevel, shr_file_setio, &
                               shr_file_freeunit
  use shr_mpi_mod      , only: shr_mpi_bcast
  use mct_mod
  use esmf

  use dead_data_mod
  use dead_mct_mod
  use dead_mod

  use seq_cdata_mod
  use seq_infodata_mod
  use seq_timemgr_mod
  use seq_comm_mct     , only: seq_comm_inst, seq_comm_name, seq_comm_suffix
  use seq_flds_mod     , only: flds_d2x => seq_flds_o2x_fields, &
                               flds_x2d => seq_flds_x2o_fields

!
! !PUBLIC TYPES:
  implicit none
  save
  private ! except

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public :: ocn_init_mct
  public :: ocn_run_mct
  public :: ocn_final_mct

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_init_mct
!
! !DESCRIPTION:
!     initialize dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

  subroutine ocn_init_mct( EClock, cdata, x2d, d2x, NLFilename )

! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Clock)         , intent(inout) :: EClock
    type(seq_cdata)          , intent(inout) :: cdata
    type(mct_aVect)          , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in) :: NLFilename ! Namelist filename

!EOP

  print *, "ALERT: OCN_INIT_MCT should not be called"


end subroutine ocn_init_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_run_mct
!
! !DESCRIPTION:
!     run method for dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------

subroutine ocn_run_mct( EClock, cdata, x2d, d2x)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock)            ,intent(inout) :: EClock
   type(seq_cdata)             ,intent(inout) :: cdata
   type(mct_aVect)             ,intent(inout) :: x2d        ! driver-> dead
   type(mct_aVect)             ,intent(inout) :: d2x        ! dead  -> driver

!EOP


end subroutine ocn_run_mct

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: ocn_final_mct
!
! !DESCRIPTION:
!     finalize method for dead ocn model
!
! !REVISION HISTORY:
!
! !INTERFACE: ------------------------------------------------------------------
!
subroutine ocn_final_mct(EClock, cdata, x2d, d2x)

    implicit none

    !----- arguments -----

    type(ESMF_Clock)            ,intent(inout) :: EClock     ! clock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

!EOP

      
 end subroutine ocn_final_mct
!===============================================================================

end module ocn_comp_mct
