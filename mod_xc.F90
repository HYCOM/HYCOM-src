      module mod_xc
      use mod_dimensions  !include 'dimensions.h'
      implicit none
!
! --- HYCOM communication interface.
! --- see README.src.mod_xc for more details.
!
      include 'unit_offset.h'
!
! --- tile dimensions and tile numbers (counting from 1), see xcspmd
      integer, public, save      :: ipr,  jpr,  ijpr, &
                                    mproc,nproc,mnproc, &
                                    mp_1st
#if defined(MPI)
!
! --- needed for some versions of mod_za
      integer, public, save      :: group_1st_in_row
#if defined(OCEANS2)
!
! --- two copies of HYCOM on distinct MPI tasks, master/slave via mod_pipe.
      integer, public, save      :: nocean,idp1_1,idp1_2
#endif
#endif
!
! --- region type (-1==unknown,
! ---               0==  closed/closed,
! ---               1==periodic/closed,
! ---               2==periodic/arctic,
! ---               3==periodic/fplane
! ---               4==  closed/fplane)
      integer, public, save      :: nreg
!
! --- timers on, usually and default .true.
      logical, public, save      :: timer_on=.true.
!
! --- fill value for land, usually 0.0
      real,    public, save      :: vland
      real*4,  public, save      :: vland4  !xcget4 only
!
! --- xctilr halo options
      integer, public, parameter :: halo_ps=1, halo_pv=11, &
                                    halo_qs=2, halo_qv=12, &
                                    halo_us=3, halo_uv=13, &
                                    halo_vs=4, halo_vv=14
!
! --- xcsync stdout flushing options
      logical, public, parameter :: flush_lp=.true., &
                                    no_flush=.false.
!
! --- generic subroutine names
      interface xcmaxr
         module procedure xcmaxr_0  ! rank 0 array (i.e. scalar)
         module procedure xcmaxr_1  ! rank 1 array
         module procedure xcmaxr_0o ! rank 0 array, old interface
         module procedure xcmaxr_1o ! rank 1 array, old interface
      end interface

      interface xcminr
         module procedure xcminr_0  ! rank 0 array (i.e. scalar)
         module procedure xcminr_1  ! rank 1 array
         module procedure xcminr_0o ! rank 0 array, old interface
         module procedure xcminr_1o ! rank 1 array, old interface
      end interface

      interface xcsumr
         module procedure xcsumr_0  ! rank 0 array (i.e. scalar)
         module procedure xcsumr_1  ! rank 1 array
      end interface
#if defined(USE_ESMF4) || defined(ESPC_COUPLE)
!
! --- public data structures for ESMF, see xcspmd
      integer, public,  save :: deBlockList(2,2,iqr*jqr)
#endif
!
! --- private subroutines
      private xciget_sm
!
! --- private timer variables, see xctmri
      character*6, private, dimension(97), save :: cc
      integer,     private,                save :: nxc
      integer,     private, dimension(97), save :: nc,nca,nci
      real*8,      private, dimension(97), save :: tc,t0
      real*8,      private, dimension(2),  save :: tcxc,tcxl
#if defined(MPI) || defined(SHMEM)
!
! --- private message passing data structures, see xcspmd
#if defined(RELO)
      integer, private, save, allocatable, dimension(:,:) :: mpe_i
      integer, private, save, allocatable, dimension(:)   :: npe_j
#else
      integer, private, save :: mpe_i(0:itdm+1,0:jqr),npe_j(0:jtdm+1)
#endif
      integer, private, save :: idproc( 0: iqr+1,0:jqr+1), &
                                idproc1(0:ijqr+1),idhalo(2), &
                                i0_pe(iqr,jqr),ii_pe(iqr,jqr), &
                                j0_pe(iqr,jqr),jj_pe(iqr,jqr), &
                                mpe_1(jqr), &
                                mpe_e(jqr)
      integer, private, save :: i1sum(iqr,jqr),iisum(iqr,jqr)
      integer, private, save :: m0_top,i0_st(iqr),ii_st(iqr), &
                                mm_top,i0_gt(iqr),ii_gt(iqr), &
                                m0_bot,i0_sb(iqr),ii_sb(iqr), &
                                mm_bot,i0_gb(iqr),ii_gb(iqr)
      integer, private, save :: null_tile
#endif
#if defined(MPI)
      integer, private, save :: mpi_comm_hycom
#endif
#if defined(MPI) || defined(SHMEM)
!
! --- private message passing buffers, used by more than one routine
#if defined(RELO)
      real,   save, private, allocatable, dimension(:,:) :: &
         al    ! xcaget and xcaput
      real,   save, private, allocatable, dimension(:) :: &
         at    ! xcaget and xcaput
      real*4, save, private, allocatable, dimension(:,:) :: &
         al4,   & ! xcaget4 and xcaput4
         alt4  ! xcaget4 and xcaput4
      real*4, save, private, allocatable, dimension(:) :: &
         at4   ! xcaget4 and xcaput4
      integer, save, private :: &
         mxsum
      real*8, save, private, allocatable, dimension(:) :: &
        sum8t,  & ! xcsum and xcsumj
        sum8j  ! xcsum and xcsumj
      real*8, save, private :: &
        sum8s
#else
      real,   save, private, dimension(itdm,jdm) :: &
         al    ! xcaget and xcaput
      real,   save, private, dimension(idm*jdm) :: &
         at    ! xcaget and xcaput
      real*4, save, private, dimension(itdm,jdm) :: &
         al4   ! xcaget4 and xcaput4
      real*4, save, private, dimension(idm*jdm,jpr) :: &
         alt4  ! xcaget4 and xcaput4
      real*4, save, private, dimension(idm*jdm) :: &
         at4   ! xcaget4 and xcaput4
      integer, private, parameter :: mxsum=(idm+4*nbdy)/(2*nbdy+1)
      real*8, save, private, dimension(mxsum*jdm) :: &
        sum8t  ! xcsum and xcsumj
      real*8, save, private, dimension(jdm) :: &
        sum8j  ! xcsum and xcsumj
      real*8, save, private :: &
        sum8s  ! xcsum and xcsumj
#endif
      integer, private, parameter :: nmax=1024
      real,   save, private, dimension(nmax) :: &
        b,c    ! xcastr, xcmaxr, xcminr, xcsumr
#endif
!
! --- actual module subroutines
      contains
#if defined(MPI) || defined(SHMEM)
# include "mod_xc_mp.h"
#else
# include "mod_xc_sm.h"
#endif
      end module mod_xc
