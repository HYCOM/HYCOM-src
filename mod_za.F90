      module mod_za
      use mod_xc  ! HYCOM communication API
!
      implicit none
!
! --- HYCOM I/O interface.
!
! --- See README.src.mod_za for more details.
!
#if defined(MPI) && ! defined(SERIAL_IO)
      integer, save, private :: &
        file_info_zaiord,file_count_zaiord, &
        file_info_zaiowr,file_count_zaiowr
#else
      private zaiordd,zaiowrd
#endif
      private ztiowrd
!
      integer, private, parameter :: uaoff = 1000 + uoff  !array I/O unit offset
#if defined(RELO)
!
! --- initialized in zaiost.
      integer, private, save :: n2drec
      real*4,  private, save, allocatable, dimension(:) :: &
                 wminy,wmaxy,htmp
#else
!
      integer, parameter :: n2drec=((itdm*jtdm+4095)/4096)*4096
      real*4,  private, save, dimension(jtdm)    :: wminy,wmaxy
      real*4,  private, save, dimension(idm*jdm) :: htmp
#endif
!
      real*4,  private, save, allocatable, dimension(:) :: w
!
      contains
#if ! defined(MPI) && ! defined(SHMEM)
#       include "mod_za_sm.h"
#elif defined(SERIAL_IO)
#       include "mod_za_mp1.h"
#else
#       include "mod_za_mp.h"
#endif
#      include "mod_za_zt.h"
      end module mod_za

#if defined(ENDIAN_IO_F90)  /* see mach_i.c for new C version */
      subroutine zaio_endian(a,n)
      implicit none
!
      integer,         intent(in)    :: n
      integer(kind=4), intent(inout) :: a(n)  ! 4-bytes
!
!**********
!*
! 1)  swap the endian-ness of the array.
!
! 2)  assumes integer(kind=1) and integer(kind=4) ocupy one and four
!     bytes respectively.
!*
!**********
!
      integer k
!
      integer(kind=4) ii4,   io4     ! 4-bytes
      common/czioxe/  ii4,   io4     ! helps prevent unwanted optimization
      save  /czioxe/
!
      integer(kind=1) ii1(4),io1(4)  ! 1-byte
      equivalence    (ii4,ii1(1)), (io4,io1(1))  ! non-standard f90
!
      do k= 1,n
        ii4 = a(k)
        io1(1) = ii1(4)
        io1(2) = ii1(3)
        io1(3) = ii1(2)
        io1(4) = ii1(1)
        a(k) = io4
      enddo
      return
      end subroutine zaio_endian
#endif /* ENDIAN_IO_F90 */
