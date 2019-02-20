!
! --- machine-specific Fortran routines
!
      subroutine machine()
!
! --- always called once at the start of the program.
!
#if defined(SGI)
      call zunder()  ! C-wrapper to flush underflow to zero on R10000
#endif
      end
#if defined(AIX)
      subroutine flush(iunit)
      implicit none
      integer iunit
!
! --- wrapper for flush system call under AIX.
!
      integer*4 iunit4
!
      iunit4=iunit
      call flush_(iunit4)
      return
      end
#endif /* AIX */
#if defined(X1)
      subroutine x1flush(iunit)
      implicit none
      integer iunit
!
! --- wrapper for flush system call on the Cray X1.
!
      integer ierr
!
      call FLUSH(iunit,ierr)
      return
      end
#endif /* X1 */
#if defined(IFC)
      subroutine flush(iunit)
      implicit none
      integer iunit
!
! --- disable the flush system call under Intel's IFC compiler.
!
      return
      end
#endif /* IFC */
#if defined(SUN)
      subroutine ieee_retrospective()
!
!     dummy routine to turn off ieee warning messages on a Sun.
!
      end
#endif /* SUN */
#if defined(T3E) || defined(YMP) || defined(X1)
      subroutine getenv(cname, cvalue)
      implicit none
!
      character*(*) cname,cvalue
!
!     this subroutine provides getenv functionality
!     on the t3e, using pxfgetenv.
!
      integer iname,ivalue,ierr
!
      iname = 0
      ierr  = 0
      call pxfgetenv(cname,iname, cvalue,ivalue, ierr)
      if     (ierr.ne.0) then
        cvalue = ' '
      endif
      return
!     end of getenv.
      end
#endif /* T3E || YMP || X1 */
#if defined(NAGFOR)
      subroutine flush(iunit)
      implicit none
      integer iunit
!
! --- a wrapper for Fortran 2003's FLUSH statement
!
      flush(iunit)
      return
      end
      subroutine getenv(cname, cvalue)
      implicit none
!
      character*(*) cname,cvalue
!
! --- a wrapper for Fortran 2003's GET_ENVIRONMENT_VARIABLE
!
      call get_environment_variable(cname,cvalue)
      return
      end
#endif /* NAGFOR */
