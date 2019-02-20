#if defined(NAN2003)
      logical function hycom_isnaninf(a)
      use ieee_arithmetic, only : ieee_is_finite
      implicit none
!
      real a
!
!**********
!*
! 1)  return .true. if a is NaN or +Inf or -Inf.
!
! 2)  version with ieee_arithmetic intrinsic module
!*
!**********
!
      hycom_isnaninf = .not. ieee_is_finite(a)
      end function hycom_isnaninf
#else
      logical function hycom_isnaninf(a)
      implicit none
!
      real a
!
!**********
!*
! 1)  return .true. if a is NaN or +Inf or -Inf.
!*
!**********
!
      hycom_isnaninf = .not. (a.ge.-huge(a) .and. a.le.huge(a))
      end
#endif
!
!> Revision history:
!>
!> Jan. 2014 - ieee_arithmetic intrinsic module version via /* NAN2003 */ macro
