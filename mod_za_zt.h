!
!-----------------------------------------------------------------------
!
!     machine dependent I/O routines.
!     per tile version, contained in mod_za.
!
!     author:  Alan J. Wallcraft,  NRL.
!
!-----------------------------------------------------------------------
!
      subroutine ztiopf(cfile,cstat, iaunit)
      implicit none
!
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cfile,cstat
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for opening a file for array i/o.
!
!     must call zaiost before first call to ztiopf.
!     see also 'ztiopn' and 'ztiope'.
!
!  2) the filename is taken from 'cfile'.
!
!     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     cstat indicates the file type, it can be 'scratch', 'old', or
!      'new'.
!     all i/o to iaunit must be performed by ztiowr.
!      arrays passed to these routines must conform to 'h'.
!     the file should be closed using ztiocl.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
! --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
!
      integer   ios,nrecl
      character cact*9
#if defined(TIMER)
!
      call xctmr0(16)
#endif
!
!     test file state.
!
      if     (iarec(iaunit).ne.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        call xchalt('(ztiopf)')
               stop '(ztiopf)'
      endif
!     write(lp,*) 'ztiopf - iaunit = ',iaunit
!     call flush(lp)
!
!     open file.
!
      inquire(iolength=nrecl) w(1:ii*jj)
!
      if     (cstat.eq.'OLD' .or. &
              cstat.eq.'old'     ) then
        cact = 'READ'
      elseif (cstat.eq.'NEW' .or. &
              cstat.eq.'new'     ) then
        cact = 'WRITE'
      else
        cact = 'READWRITE'
      endif
!
#if defined(X1)
      call asnunit(iaunit+uaoff,'-F event,cachea:4096:4:2 -B on',ios)
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios   = ',ios
        call flush(lp)
        call xchalt('(ztiopf)')
               stop '(ztiopf)'
      endif !ios
#endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+uaoff,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+uaoff,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit,trim(cfile)
        write(lp,*) 'ios   = ',ios
        call flush(lp)
        call xchalt('(ztiopf)')
               stop '(ztiopf)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile,  &
           form='unformatted', status=cstat, &
           access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        call flush(lp)
        call xchalt('(ztiopf)')
               stop '(ztiopf)'
      endif !ios
      iarec(iaunit) = 0
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in ztiopf -  array I/O unit ', &
         i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in ztiopf -  can''t asnunit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in ztiopf -  can''t open unit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
      end subroutine ztiopf

      subroutine ztiocl(iaunit)
      implicit none
!
      integer, intent(in)    :: iaunit
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for array i/o file closing.
!
!     must call ztiopn for this array unit before calling ztiocl.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!*
!**********
!
      integer ios
#if defined(TIMER)
!
      call xctmr0(16)
#endif
!
!     write(lp,*) 'ztiocl - iaunit = ',iaunit
!     call flush(lp)
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        call xchalt('(ztiocl)')
               stop '(ztiocl)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        close(unit=iaunit+uaoff, status='keep')
#if defined(T3E) || defined(YMP) || defined(X1)
        call asnunit(iaunit+uaoff,'-R',ios)
#endif
      endif
      iarec(iaunit) = -1
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in ztiocl -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
      end subroutine ztiocl

      subroutine ztiowr3(h, l, mask,lmask, hmin,hmax, iaunit, lreal4)
      implicit none
!
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l), &
               intent(inout) :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l), &
               intent(inout) :: h
#endif
!
!**********
!*
!  1) machine specific routine for 3-d array writing.
!
!     must call ztiopn for this array unit before calling ztiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to ztiopn.
!
!  4) hmin,hmax are returned as the minimum and maximum value in the array.
!     if lmask==.true. the range is only where mask.ne.0, with all other
!     values output as 2.0**100.
!
!  5) If lreal4==.true. then h is overwritten on exit with real*4 version
!     of the same array.  This is typically used for reproducability on
!     restart.
!*
!**********
!
!     this version just calls ztiowr l times.
!
      integer k
!
      do k= 1,l
        call ztiowr(h(1-nbdy,1-nbdy,k), mask,lmask, &
                    hmin(k),hmax(k), iaunit, lreal4)
      enddo
      return
      end subroutine ztiowr3

      subroutine ztiowr(h, mask,lmask, hmin,hmax,  iaunit, lreal4)
      implicit none
!
      logical, intent(in)    :: lmask,lreal4
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: h
#endif
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for array writing.
!
!     must call ztiopn for this array unit before calling ztiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to ztiopn.
!
!  4) hmin,hmax are returned as the minimum and maximum value in the array.
!     if lmask==.true. the range is only where mask.ne.0, with all other
!     values output as 2.0**100.
!
!  5) If lreal4==.true. then h is overwritten on exit with real*4 version
!     of the same array.  This is typically used for reproducability on
!     restart.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
! --- n2drec = size of output 2-d array, multiple of 4096
      real*4     spval
      parameter (spval=2.0**100)
!
      character cfile*256
      integer   ios, i,j
#if defined(TIMER)
!
      call xctmr0(18)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        call xchalt('(ztiowr)')
               stop '(ztiowr)'
      endif
!
      if     (lreal4) then
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                w(i+(j-1)*ii) = h(i,j)
                wminy(j) = min( wminy(j), w(i+(j-1)*ii) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*ii) )
              else
                w(i+(j-1)*ii) = spval
              endif
#if defined(REAL4)
! ---         h(i,j) = w(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = w(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,ii
              w(i+(j-1)*ii) = h(i,j)
              if     (w(i+(j-1)*ii).ne.spval) then
                wminy(j) = min( wminy(j), w(i+(j-1)*ii) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*ii) )
              endif
#if defined(REAL4)
! ---         h(i,j) = w(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = w(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        endif
      else
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                w(i+(j-1)*ii) = h(i,j)
                wminy(j) = min( wminy(j), w(i+(j-1)*ii) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*ii) )
              else
                w(i+(j-1)*ii) = spval
              endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,ii
              w(i+(j-1)*ii) = h(i,j)
              if     (w(i+(j-1)*ii).ne.spval) then
                wminy(j) = min( wminy(j), w(i+(j-1)*ii) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*ii) )
              endif
            enddo
          enddo
        endif
      endif
      hmin = minval(wminy(1:jj))
      hmax = maxval(wmaxy(1:jj))
      iarec(iaunit) = iarec(iaunit) + 1
      call ztiowrd(w,ii*jj, iaunit+uaoff,iarec(iaunit),ios)
      if     (ios.ne.0) then
        write(lp,9100) iarec(iaunit),iaunit
        call flush(lp)
        cfile = ' '
        inquire(unit=iaunit+uaoff,name=cfile)
        write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
        call flush(lp)
        call xchalt('(ztiowr)')
               stop '(ztiowr)'
      endif !ios
#if defined(TIMER)
!
      call xctmr1(18)
#endif
      return
!
 9000 format(/ /10x,'error in ztiowr -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in ztiowr -  can''t write record', &
         i4,' on array I/O unit ',i3,'.'/ /)
      end subroutine ztiowr

      subroutine ztiowrd(a,n, iunit,irec,ios)
      implicit none
!
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(inout) :: a(n)  !needed if zaio_endian is called
!
!**********
!*
! 1)  direct access write a single record.
!
! 2)  expressed as a subroutine because i/o with 
!     implied do loops can be slow on some machines.
!*
!**********
!
#if defined(ENDIAN_IO)
      call zaio_endian(a,n)  ! overwrites a
#endif
      write(unit=iunit, rec=irec, iostat=ios) a
      return
      end subroutine ztiowrd
