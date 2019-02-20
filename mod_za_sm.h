!
!-----------------------------------------------------------------------
!
!     machine dependent I/O routines.
!     single processor version, contained in mod_za.
!
!     author:  Alan J. Wallcraft,  NRL.
!
!-----------------------------------------------------------------------
!
      subroutine zagetc(cline,ios, iunit)
      implicit none
!
      character*80, intent(out)   :: cline
      integer,      intent(out)   :: ios
      integer,      intent(in)    :: iunit
!
!**********
!*
!  1) machine specific routine for reading one text line from a file.
!*
!**********
!
      read(iunit,'(a)',iostat=ios) cline
      return
      end subroutine zagetc

      subroutine zaiopn(cstat, iaunit)
      implicit none
!
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cstat
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for opening a file for array i/o.
!
!     must call zaiost before first call to zaiopn.
!     see also 'zaiope' and 'zaiopf'.
!
!  2) the filename is taken from the environment variable FORxxxA,
!       where xxx = iaunit, with default fort.xxxa.
!
!     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     cstat indicates the file type, it can be 'scratch', 'old', or
!      'new'.
!     all i/o to iaunit must be performed by zaiord and zaiowr.
!     the file should be closed using zaiocl.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
!
      integer   ios,nrecl
      character cfile*256,cenv*7
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
        stop '(zaiopn)'
      endif
!
!     get filename.
!
      write(cenv,"('FOR',i3.3,'A')") iaunit
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(cfile,"('fort.',i3.3,'a')") iaunit
      endif
!     write(lp,*) 'zaiopn - iaunit = ',iaunit
!     call flush(lp)
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
!
!     open file.
!
      inquire(iolength=nrecl) w
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
#if defined(X1)
      call asnunit(iaunit+uaoff,'-F event,cachea:4096:4:2 -B on',ios)
      if     (ios.ne.0) then
        write(lp,9050) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        stop '(zaiopn)'
      endif !ios
#endif
#if defined(YMP)
      if     (mod(nrecl,16384).eq.0 .and. nrecl.gt.16384*4) then
       call asnunit(iaunit+uaoff,'-F syscall -N ieee',ios)
      else
        call asnunit(iaunit+uaoff,'-F cachea:8:16:2 -N ieee',ios)
      endif
      if     (ios.ne.0) then
        write(lp,9050) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        stop '(zaiopn)'
      endif !ios
#endif
      if     (cstat.eq.'scratch' .or. &
              cstat.eq.'SCRATCH'     ) then
        open(unit=iaunit+uaoff,              &
             form='unformatted', status='scratch', &
             access='direct', recl=nrecl, action=cact, iostat=ios)
      else
        open(unit=iaunit+uaoff, file=cfile,  &
             form='unformatted', status=cstat, &
             access='direct', recl=nrecl, action=cact, iostat=ios)
      endif
      if     (ios.ne.0) then
        write(lp,9100) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        if     (cstat.eq.'scratch' .or. &
                cstat.eq.'SCRATCH'     ) then
          write(lp,'(a)')  "status='SCRATCH'"
        else
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
        endif
        call flush(lp)
        stop '(zaiopn)'
      endif !ios
      iarec(iaunit) = 0
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiopn -  array I/O unit ', &
         i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiopn -  can''t asnunit ',i3, &
         ', for array I/O.'/ /)
#endif
 9100 format(/ /10x,'error in zaiopn -  can''t open unit ',i3, &
         ', for array I/O.'/ /)
      end subroutine zaiopn

      subroutine zaiope(cenv,cstat, iaunit)
      implicit none
!
      integer,       intent(in)    :: iaunit
      character*(*), intent(in)    :: cenv,cstat
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for opening a file for array i/o.
!
!     must call zaiost before first call to zaiope.
!     see also 'zaiopn' and 'zaiopf'.
!
!  2) the filename is taken from environment variable 'cenv'.
!
!     array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     cstat indicates the file type, it can be 'scratch', 'old', or
!      'new'.
!     all i/o to iaunit must be performed by zaiord and zaiowr.
!      arrays passed to these routines must conform to 'h'.
!     the file should be closed using zaiocl.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
!
      integer   ios,nrecl
      character cfile*256
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
        stop '(zaiope)'
      endif
!
!     get filename.
!
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(lp,9300) trim(cenv)
        write(lp,*) 'iaunit = ',iaunit
        call flush(lp)
        stop '(zaiope)'
      endif
!     write(lp,*) 'zaiope - iaunit = ',iaunit
!     call flush(lp)
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
!
!     open file.
!
      inquire(iolength=nrecl) w
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
        write(lp,*) 'ios = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        stop '(zaiope)'
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
        write(lp,*) 'ios = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        stop '(zaiope)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile,  &
           form='unformatted', status=cstat, &
           access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        stop '(zaiope)'
      endif !ios
      iarec(iaunit) = 0
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiope -  array I/O unit ', &
         i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiope -  can''t asnunit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiope -  can''t open unit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
 9300 format(/ /10x,'error in zaiope -  environment variable ',a, &
         ' not defined'/ /)
      end subroutine zaiope

      subroutine zaiopf(cfile,cstat, iaunit)
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
!     must call zaiost before first call to zaiopf.
!     see also 'zaiopn' and 'zaiope'.
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
!     all i/o to iaunit must be performed by zaiord and zaiowr.
!      arrays passed to these routines must conform to 'h'.
!     the file should be closed using zaiocl.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
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
        stop '(zaiopf)'
      endif
!     write(lp,*) 'zaiopf - iaunit = ',iaunit
!     call flush(lp)
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
#if defined(TIMER)
        call xctmr1(16)
#endif
        return
      endif
!
!     open file.
!
      inquire(iolength=nrecl) w
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
        stop '(zaiopf)'
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
        stop '(zaiopf)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile,  &
           form='unformatted', status=cstat, &
           access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        call flush(lp)
        stop '(zaiopf)'
      endif !ios
      iarec(iaunit) = 0
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiopf -  array I/O unit ', &
         i3,' is not marked as available.'/ /)
#if defined(YMP) || defined(X1)
 9050 format(/ /10x,'error in zaiopf -  can''t asnunit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
#endif
 9100 format(/ /10x,'error in zaiopf -  can''t open unit ',i3, &
         ', for array I/O.' / &
         10x,'cfile = ',a/ /)
      end subroutine zaiopf

      subroutine zaiopi(lopen, iaunit)
      implicit none
!
      logical, intent(out)   :: lopen
      integer, intent(in)    :: iaunit
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) is an array i/o unit open?
!
!  2) must call zaiost before first call to zaiopi.
!*
!**********
!
      lopen = iarec(iaunit).ne.-1
      return
      end subroutine zaiopi

      subroutine zaiost
      implicit none
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for initializing array i/o.
!
!  2) see also zaiopn, zaiord, zaiowr, and zaiocl.
!*
!**********
!
      integer i
! 
      write(lp,'(/a/)') &
          'zaiost - Array I/O is Fortran DA I/O'
!
      do i= 1,999
        iarec(i) = -1
      enddo
#if defined(RELO)
!
! --- n2drec = size of output 2-d array, multiple of 4096
      n2drec = ((itdm*jtdm+4095)/4096)*4096
!
      allocate( w(n2drec),wminy(jtdm),wmaxy(jtdm) )
      call mem_stat_add( (n2drec+2*jtdm)/2) !real*4, so /2
#endif
#if defined(TIMER)
!
!     initialize timers.
!
      call xctmrn(16,'zaio**')
      call xctmrn(17,'zaiord')
      call xctmrn(18,'zaiowr')
#endif
      return
      end subroutine zaiost

      subroutine zaiocl(iaunit)
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
!     must call zaiopn for this array unit before calling zaiocl.
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
!     write(lp,*) 'zaiocl - iaunit = ',iaunit
!     call flush(lp)
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiocl)'
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
 9000 format(/ /10x,'error in zaiocl -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
      end subroutine zaiocl

      subroutine zaiofl(iaunit)
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
!  1) machine specific routine for array i/o buffer flushing.
!
!     must call zaiopn for this array unit before calling zaiocl.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!*
!**********
!
      integer   irlen
      character cfile*256
#if defined(TIMER)
!
      call xctmr0(16)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiofl)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        inquire(unit=iaunit+uaoff, name=cfile, recl=irlen)
        close(  unit=iaunit+uaoff, status='keep')
        open(   unit=iaunit+uaoff, file=cfile, form='unformatted',  &
                access='direct', recl=irlen)
      endif
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiofl -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
      end subroutine zaiofl

      subroutine zaioiq(iaunit, irec)
      implicit none
!
      integer, intent(in)    :: iaunit
      integer, intent(out)   :: irec
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for array i/o inquiry.
!
!  2) returns the number of records processed, or -1 for a closed file.
!*
!**********
!
      irec = iarec(iaunit)
      return
      end subroutine zaioiq

      subroutine zaiorw(iaunit)
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
!  1) machine specific routine for array i/o file rewinding.
!
!     must call zaiopn for this array unit before calling zaiocl.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!*
!**********
#if defined(TIMER)
!
      call xctmr0(16)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiorw)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = 0
      endif
!     write(lp,*) 'zaiorw - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiorw -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
      end subroutine zaiorw

      subroutine zaiord3(h, l, mask,lmask, hmin,hmax,  iaunit)
      implicit none
!
      logical, intent(in)    :: lmask
      integer, intent(in)    :: l,iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin(l),hmax(l)
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l), &
               intent(out)   :: h
#else
      real,    intent(out)   :: hmin(l),hmax(l)
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l), &
               intent(out)   :: h
#endif
!
!**********
!*
!  1) machine specific routine for 3-d array reading.
!
!     must call zaiopn for this array unit before calling zaiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
!
!  4) hmin,hmax are returned as the minimum and maximum value in the 
!     array, ignoring array elements set to 2.0**100.  
!     if lmask==.true. the range is calculated only where mask.ne.0,
!     with all other values unchanged in h on exit.  It is then an
!     error if mask.ne.0 anywhere the input is 2.0**100.
!*
!**********
!
!     this version just calls zaiord l times.
!
      integer k
!
      do k= 1,l
        call zaiord(h(1-nbdy,1-nbdy,k), mask,lmask, &
                    hmin(k),hmax(k), iaunit)
      enddo
!
      return
      end subroutine zaiord3

      subroutine zaiord(h, mask,lmask, hmin,hmax,  iaunit)
      implicit none
!
      logical, intent(in)    :: lmask
      integer, intent(in)    :: iaunit
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: mask
#if defined(REAL4)
      real*4,  intent(out)   :: hmin,hmax
      real*4,  dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(out)   :: h
#else
      real,    intent(out)   :: hmin,hmax
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(out)   :: h
#endif
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for array reading.
!
!     must call zaiopn for this array unit before calling zaiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
!
!  4) hmin,hmax are returned as the minimum and maximum value in the 
!     array, ignoring array elements set to 2.0**100.  
!     if lmask==.true. the range is calculated only where mask.ne.0,
!     with all other values unchanged in h on exit.  It is then an
!     error if mask.ne.0 anywhere the input is 2.0**100.
!*
!**********
!
! --- spval  = data void marker, 2^100 or about 1.2676506e30
      real*4     spval
      parameter (spval=2.0**100)
!
      character cfile*256
      integer   ios, i,j
#if defined(TIMER)
!
      call xctmr0(17)
#endif
!
!     write(lp,*) 'zaiord - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiord)'
      endif
!
      if     (iarec(iaunit).eq.-99) then  !dummy I/O
        w(1:n2drec) = 0.0
      else !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
        call zaiordd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          write(lp,*) 'ios = ',ios
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          stop '(zaiord)'
        endif !ios
      endif
      if     (lmask) then
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1,jtdm
          wminy(j) =  spval  !simplifies OpenMP parallelization
          wmaxy(j) = -spval  !simplifies OpenMP parallelization
          do i= 1,itdm
            if     (mask(i,j).ne.0) then
              h(i,j)   =                w(i+(j-1)*itdm)
              wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
              wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
            endif
          enddo
        enddo
      else
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1,jtdm
          wminy(j) =  spval  !simplifies OpenMP parallelization
          wmaxy(j) = -spval  !simplifies OpenMP parallelization
          do i= 1,itdm
            h(i,j) = w(i+(j-1)*itdm)
            if     (w(i+(j-1)*itdm).ne.spval) then
              wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
              wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
            endif
          enddo
        enddo
      endif
      hmin = minval(wminy(1:jtdm))
      hmax = maxval(wmaxy(1:jtdm))
!
      if     (lmask .and. hmax.eq.spval) then
        write(lp,9200) iarec(iaunit),iaunit
        call flush(lp)
        cfile = ' '
        inquire(unit=iaunit+uaoff,name=cfile)
        write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
        call flush(lp)
        stop '(zaiord)'
      endif
!
#if defined(TIMER)
!
      call xctmr1(17)
#endif
      return
!
 9000 format(/ /10x,'error in zaiord -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiord -  can''t read record', &
         i4,' on array I/O unit ',i3,'.'/ /)
 9200 format(/ /10x,'error in zaiord -  record', &
         i4,' on array I/O unit ',i3, &
         ' has 2.0**100 outside masked region.'/ /)
      end subroutine zaiord

      subroutine zaiordg(hh4, iaunit)
      implicit none
!
      integer, intent(in)    :: iaunit
      real*4,  dimension (itdm,jtdm), &
               intent(out)   :: hh4
!
      integer        iarec
      common/czioxx/ iarec(999)
      save  /czioxx/
!
!**********
!*
!  1) machine specific routine for array reading.
!     returns non-tiled real*4 array.
!
!     must call zaiopn for this array unit before calling zaiordg.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
!*
!**********
!
      character cfile*256
      integer   ios, i,j
#if defined(TIMER)
!
      call xctmr0(17)
#endif
!
!     write(lp,*) 'zaiordg - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiordg)'
      endif
!
      if     (iarec(iaunit).eq.-99) then  !dummy I/O
        w(1:n2drec) = 0.0
      else !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
        call zaiordd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          write(lp,*) 'ios = ',ios
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          stop '(zaiord)'
        endif !ios
      endif
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j= 1,jtdm
        do i= 1,itdm
          hh4(i,j) = w(i+(j-1)*itdm)
        enddo
      enddo
!
#if defined(TIMER)
!
      call xctmr1(17)
#endif
      return
!
 9000 format(/ /10x,'error in zaiordg -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiordg -  can''t read record', &
         i4,' on array I/O unit ',i3,'.'/ /)
      end subroutine zaiordg

      subroutine zaiordd(a,n, iunit,irec,ios)
      implicit none
!
      integer, intent(in)    :: n,iunit,irec
      integer, intent(out)   :: ios
      real*4,  intent(out)   :: a(n)
!
!**********
!*
! 1)  direct access read a single record.
!
! 2)  expressed as a subroutine because i/o with 
!     implied do loops can be slow on some machines.
!*
!**********
!
      read(unit=iunit, rec=irec, iostat=ios) a
#if defined(ENDIAN_IO)
      call zaio_endian(a,n)
#endif
      return
      end subroutine zaiordd

      subroutine zaiosk(iaunit)
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
!  1) machine specific routine for skipping an array read.
!
!     must call zaiopn for this array unit before calling zaiosk.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
!*
!**********
#if defined(TIMER)
!
      call xctmr0(16)
#endif
!
!     write(lp,*) 'zaiosk - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call flush(lp)
        stop '(zaiosk)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
#if defined(TIMER)
!
      call xctmr1(16)
#endif
      return
!
 9000 format(/ /10x,'error in zaiosk -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
      end subroutine zaiosk

      subroutine zaiowr3(h, l, mask,lmask, hmin,hmax, iaunit, lreal4)
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
!     must call zaiopn for this array unit before calling zaiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
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
!     this version just calls zaiowr l times.
!
      integer k
!
      do k= 1,l
        call zaiowr(h(1-nbdy,1-nbdy,k), mask,lmask, &
                    hmin(k),hmax(k), iaunit, lreal4)
      enddo
      return
      end subroutine zaiowr3

      subroutine zaiowr(h, mask,lmask, hmin,hmax,  iaunit, lreal4)
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
!     must call zaiopn for this array unit before calling zaiord.
!
!  2) array i/o is fortran real*4 direct access i/o to unit iaunit+uaoff.
!
!  3) iaunit+uaoff is the i/o unit used for arrays.  array i/o might not
!      use fortran i/o units, but, for compatability, assume that
!      iaunit+uaoff refers to a fortran i/o unit anyway.
!     the array, 'h',  must conform to that passed in the associated
!      call to zaiopn.
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
        stop '(zaiowr)'
      endif
!
      if     (lreal4) then
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,itdm
              if     (mask(i,j).ne.0) then
                w(i+(j-1)*itdm) = h(i,j)
                wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
              else
                w(i+(j-1)*itdm) = spval
              endif
#if defined(REAL4)
! ---         h(i,j) = w(i+(j-1)*itdm)  ! h is already real*4
#else
              h(i,j) = w(i+(j-1)*itdm)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,itdm
              w(i+(j-1)*itdm) = h(i,j)
              if     (w(i+(j-1)*itdm).ne.spval) then
                wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
              endif
#if defined(REAL4)
! ---         h(i,j) = w(i+(j-1)*itdm)  ! h is already real*4
#else
              h(i,j) = w(i+(j-1)*itdm)  ! h is not real*4, so update it
#endif
            enddo
          enddo
        endif
      else
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,itdm
              if     (mask(i,j).ne.0) then
                w(i+(j-1)*itdm) = h(i,j)
                wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
              else
                w(i+(j-1)*itdm) = spval
              endif
            enddo
          enddo
        else
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jtdm
            wminy(j) =  spval  !simplifies OpenMP parallelization
            wmaxy(j) = -spval  !simplifies OpenMP parallelization
            do i= 1,itdm
              w(i+(j-1)*itdm) = h(i,j)
              if     (w(i+(j-1)*itdm).ne.spval) then
                wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
                wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
              endif
            enddo
          enddo
        endif
      endif
      do i= itdm*jtdm+1,n2drec
        w(i) = spval
      enddo
      hmin = minval(wminy(1:jtdm))
      hmax = maxval(wmaxy(1:jtdm))
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
        call zaiowrd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          stop '(zaiowr)'
        endif !ios
      endif
#if defined(TIMER)
!
      call xctmr1(18)
#endif
      return
!
 9000 format(/ /10x,'error in zaiowr -  array I/O unit ', &
         i3,' is not marked as open.'/ /)
 9100 format(/ /10x,'error in zaiowr -  can''t write record', &
         i4,' on array I/O unit ',i3,'.'/ /)
      end subroutine zaiowr
      subroutine zaiowrd(a,n, iunit,irec,ios)
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
      end subroutine zaiowrd
!
!
!> Revision history:
!>
!> July 2017 - added zaiordg
