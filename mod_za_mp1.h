!
!-----------------------------------------------------------------------
!
!     machine dependent I/O routines.
!     message passing version, with all I/O from first processor.
!     contained in module mod_za.
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
!
!  2) The read is performed on the first processor only.
!*
!**********
!
      integer        iline,ibuf
      common/czgetc/ iline(81,0:1),ibuf
      save  /czgetc/
!
      integer i
!
! --- I/O from first processor only
!
      ibuf = mod(ibuf+1,2)
!
      if     (mnproc.eq.1) then
        read(iunit,'(a)',iostat=ios) cline
        do i= 1,80
          iline(i,ibuf) = ichar(cline(i:i))
        enddo
        iline(81,ibuf) = ios
      endif
!
!     broadcast to all other processors
!
      call xcgetc(iline(:,ibuf))
      do i= 1,80
        cline(i:i) = char(iline(i,ibuf))
      enddo
      ios = iline(81,ibuf)  ! iostat value
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
! --- n2drec = size of output 2-d array, multiple of 4096
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
        call xcstop('(zaiopn)')
               stop '(zaiopn)'
      endif
!
      iarec(iaunit) = 0
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiopn - iaunit = ',iaunit
!     call flush(lp)
!
!     get filename.
!
      write(cenv,"('FOR',i3.3,'A')") iaunit
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(cfile,"('fort.',i3.3,'a')") iaunit
      endif
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
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
        call xchalt('(zaiopn)')
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
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif !ios
#endif
      if     (cstat.eq.'scratch' .or. &
              cstat.eq.'SCRATCH'     ) then
        open(unit=iaunit+uaoff, &
             form='unformatted', status='scratch', &
             access='direct', recl=nrecl, action=cact, iostat=ios)
      else
        open(unit=iaunit+uaoff, file=cfile, &
             form='unformatted', status=cstat, &
             access='direct', recl=nrecl, action=cact, iostat=ios)
      endif
      if     (ios.ne.0) then
        write(lp,9100) iaunit
        write(lp,*) 'ios = ',ios
        call flush(lp)
        if     (cstat.eq.'scratch' .or. &
                cstat.eq.'SCRATCH'     ) then
          write(6,'(a)')  "status='SCRATCH'"
        else
          write(6,'(3a)') 'FILENAME="',trim(cfile),'"'
        endif
        call flush(lp)
        call xchalt('(zaiopn)')
               stop '(zaiopn)'
      endif !ios
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
! --- n2drec = size of output 2-d array, multiple of 4096
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
        call xcstop('(zaiope)')
               stop '(zaiope)'
      endif
!
      iarec(iaunit) = 0
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiope - iaunit = ',iaunit
!     call flush(lp)
!
!     get filename.
!
      cfile = ' '
      call getenv(cenv,cfile)
      if     (cfile.eq.' ') then
        write(lp,9300) trim(cenv)
        write(lp,*) 'iaunit = ',iaunit
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
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
        call xchalt('(zaiope)')
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
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile, &
           form='unformatted', status=cstat, &
           access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        write(lp,*) 'cenv = ',trim(cenv)
        call flush(lp)
        call xchalt('(zaiope)')
               stop '(zaiope)'
      endif !ios
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
        call xcstop('(zaiopf)')
               stop '(zaiopf)'
      endif
!
      iarec(iaunit) = 0
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiopf - iaunit = ',iaunit
!     call flush(lp)
!
!     dummy I/O?
!
      if     (cfile.eq.'/dev/null' .or. &
              cfile.eq.'/dev/zero'     ) then
        iarec(iaunit) = -99
        call xcsync(no_flush)
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
        call xchalt('(zaiopf)')
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
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif !ios
#endif
      open(unit=iaunit+uaoff, file=cfile, &
           form='unformatted', status=cstat, &
           access='direct', recl=nrecl, action=cact, iostat=ios)
      if     (ios.ne.0) then
        write(lp,9100) iaunit,trim(cfile)
        write(lp,*) 'ios  = ',ios
        call flush(lp)
        call xchalt('(zaiopf)')
               stop '(zaiopf)'
      endif !ios
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
      if     (mnproc.eq.1) then
        write(lp,'(/a/)') &
          'zaiost - Array I/O is Fortran DA I/O from the 1st task'
      endif
!
      do i= 1,999
        iarec(i) = -1
      enddo
!
! --- w on 1st tile only
!
#if defined(RELO)
!
! --- n2drec = size of output 2-d array, multiple of 4096
      n2drec = ((itdm*jtdm+4095)/4096)*4096
!
      if     (mnproc.eq.1) then
        allocate( w(n2drec),wminy(jtdm),wmaxy(jtdm),htmp(idm*jdm) )
        call mem_stat_add( (n2drec+2*jtdm+idm*jdm)/2) !real*4, so /2
      else
        allocate( w(1),     wminy(jtdm),wmaxy(jtdm),htmp(idm*jdm) )
        call mem_stat_add( (     1+2*jtdm+idm*jdm)/2) !real*4, so /2
      endif !1st tile:else
#else
      if     (mnproc.eq.1) then
        allocate( w(n2drec) )
      else
        allocate( w(1) )
      endif !1st tile:else
#endif
#if defined(TIMER)
!
!     initialize timers.
!
      call xctmrn(16,'zaio**')
      call xctmrn(17,'zaiord')
      call xctmrn(18,'zaiowr')
      call xctmrn(19,'zaioIO')
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
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiocl)')
               stop '(zaiocl)'
      endif
!
      iarec(iaunit) = -1
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiocl - iaunit = ',iaunit
!     call flush(lp)
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        close(unit=iaunit+uaoff, status='keep')
#if defined(T3E) || defined(YMP) || defined(X1)
        call asnunit(iaunit+uaoff,'-R',ios)
#endif
      endif
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
        call xcstop('(zaiofl)')
               stop '(zaiofl)'
      endif
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        inquire(unit=iaunit+uaoff, name=cfile, recl=irlen)
        close(  unit=iaunit+uaoff, status='keep')
        open(   unit=iaunit+uaoff, file=cfile, form='unformatted', &
                access='direct', recl=irlen)
      endif
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
        call xcstop('(zaiorw)')
               stop '(zaiorw)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = 0
      endif
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiorw - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
!
!  5) Optimized by Dan Moore, Planning Systems Inc., August 2005.
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
      real*4    wmin,wmax
      real      rmin(1),rmax(1)
#if defined(TIMER)
!
      call xctmr0(17)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiord)')
               stop '(zaiord)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
!
      wmin =  spval
      wmax = -spval
!
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*)'In initial proc 1 section'
!     write(lp,*) 'zaiord - iaunit,rec = ',iaunit,iarec(iaunit)
!     write(lp,*) 'zaiord - mask.1,1    = ',amsk(1,1)
!     write(lp,*) 'zaiord - h.1,1       = ',atmp(1,1)
!     call flush(lp)
!
      if     (iarec(iaunit).eq.-99) then  !dummy I/O
        w(1:n2drec) = 0.0
      else !standard I/O
        call zaiordd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          write(lp,*) 'ios = ',ios
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          call xchalt('(zaiord)')
                 stop '(zaiord)'
        endif !ios
       endif ! dummy I/O
!      write(lp,*)'I/O on proc 1 finished?'
!      write(lp,*)'lmask = ',lmask
!
       if     (.not.lmask)then
!        Get global min and max on first processor.
!        must be done here because tiles need not cover the full domain.
!$OMP    PARALLEL DO PRIVATE(j,i) &
!$OMP             SCHEDULE(STATIC,jblk)
         do j= 1,jtdm
#if defined(ENDIAN_IO)
           call zaio_endian(w(1+(j-1)*itdm),itdm)  !swap to big-endian
#endif
           wminy(j) =  spval
           wmaxy(j) = -spval
           do i= 1,itdm
             if     (w(i+(j-1)*itdm).ne.spval) then
               wminy(j) = min( wminy(j), w(i+(j-1)*itdm) )
               wmaxy(j) = max( wmaxy(j), w(i+(j-1)*itdm) )
             endif
           enddo !i
         enddo !j
         wmin = minval(wminy(1:jtdm))
         wmax = maxval(wmaxy(1:jtdm))
       endif !Not Lmask (global min,max on processor #1)
!
       endif !end I/O from first processor only
!
! --- put field from 1st processor to all tiles
      call xcaput4(w,htmp, 1)  !w cast to a 2-d array
!
! --- Each processor loads h from htmp (where mask = 1)
! --- Each processor does local min max if lmask is true
!
      if     (lmask) then
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1,jj
#if defined(ENDIAN_IO)
          call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          wminy(j) =  spval
          wmaxy(j) = -spval
          do i= 1,ii
            if     (mask(i,j).ne.0) then
              h(i,j) =                  htmp(i+(j-1)*ii)
              wminy(j) = min( wminy(j), htmp(i+(j-1)*ii) )
              wmaxy(j) = max( wmaxy(j), htmp(i+(j-1)*ii) )
            endif
          enddo !i
        enddo !j
        wmin = minval(wminy(1:jj))
        wmax = maxval(wmaxy(1:jj))
      else
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j= 1,jj
          do i= 1,ii
            h(i,j) = htmp(i+(j-1)*ii)
          enddo !i
        enddo !j
      endif !lmask:else
!
! --- Min/Max broadcast/gather
!
      rmin(1) = wmin
      rmax(1) = wmax
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
!
      if     (lmask .and. hmax.eq.spval) then
        if     (mnproc.eq.1) then
        write(lp,9200) iarec(iaunit),iaunit
        call flush(lp)
        cfile = ' '
        inquire(unit=iaunit+uaoff,name=cfile)
        write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
        call flush(lp)
        endif !master
        call xcstop('(zaiord)')
               stop '(zaiord)'
      endif
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
! --- n2drec = size of output 2-d array, multiple of 4096
!
      character cfile*256
      integer   ios, i,j
#if defined(TIMER)
!
      call xctmr0(17)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiordg)')
               stop '(zaiordg)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
!
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*)'In initial proc 1 section'
!     write(lp,*) 'zaiordg - iaunit,rec = ',iaunit,iarec(iaunit)
!     write(lp,*) 'zaiordg - h.1,1       = ',atmp(1,1)
!     call flush(lp)
!
      if     (iarec(iaunit).eq.-99) then  !dummy I/O
        w(1:n2drec) = 0.0
      else !standard I/O
        call zaiordd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
        if     (ios.ne.0) then
          write(lp,9100) iarec(iaunit),iaunit
          write(lp,*) 'ios = ',ios
          call flush(lp)
          cfile = ' '
          inquire(unit=iaunit+uaoff,name=cfile)
          write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
          call flush(lp)
          call xchalt('(zaiordg)')
                 stop '(zaiordg)'
        endif !ios
#if defined(ENDIAN_IO)
        call zaio_endian(w,itdm*jtdm)  !swap to big-endian
#endif

       endif ! dummy:standard I/O
!      write(lp,*)'I/O on proc 1 finished?'
!
       endif !end I/O from first processor only
!
! --- put field from 1st processor to all tiles
      call xcgput4(w,hh4, 1)  !w cast to a 2-d array
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
#if defined(ENDIAN_IO)
!     on exit, a is still little-endian
#endif
!
! 2)  expressed as a subroutine because i/o with
!     implied do loops can be slow on some machines.
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(19)
#endif
      read(unit=iunit, rec=irec, iostat=ios) a
#if defined(TIMER)
!
      call xctmr1(19)
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
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiosk)')
               stop '(zaiosk)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
!
! --- I/O from first processor only
!
      if     (mnproc.eq.1) then
!
!     write(lp,*) 'zaiosk - iaunit,rec = ',iaunit,iarec(iaunit)
!     call flush(lp)
!
      endif  ! I/O from first processor only
      call xcsync(no_flush)
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
!
!  6) Optimized by Dan Moore, Planning Systems Inc., August 2005.
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
      real      rmin(1),rmax(1)
      real*4    data_void(1),vsave4
#if defined(TIMER)
!
      call xctmr0(18)
#endif
!
      if     (iarec(iaunit).eq.-1) then
        write(lp,9000) iaunit
        call xcstop('(zaiowr)')
               stop '(zaiowr)'
      endif
!
      if     (iarec(iaunit).ne.-99) then  !standard I/O
        iarec(iaunit) = iarec(iaunit) + 1
      endif
!
      data_void(1) = spval
#if defined(ENDIAN_IO)
      call zaio_endian(data_void,1) !swap to big-endian
#endif
!
! --- Copy into real*4 buffer, and find Min,Max
!
      if     (lreal4) then
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                htmp(i+(j-1)*ii) = h(i,j)
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              else
                htmp(i+(j-1)*ii) = spval
              endif
#if defined(REAL4)
! ---         h(i,j) = htmp(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = htmp(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        else !.not.lmask
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              htmp(i+(j-1)*ii) = h(i,j)
              if     (htmp(i+(j-1)*ii).ne.spval) then
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              endif
#if defined(REAL4)
! ---         h(i,j) = htmp(i+(j-1)*ii)  ! h is already real*4
#else
              h(i,j) = htmp(i+(j-1)*ii)  ! h is not real*4, so update it
#endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        endif !lmask:else
      else !.not.lreal4
        if     (lmask) then
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              if     (mask(i,j).ne.0) then
                htmp(i+(j-1)*ii) = h(i,j)
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              else
                htmp(i+(j-1)*ii) = spval
              endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        else !.not.lmask
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j= 1,jj
            wminy(j)= spval
            wmaxy(j)=-spval
            do i= 1,ii
              htmp(i+(j-1)*ii) = h(i,j)
              if     (htmp(i+(j-1)*ii).ne.spval) then
                wminy(j)=min(wminy(j),htmp(i+(j-1)*ii))
                wmaxy(j)=max(wmaxy(j),htmp(i+(j-1)*ii))
              endif
            enddo !i
#if defined(ENDIAN_IO)
            call zaio_endian(htmp(1+(j-1)*ii),ii)  !swap to big-endian
#endif
          enddo !j
        endif !lmask:else
      endif !lreal4:else
!
      rmin(1) = minval(wminy(1:jj))
      rmax(1) = maxval(wmaxy(1:jj))
      call xcminr(rmin)
      call xcmaxr(rmax)
      hmin = rmin(1)
      hmax = rmax(1)
!
! --- I/O from first processor only
#if defined(ENDIAN_IO)
! --- htmp and data_void are already big-endian
#endif
!
      vsave4 = vland4
      vland4 = data_void(1)
      call xcaget4(w,htmp, 1)  !htmp to w (as a 2-d array) on 1st processor
      vland4 = vsave4
!
      if     (mnproc.eq.1) then
        do i= itdm*jtdm+1,n2drec
          w(i) = data_void(1)
        enddo
        if     (iarec(iaunit).ne.-99) then  !standard I/O
          call zaiowrd(w,n2drec, iaunit+uaoff,iarec(iaunit),ios)
          if     (ios.ne.0) then
            write(lp,9100) iarec(iaunit),iaunit
            call flush(lp)
            cfile = ' '
            inquire(unit=iaunit+uaoff,name=cfile)
            write(lp,'(3a)') 'FILENAME="',trim(cfile),'"'
            call flush(lp)
            call xchalt('(zaiowr)')
                   stop '(zaiowr)'
          endif !ios
        endif !standard I/O
      endif !I/O from first processor only
!
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
      real*4,  intent(in)    :: a(n)
!
!**********
!*
! 1)  direct access write a single record.
#if defined(ENDIAN_IO)
!     a is already big-endian.
#endif
!
! 2)  expressed as a subroutine because i/o with
!     implied do loops can be slow on some machines.
!*
!**********
!
#if defined(TIMER)
!
      call xctmr0(19)
#endif
      write(unit=iunit, rec=irec, iostat=ios) a
#if defined(TIMER)
!
      call xctmr1(19)
#endif
      return
      end subroutine zaiowrd
!
!
!> Revision history:
!>
!> Feb  2015 - w on 1st tile only
!> July 2017 - added zaiordg
