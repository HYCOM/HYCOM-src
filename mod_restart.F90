      module mod_restart
      use mod_xc         ! HYCOM communication interface
      use mod_za         ! HYCOM I/O interface

      implicit none
     
!
! --- module for restart and related routines
!
      private !! default is private
      public  :: restart_in, restart_out, restart_zero

      interface restart_in
        module procedure restart_in_standalone
        module procedure restart_in_coupled
      end interface restart_in 

      interface restart_out
        module procedure restart_out_standalone
        module procedure restart_out_coupled
      end interface restart_out

      contains

      subroutine restart_in_standalone(nstep0, dtime0, flnmra,flnmrb)

      implicit none
      integer       nstep0
      real*8        dtime0
      character*(*) flnmra,flnmrb
!
      call restart_in_coupled(nstep0, dtime0, flnmra,flnmrb,.false.)

      end subroutine restart_in_standalone

      subroutine restart_out_standalone(nstepx, dtimex, flnmra,flnmrb, last)
      implicit none
!
      logical last
      integer nstepx
      real*8  dtimex
      character*(*) flnmra,flnmrb
!
      call restart_out_coupled(nstepx, dtimex, flnmra,flnmrb, last, .false.)

      end subroutine restart_out_standalone

      subroutine restart_in_coupled(nstep0, dtime0, flnmra,flnmrb, restart_cpl)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_tides      ! HYCOM tides
      implicit none
!
      integer       nstep0
      real*8        dtime0
      character*(*) flnmra,flnmrb
      logical, intent(in) :: restart_cpl ! coupled case
!
!     read in a restart file.
!     flnmra is the ".a" file, and flnmrb is the ".b" file.
!
      logical   lmyin,ltidin,lold
      integer   i,ios,j,k,kskip,ktr
      character cline*80
!
# include "stmt_fns.h"
!
      call zaiopf(flnmra,'old', 11)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+11,file=flnmrb, &
              status='old',action='read',form='formatted')
      endif
      call zagetc(cline,ios, uoff+11)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'I/O error from zagetc, iunit,ios = ',uoff+11,ios
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      if     (mnproc.eq.1) then
      write(lp,'(a)') trim(cline)
      endif !1st tile
      if     (cline(1:9).eq.'RESTART: ') then
        lold = .true.  !original, larger, restart file
      elseif (cline(1:9).eq.'RESTART2:') then
        lold = .false.
      else
        if     (mnproc.eq.1) then
        write(lp,'(/ a /)') 'error in hycom - unknown restart type'
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      call zagetc(cline,ios, uoff+11)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'I/O error from zagetc, iunit,ios = ',uoff+11,ios
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      if     (mnproc.eq.1) then
      write(lp,'(a)') trim(cline)
      call flush(lp)
      endif !1st tile
      i = index(cline,'=')
      read(cline(i+1:),*) nstep0,dtime0
!
      call restart_in3d(u,     2*kdm, iu, 'u       ')
      call restart_in3d(v,     2*kdm, iv, 'v       ')
      call restart_in3d(dp,    2*kdm, ip, 'dp      ')
      call restart_in3d(temp,  2*kdm, ip, 'temp    ')
      call restart_in3d(saln,  2*kdm, ip, 'saln    ')
      if     (lold) then
        call restart_in3d(th3d,  2*kdm, ip, 'th3d    ')
      else
        do k= 1,kdm
          do j= 1,jj
            do i= 1,ii
              if     (ip(i,j).eq.1) then
                th3d(i,j,k,1)=sig(temp(i,j,k,1),saln(i,j,k,1))-thbase
                th3d(i,j,k,2)=sig(temp(i,j,k,2),saln(i,j,k,2))-thbase
              else
                th3d(i,j,k,1) = 0.0
                th3d(i,j,k,2) = 0.0
              endif
            enddo !i
          enddo !j
        enddo !k
      endif
!
!     do we have MY2.5 arrays in the restart file?
!
      call zagetc(cline,ios, uoff+11)
      if     (ios.ne.0) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4,i9 /)') &
            'I/O error from zagetc, iunit,ios = ',uoff+11,ios
        endif !1st tile
        call xcstop('(restart_in)')
               stop '(restart_in)'
      endif
      if     (lold) then
        kskip = 12*kdm + 2
      else
        kskip = 10*kdm + 2
      endif
      call restart_inrw(kskip)
!
      lmyin = cline(1:8).eq.'q2      '
!
      if     (lmyin) then
!
!       MY2.5 in restart file.
!
        if (mxlmy) then
          call restart_in3d(q2    ,2*kdm+4, ip, 'q2      ')
          call restart_in3d(q2l   ,2*kdm+4, ip, 'q2l     ')
          call restart_in3d(vctymy,  kdm+2, ip, 'vctymy  ')
          call restart_in3d(difqmy,  kdm+2, ip, 'difqmy  ')
          call restart_in3d(diftmy,  kdm+2, ip, 'diftmy  ')
        else
          if     (mnproc.eq.1) then
          write(lp,'(a)') 'RESTART: skipping MY2.5 input fields'
          call flush(lp)
          endif !1st tile
          kskip = 7*kdm+14
          do k= 1,kskip
            call zagetc(cline,ios, uoff+11)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)') &
                  'I/O error from zagetc, iunit,ios = ',uoff+11,ios
              endif !1st tile
              call xcstop('(restart_in)')
                     stop '(restart_in)'
            endif
!           if     (mnproc.eq.1) then
!           write(lp,'(a)') cline
!           endif !1st tile
            call zaiosk(11)
          enddo !k
        endif !mxlmy:else
!
!       do we have DETIDE arrays in the restart file?
!
        call zagetc(cline,ios, uoff+11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'I/O error from zagetc, iunit,ios = ',uoff+11,ios
          endif !1st tile
          call xcstop('(restart_in)')
                 stop '(restart_in)'
        endif
        if     (lold) then
          kskip = 12*kdm + 2 + 7*kdm+14
        else
          kskip = 10*kdm + 2 + 7*kdm+14
        endif
        call restart_inrw(kskip)
      elseif (mxlmy) then
        if     (mnproc.eq.1) then
        write(lp,'(a)') 'RESTART: no MY2.5 fields input'
        call flush(lp)
        endif !1st tile
      endif !lmyin:mxlmy
!
!     do we have DETIDE arrays in the restart file?
!
      ltidin = cline(1:8).eq.'uhrly   '
!
      if     (ltidin) then
!
!       DETIDE in restart file.
!       25 or 49 hrs present?
!
        call restart_inrw(kskip+25)
        call zagetc(cline,ios, uoff+11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'I/O error from zagetc, iunit,ios = ',uoff+11,ios
          endif !1st tile
          call xcstop('(restart_in)')
                 stop '(restart_in)'
        endif
        if     (cline(1:8).eq.'uhrly   ') then !49 hrs
          nhrly = 49  ![uv]ntide will be initialised in tides_set (tides_detide)
        else
          nhrly = 25  ![uv]ntide will be initialised in tides_set (tides_detide)
        endif
        call restart_inrw(kskip)
!
        if (tidflg.gt.0) then
          if     (.not.allocated(uhrly)) then
             allocate(  uhrly(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,49), &
                        vhrly(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,49), &
                       untide(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                       vntide(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)    )
            call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy)*50 )
          else
            if     (mnproc.eq.1) then
              write(lp,'(/ a /)') &
                'error - uhrly already allocated'
            endif !1st tile
            call xcstop('(restart_in)')
                   stop '(restart_in)'
          endif !allocated:else
          call restart_in3d(uhrly ,nhrly, iu, 'uhrly   ')
          call restart_in3d(vhrly ,nhrly, iv, 'vhrly   ')
        else
          if     (mnproc.eq.1) then
          write(lp,'(a)') 'RESTART: skipping DETIDE input fields'
          call flush(lp)
          endif !1st tile
          do k= 1,2*nhrly
            call zagetc(cline,ios, uoff+11)
            if     (ios.ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(/ a,i4,i9 /)') &
                  'I/O error from zagetc, iunit,ios = ',uoff+11,ios
              endif !1st tile
              call xcstop('(restart_in)')
                     stop '(restart_in)'
            endif
!           if     (mnproc.eq.1) then
!           write(lp,'(a)') cline
!           endif !1st tile
            call zaiosk(11)
          enddo !k
        endif !tidflg:else
      elseif (tidflg.gt.0) then
! ---   [uv]hrly & [uv]ntide will be allocated/initialized in tides_set
        if     (mnproc.eq.1) then
        write(lp,'(a)') 'RESTART: no DETIDE fields input'
        call flush(lp)
        endif !1st tile
      endif !ltidin:tidflg
!
      call restart_in3d(ubavg,     3, iu, 'ubavg   ')
      call restart_in3d(vbavg,     3, iv, 'vbavg   ')
      call restart_in3d(pbavg,     3, ip, 'pbavg   ')
      call restart_in3d(pbot,      1, ip, 'pbot    ')
      call restart_in3d(psikk,kapnum, ip, 'psikk   ')  !kapnum 1 or 2
      call restart_in3d(thkk, kapnum, ip, 'thkk    ')  !kapnum 1 or 2
      call restart_in3d(dpmixl,    2, ip, 'dpmixl  ')
!
      if (sshflg.eq.2) then
! ---   montg_c makes montg1 steric SSH, assumes long term mean SSH is all steric
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              pbavg(i,j,1) = pbavg(i,j,1) - montg_c(i,j)*rhoref
              pbavg(i,j,2) = pbavg(i,j,2) - montg_c(i,j)*rhoref
              pbavg(i,j,3) = pbavg(i,j,3) - montg_c(i,j)*rhoref
            endif ! ip
          enddo !i
        enddo !j
      endif ! sshflg==2
!
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
             oneta(i,j,1) = max( oneta0, 1.0 + pbavg(i,j,1)/pbot(i,j) )
             oneta(i,j,2) = max( oneta0, 1.0 + pbavg(i,j,2)/pbot(i,j) )
          else
             oneta(i,j,1) = 1.0
             oneta(i,j,2) = 1.0
          endif !ip
        enddo !i
      enddo !j
      vland = 1.0
      call xctilr(oneta,1,2, nbdy,nbdy, halo_ps)
      vland = 0.0
!
      if (icegln) then
        call zagetc(cline,ios, uoff+11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'I/O error from zagetc, iunit,ios = ',uoff+11,ios
          endif !1st tile
          call xcstop('(restart_in)')
                 stop '(restart_in)'
        endif
        if     (ios.ne.0 .or. cline(1:8).ne.'temice  ') then
!
! ---     assume this is an addition of ice to the simulation.
!
          if     (mnproc.eq.1) then
          write(lp,'(/ a /)') 'adding ice to the simulation.'
          call flush(lp)
          endif !1st tile
!
          do j= 1,jj
            do i= 1,ii
              temice(i,j) = temp(i,j,1,1)
              covice(i,j) = 0.0
              thkice(i,j) = 0.0
            enddo
          enddo
          if     (trcrin .and. cline(1:8).eq.'tracer  ') then
! ---       reposition file for tracer input
            if     (lold) then
              kskip = 12*kdm+14+2*kapnum
            else
              kskip = 10*kdm+14+2*kapnum
            endif
            if     (lmyin) then
              kskip = kskip + 7*kdm+14
            endif
            if     (ltidin) then
              kskip = kskip + 2*nhrly
            endif
            call restart_inrw(kskip)
          endif
        else
!
! ---     reposition file for ice input
!
          if     (lold) then
            kskip = 12*kdm+14+2*kapnum
          else
            kskip = 10*kdm+14+2*kapnum
          endif
          if     (lmyin) then
            kskip = kskip + 7*kdm+14
          endif
          if     (ltidin) then
            kskip = kskip + 2*nhrly
          endif
          call restart_inrw(kskip)
!
          call restart_in3d(temice,    1, ip, 'temice  ')
          call restart_in3d(covice,    1, ip, 'covice  ')
          call restart_in3d(thkice,    1, ip, 'thkice  ')
        endif  !new ice:read ice
      else
! ---   no sea ice, but still need covice
        do j= 1,jj
          do i= 1,ii
            covice(i,j) = 0.0
          enddo
        enddo
      endif  !icegln:else
      if (trcrin) then
        do ktr= 1,ntracr
          call restart_in3d(tracer(1-nbdy,1-nbdy,1,1,ktr), &
                                   2*kdm, ip, 'tracer  ')
        enddo
      endif

#if defined (USE_NUOPC_CESMBETA)
      if (restart_cpl) then
        call zagetc(cline,ios, uoff+11)
        if (cline(1:8).eq. 'tml     ') then
!!Alex average export fields
! ---  reposition file for coupled input
!
          if     (lold) then
            kskip = 12*kdm+14+2*kapnum
          else
            kskip = 10*kdm+14+2*kapnum
          endif
          if     (lmyin) then
            kskip = kskip + 7*kdm+14
          endif
          if     (ltidin) then
            kskip = kskip + 2*nhrly
          endif
          call restart_inrw(kskip)
          call restart_in3d(tml ,    1, ip, 'tml     ')
          call restart_in3d(sml ,    1, ip, 'sml     ')
          call restart_in3d(umxl,    1, ip, 'umxl    ')
          call restart_in3d(vmxl,    1, ip, 'vmxl    ')
        else
! ---   reposition file for coupled input
!
          if     (lold) then
            kskip = 12*kdm+14+2*kapnum
          else
            kskip = 10*kdm+14+2*kapnum
          endif
          if     (lmyin) then
            kskip = kskip + 7*kdm+14
          endif
          if     (ltidin) then
            kskip = kskip + 2*nhrly
          endif
          call restart_inrw(kskip)
        endif !cline
      endif ! restart_cpl = .true.
#endif  /* USE_NUOPC_CESMBETA */

      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        close (unit=uoff+11)
      endif
      call zaiocl(11)
!
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          srfhgt(i,j) = 0.0 !for pipe_compareall
          montg1(i,j) = 0.0 !for pipe_compareall
            dpbl(i,j) = 0.0 !for pipe_compareall
           klist(i,j) = kk  !for MY2.5 mixed layer
        enddo
      enddo
      return
      end subroutine restart_in_coupled

      subroutine restart_in3d(field,l, mask, cfield)
!
      integer   l
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l) :: &
       field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       mask
      character cfield*8
!
! --- read a single restart 3-d array field.
!
      integer   i,ios,layer,level,k
      real      hmina(2*kdm+49),hminb,hmaxa(2*kdm+49),hmaxb  !+49 for [uv]hrly
      character cline*80
!
      if     (mnproc.eq.1) then
      write(lp,'(a,i3,2x,a)') 'restart_in3d - l,cfield = ',l,cfield
      call flush(lp)
      endif !1st tile
      call zaiord3(field,l, mask,.false., hmina,hmaxa, 11)
!
      do k= 1,l
        call zagetc(cline,ios, uoff+11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'I/O error from zagetc, iunit,ios = ',uoff+11,ios
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
        if     (mnproc.eq.1) then
        write (lp,'(a)')  trim(cline)
        endif !1st tile
        if     (cline(1:8).ne.cfield) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') trim(cline), &
                 'error in restart_in3d - expected ',cfield
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
        i = index(cline,'=')
        read (cline(i+1:),*) layer,level,hminb,hmaxb
        if     (abs(hmina(k)-hminb).gt.abs(hminb)*1.e-4 .or. &
                abs(hmaxa(k)-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,3i3 / a / a,1p3e14.6 / a,1p3e14.6 /)') &
            'error - .a and .b files not consistent:', &
            'iunit,k,l = ',11,k,l, &
            cline, &
            '.a,.b min = ',hmina(k),hminb,hmina(k)-hminb, &
            '.a,.b max = ',hmaxa(k),hmaxb,hmaxa(k)-hmaxb
          endif !1st tile
          call xcstop('(restart_in3d)')
                 stop '(restart_in3d)'
        endif
      enddo
!
      return
      end subroutine restart_in3d

      subroutine restart_inrw(kline)
!
      integer   kline
!
!     reposition the input restart .b file at line kline.
!
      integer   ios,k
      character cline*80
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        rewind(uoff+11)
      endif
      do k= 1,kline
        call zagetc(cline,ios, uoff+11)
        if     (ios.ne.0) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,i4,i9 /)') &
              'I/O error from zagetc, iunit,ios = ',uoff+11,ios
          endif !1st tile
          call xcstop('(restart_inrw)')
                 stop '(restart_inrw)'
        endif
      enddo !k
      if     (mnproc.eq.1) then
      write(lp,'(a,i5)') 'restart_inrw, kline =',kline
      write(lp,'(a)')    trim(cline)
      call flush(lp)
      endif !1st tile
      return
      end subroutine restart_inrw

      subroutine restart_out_coupled(nstepx, dtimex, flnmra,flnmrb, last, restart_cpl)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_tides      ! HYCOM tides
      implicit none

      logical last
      integer nstepx
      real*8  dtimex
      character*(*) flnmra,flnmrb
      logical, intent(in) :: restart_cpl ! coupled

!
!     write out in a restart file on unit 12 or 22 (and a flux file on 25).
!
!     flnmra is the ".a" file (usually without the .a, and
!     flnmrb is the ".b" file (usually without the .b).
!     Usually flnmra == flnmrb, and there are standard and backup restarts;
!     otherwise the restart is unique and flnmra and flnmrb are the complete
!     filenames (including any .a and .b).
!
      logical   lopen
      integer   i,iunit,iunta,j,k,ktr,l
      real      xmin(2*kdm+49),xmax(2*kdm+49)  !+49 for [uv]hrly
      character cline*80
!
      integer, save :: icount = 0
!
# include "stmt_fns.h"
!
      icount = icount + 1
!
      if     (flnmra.ne.flnmrb .or. last .or. mod(icount,2).eq.0) then
        iunta = 12  ! standard restart file
      else
        iunta = 22  ! backup   restart file
      endif
      iunit = uoff+iunta
!!Alex add unit for CESM restart
      if (restart_cpl) iunta = 15
!
      call zaiopi(lopen, iunta)
      if (.not.lopen) then
        if     (flnmra.ne.flnmrb) then
          call zaiopf(trim(flnmra),       'new', iunta)  !unique
        elseif (iunta.eq.12) then
          call zaiopf(trim(flnmra)//'.a', 'new', iunta)  !standard
        elseif (iunta.eq.15) then
          call zaiopf(trim(flnmra)//'.a', 'new', iunta)  !standard coupled !!Alex
        else
          call zaiopf(trim(flnmra)//'1.a','new', iunta)  !backup
        endif
        if     (mnproc.eq.1) then
          if     (flnmra.ne.flnmrb) then
            open (unit=iunit,file=trim(flnmrb),   & !12
                  status='new',action='write',form='formatted')
            write(lp,'(a)') ' creating a new unique   restart file'
          elseif (iunta.eq.12) then
            open (unit=iunit,file=trim(flnmra)//'.b',    & !12
                  status='new',action='write',form='formatted')
            write(lp,'(a)') ' creating a new standard restart file'
          elseif (iunta.eq.15) then
            open (unit=iunit,file=trim(flnmra)//'.b',    &!15 !!Alex
                  status='new',action='write',form='formatted')
            write(lp,'(a)') ' creating a new standard restart file'
          else
            open (unit=iunit,file=trim(flnmra)//'1.b',   & !22
                  status='new',action='write',form='formatted')
            write(lp,'(a)') ' creating a new backup   restart file'
          endif
          call flush(lp)
        endif !1st tile
      elseif (flnmra.ne.flnmrb) then
        if     (mnproc.eq.1) then
          write(lp,'(a)')  &
            ' error - (unique) restart file already exists.'
          write(lp,'(a,a)')  &
            ' flnmra = ',trim(flnmra)
          write(lp,'(a,a)')  &
            ' flnmrb = ',trim(flnmrb)
        endif !1st tile
        call xcstop('(restart_out)')
               stop '(restart_out)'
      else
        call zaiorw(iunta)
        if     (mnproc.eq.1) then
          rewind(unit=iunit)
          if     (iunta.eq.12) then
            write(lp,'(a)')  &
              ' over-writing any previous standard restart'
          else
            write(lp,'(a)')  &
              ' over-writing any previous backup   restart'
          endif
          call flush(lp)
        endif !1st tile
      endif
!
      if     (mnproc.eq.1) then
      write(iunit,'(a,4i6)') 'RESTART2: iexpt,iversn,yrflag,sigver = ', &
                                    iexpt,iversn,yrflag,sigver
      write(cline,*)                nstepx,dtimex,thbase
      write(iunit,'(a,a)')   'RESTART2: nstep,dtime,thbase = ', &
                                    trim(cline)
      call flush(iunit)
      endif !1st tile
!
      call zaiowr3(u,      2*kdm, iu,.false., xmin,xmax, iunta,.true.)
      call xctilr( u,    1,2*kdm, nbdy,nbdy, halo_uv)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(iunit,4100) 'u       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(v,      2*kdm, iv,.false., xmin,xmax, iunta,.true.)
      call xctilr( v,    1,2*kdm, nbdy,nbdy, halo_vv)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(iunit,4100) 'v       ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(dp,     2*kdm, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( dp,   1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(iunit,4100) 'dp      ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(temp,   2*kdm, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( temp, 1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(iunit,4100) 'temp    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(saln,   2*kdm, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( saln, 1,2*kdm, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 0,1
        do k= 1,kdm
          write(iunit,4100) 'saln    ',k,l+1,xmin(k+l*kdm),xmax(k+l*kdm)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile

!
! --- temp and saln may have been changed, so update th3d
      do k= 1,kdm
        do j= 1-nbdy,jj+nbdy
          do i= 1-nbdy,ii+nbdy
            if     (ip(i,j).eq.1) then
              th3d(i,j,k,1)=sig(temp(i,j,k,1),saln(i,j,k,1))-thbase
              th3d(i,j,k,2)=sig(temp(i,j,k,2),saln(i,j,k,2))-thbase
            else
              th3d(i,j,k,1) = 0.0
              th3d(i,j,k,2) = 0.0
            endif
          enddo !i
        enddo !j
      enddo !k
!
      if (mxlmy) then
        call zaiowr3(q2, 2*kdm+4, ip,.false., xmin,xmax, iunta,.true.)
        if     (mnproc.eq.1) then
        do l= 0,1
          do k= 1,kdm+2
            write(iunit,4100) 'q2      ' &
                          ,k,l+1,xmin(k+l*(kdm+2)),xmax(k+l*(kdm+2))
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(q2l, 2*kdm+4, ip,.false., xmin,xmax, iunta,.true.)
        if     (mnproc.eq.1) then
        do l= 0,1
          do k= 1,kdm+2
            write(iunit,4100) 'q2l     ' &
                          ,k,l+1,xmin(k+l*(kdm+2)),xmax(k+l*(kdm+2))
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(vctymy, kdm+2, ip,.false., xmin,xmax, iunta,.true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(iunit,4100) 'vctymy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(difqmy, kdm+2, ip,.false., xmin,xmax, iunta,.true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(iunit,4100) 'difqmy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(diftmy, kdm+2, ip,.false., xmin,xmax, iunta,.true.)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 1,kdm+2
            write(iunit,4100) 'diftmy  ',k,l,xmin(k),xmax(k)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
      endif !mxlmy
!
      if (tidflg.gt.0) then
        call zaiowr3(uhrly,   49, iu,.false., xmin,xmax, iunta,.true.)
        call xctilr( uhrly, 1,49, nbdy,nbdy, halo_uv)
        if     (mnproc.eq.1) then
        do l= 1,49
          do k= 0,0
            write(iunit,4100) 'uhrly   ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(vhrly,   49, iv,.false., xmin,xmax, iunta,.true.)
        call xctilr( vhrly, 1,49, nbdy,nbdy, halo_vv)
        if     (mnproc.eq.1) then
        do l= 1,49
          do k= 0,0
            write(iunit,4100) 'vhrly   ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
      endif !tidflg
!
      call zaiowr3(ubavg,      3, iu,.false., xmin,xmax, iunta,.true.)
      call xctilr( ubavg,    1,3, nbdy,nbdy, halo_uv)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(iunit,4100) 'ubavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(vbavg,      3, iv,.false., xmin,xmax, iunta,.true.)
      call xctilr( vbavg,    1,3, nbdy,nbdy, halo_vv)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(iunit,4100) 'vbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      if (sshflg.eq.2) then
! ---   unwind the pbavg correction for compatibility with psikk
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              pbavg(i,j,1) = pbavg(i,j,1) + montg_c(i,j)*rhoref
              pbavg(i,j,2) = pbavg(i,j,2) + montg_c(i,j)*rhoref
              pbavg(i,j,3) = pbavg(i,j,3) + montg_c(i,j)*rhoref
            endif ! ip
          enddo !i
        enddo !j
      endif ! sshflg==2
      call zaiowr3(pbavg,      3, ip,.false., xmin,xmax, iunta,.true.)
      if (sshflg.eq.2) then
! ---   re-apply the pbavg correction
        do j= 1,jj
          do i= 1,ii
            if     (ip(i,j).eq.1) then
              pbavg(i,j,1) = pbavg(i,j,1) - montg_c(i,j)*rhoref
              pbavg(i,j,2) = pbavg(i,j,2) - montg_c(i,j)*rhoref
              pbavg(i,j,3) = pbavg(i,j,3) - montg_c(i,j)*rhoref
            endif ! ip
          enddo !i
        enddo !j
      endif ! sshflg==2
      call xctilr( pbavg,    1,3, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,3
        do k= 0,0
          write(iunit,4100) 'pbavg   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(pbot,       1, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( pbot,     1,1, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,1
        do k= 0,0
          write(iunit,4100) 'pbot    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(psikk,  kapnum,ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( psikk,1,kapnum,nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,kapnum  !kapnum 1 or 2
        do k= 0,0
          write(iunit,4100) 'psikk   ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(thkk,  kapnum, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( thkk,1,kapnum, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,kapnum  !kapnum 1 or 2
        do k= 0,0
          write(iunit,4100) 'thkk    ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
      call zaiowr3(dpmixl,     2, ip,.false., xmin,xmax, iunta,.true.)
      call xctilr( dpmixl,   1,2, nbdy,nbdy, halo_ps)
      if     (mnproc.eq.1) then
      do l= 1,2
        do k= 0,0
          write(iunit,4100) 'dpmixl  ',k,l,  xmin(l),xmax(l)
        enddo
      enddo
      call flush(iunit)
      endif !1st tile
!
! --- needed because zaiowr fields have been cast to REAL*4
!
      do j= 1,jj
        do i= 1,ii
          if     (ip(i,j).eq.1) then
            oneta(i,j,1) = max( oneta0, 1.0 + pbavg(i,j,1)/pbot(i,j) )
            oneta(i,j,2) = max( oneta0, 1.0 + pbavg(i,j,2)/pbot(i,j) )
          else
            oneta(i,j,1) = 1.0
            oneta(i,j,2) = 1.0
          endif !ip
        enddo !i
      enddo !j
      vland = 1.0
      call xctilr(oneta,1,2, nbdy,nbdy, halo_ps)
      vland = 0.0

      if (icegln) then
        call zaiowr3(temice,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( temice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'temice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(covice,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( covice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'covice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(thkice,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( thkice,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'thkice  ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
      endif
      if (trcout) then
        do ktr= 1,ntracr
          call zaiowr3(tracer(1-nbdy,1-nbdy,1,1,ktr),  2*kdm, &
                       ip,.false., xmin,xmax, iunta,.true.)
          call xctilr( tracer(1-nbdy,1-nbdy,1,1,ktr),1,2*kdm, &
                        nbdy,nbdy, halo_ps)
          if     (mnproc.eq.1) then
          do l= 0,1
            do k= 1,kdm
              write(iunit,4100) 'tracer  ',k,l+1,xmin(k+l*kdm), &
                                              xmax(k+l*kdm)
            enddo
          enddo
          call flush(iunit)
          endif !1st tile
        enddo !ktr
      endif !trcout
#if defined (USE_NUOPC_CESMBETA)
      if (restart_cpl) then
!!Alex averaged export fields
        call zaiowr3(tml,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( tml,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'tml     ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(sml,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( sml,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'sml     ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(umxl,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( umxl,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'umxl    ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
        call zaiowr3(vmxl,     1, ip,.false., xmin,xmax, iunta,.true.)
        call xctilr( vmxl,   1,1, nbdy,nbdy, halo_ps)
        if     (mnproc.eq.1) then
        do l= 1,1
          do k= 0,0
            write(iunit,4100) 'vmxl    ',k,l,  xmin(l),xmax(l)
          enddo
        enddo
        call flush(iunit)
        endif !1st tile
      endif ! restart_cpl = .true.
#endif /* USE_NUOPC_CESMBETA */  
      if     (flnmra.ne.flnmrb) then  !unique restart file
        call zaiocl(iunta)
        if     (mnproc.eq.1) then
          close(unit=iunit)
          write(lp,'(a,f11.3)')  &
            ' unique restart created at model day',dtimex
          call flush(lp)
        endif !1st tile
      elseif (last) then !close all restart files
        call zaiocl(iunta) !iunta==12
        if     (mnproc.eq.1) then
          close(unit=iunit)
          write(lp,'(a,f11.3)')  &
            ' restart created & closed at model day',dtimex
          call flush(lp)
        endif
        call zaiopi(lopen, 22) !backup restart file?
        if (lopen) then
          call zaiocl(22)
          if     (mnproc.eq.1) then
            close(unit=uoff+22)
          endif
        endif
      else
        call zaiofl(iunta)
        if     (mnproc.eq.1) then
          call flush(iunit)
          write(lp,'(a,f11.3)')  &
            ' restart created at model day',dtimex
          call flush(lp)
        endif
      endif
      call xcsync(flush_lp)
!
! --- output to flux file
!
      if (.FALSE.) then ! turn on/off flux output
        call zaiopi(lopen, 25)
        if (.not.lopen) then
          call zaiopf(trim(flnmflx)//'.a','new', 25)
          if     (mnproc.eq.1) then
          open (unit=uoff+25,file=trim(flnmflx)//'.b', &
                status='new',action='write',form='formatted')
          write(uoff+25,'(a,3i6)') 'FLUXES: iexpt,iversn,yrflag = ', &
                                       iexpt,iversn,yrflag
          call flush(uoff+25)
          endif !1st tile
        endif
!
        if     (mnproc.eq.1) then
        write(cline,*)               nstepx,dtimex
        write(uoff+25,'(a,a)')   'FLUXES: nstep,dtime = ', &
                                     trim(cline)
        call flush(uoff+25)
        endif !1st tile
!
        call zaiowr3(dpav,     kdm, ip,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(uoff+25,4100) 'dpav    ',k,l, &
                                xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush(uoff+25)
        endif !1st tile
        call zaiowr3(uflxav,   kdm, iu,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(uoff+25,4100) 'uflxav  ',k,l, &
                                xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush(uoff+25)
        endif !1st tile
        call zaiowr3(vflxav,   kdm, iv,.true.,  xmin,xmax, 25, .false.)
        if     (mnproc.eq.1) then
        do l= 0,0
          do k= 1,kdm
            write(uoff+25,4100) 'vflxav  ',k,l, &
                                xmin(k+l*kdm),xmax(k+l*kdm)
          enddo
        enddo
        call flush( uoff+25)
        endif !1st tile
        call zaiofl(25)
      endif !flux output
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          srfhgt(i,j) = 0.0 !consistent with restart_in
          montg1(i,j) = 0.0 !consistent with restart_in
            dpbl(i,j) = 0.0 !consistent with restart_in
           klist(i,j) = kk  !for MY2.5 mixed layer
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      return
 4100 format(a,': layer,tlevel,range = ',i3,i3,2x,1p2e16.7)
      end subroutine restart_out_coupled

      subroutine restart_zero
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_tides      ! HYCOM tides
      implicit none

!
!     replacement for restart_in in dummy version
!     set all fields to zero
!
      call restart_zero3d(u,     2*kdm, iu, 'u       ')
      call restart_zero3d(v,     2*kdm, iv, 'v       ')
      call restart_zero3d(dp,    2*kdm, ip, 'dp      ')
      call restart_zero3d(temp,  2*kdm, ip, 'temp    ')
      call restart_zero3d(saln,  2*kdm, ip, 'saln    ')
      call restart_zero3d(th3d,  2*kdm, ip, 'th3d    ')
!
      call restart_zero3d(ubavg,     3, iu, 'ubavg   ')
      call restart_zero3d(vbavg,     3, iv, 'vbavg   ')
      call restart_zero3d(pbavg,     3, ip, 'pbavg   ')
      call restart_zero3d(pbot,      1, ip, 'pbot    ')
      call restart_zero3d(psikk,kapnum, ip, 'psikk   ')  !kapnum 1 or 2
      call restart_zero3d(thkk, kapnum, ip, 'thkk    ')  !kapnum 1 or 2
      call restart_zero3d(dpmixl,    2, ip, 'dpmixl  ')
      call restart_zero3d(temice,    1, ip, 'temice  ')
      call restart_zero3d(covice,    1, ip, 'covice  ')
      call restart_zero3d(thkice,    1, ip, 'thkice  ')
!
! --- not in a actual restart, but zero them anyway
!
      call restart_zero3d(srfhgt,    1, ip, 'srfhgt  ')
      call restart_zero3d(montg1,    1, ip, 'montg1  ')
      call restart_zero3d(dpbl,      1, ip, 'dpbl    ')
      return
      end subroutine restart_zero

      subroutine restart_zero3d(field,l, mask, cfield)
!
      integer   l
      real,    dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,l) :: &
       field
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       mask
      character cfield*8
!
! --- zero a single restart 3-d array field.
!
      if     (mnproc.eq.1) then
      write(lp,'(a,i3,2x,a)') 'restart_zero3d - l,cfield = ',l,cfield
      call flush(lp)
      endif !1st tile
      field(:,:,:) = 0.0
!
      return
      end subroutine restart_zero3d
!
      end module mod_restart

!
!> Revision history:
!>
!> May  2007 - removed th3d from the restart file
!> Mar. 2010 - removed  DETIDE from the restart file
!> Apr. 2010 - put back DETIDE into the restart file
!> Aug. 2010 - 49-hour DETIDE
!> May  2012 - added restart_zero
!> Sep. 2015 - if no sea ice, still set covice to zero
!> Aug. 2018 - added onetai
!> Dec. 2018 - module form 
!> Dec. 2018 - added onetai and oneta to allow for type casting to REAL*4
!> Feb. 2019 - onetai set to 1.0
!> Feb  2019 - montg_c correction to pbavg (see momtum for correction to psikk)
!> Feb. 2019 - replaced onetai by 1.0
