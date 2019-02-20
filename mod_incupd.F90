#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#endif
      module mod_incupd
      use mod_xc  ! HYCOM communication interface
!
      implicit none
!
! --- HYCOM incremental updating (for data assimilation).
! --- Stochastic anomaly forcing is also handled here for convenence.
!
      integer, save, public  :: &
       incflg,    & ! incremental update flag (0=no,-/+1=yes,-/+2=full-velocity)
                    !   -ve to read in difference archive for increment
       incstp,    & ! no. timesteps for full update (1=full insertion)
       incupf,    & ! number of days of incremental updating input
                    !   -ve to write a restart after each update completes
       incice       ! direct insertion of sea ice concentration flag
!
      logical, save, private ::  &
       lrestart  ! write a restart after each update completes
!
      integer, save, private ::  &
       ncount,    & ! increment time step counter
       ncountd      ! increment day counter
!
      real*8, save, private ::  &
       dtimeu    ! next days increment field
!
      real,    allocatable, dimension(:,:), &
               save, private :: &
        ubinc,    & !  ubaro increment
        vbinc,    & !  vbaro increment
       covinc,    & ! covice increment for direct insertion
       thkinc       ! thkice increment for direct insertion

!
      real,    allocatable, dimension(:,:,:), &
               save, private :: &
        tinc,     & !     t increment
        sinc,     & !     s increment
       dpinc,     & !    dp increment
        uinc,     & !     u increment
        vinc        !     v increment

      contains

      subroutine incupd_init(dtime0)
!
      real*8 dtime0
!
! --- subroutine used to calculate increment field for the incremental updating
! --- version: dec 2005
!
      integer i,j,k
      logical lopen
!
! --- allocate arrays
!
      allocate(  tinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                 sinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                dpinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                 uinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                 vinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm), &
                ubinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vbinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               covinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               thkinc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     )
      call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy)*kdm )
      call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) )
!
! --- set counter to zero
!
      ncount  = 0
      ncountd = 0
      dtimeu  = 1.d-6
!
! --- restart flag
!
      lrestart = incupf.lt.0
      incupf   = abs(incupf)
!
! --- read the target fields, and initialize the "inc" arrays.
!
      if     (incflg.gt.0) then
        call incupd_read_full(dtime0)
      else  !(incflg.le.0)
        call incupd_read_diff(dtime0)
      endif
!
      return
      end subroutine incupd_init

      subroutine incupd_rd(dtime0)
!
      real*8 dtime0
!
! --- subroutine used to calculate increment field for the incremental updating
! --- version: dec 2005
!
      integer i,j,k
      logical lopen
!
      if     (ncountd.gt.incupf) then
        if     (ncountd.eq.incupf+1) then
!         should never get here, see ncountd=incupf+99 in incupd_read_*
          if (mnproc.eq.1) then
          write(lp,*) '... ended updating fields with increments ....'
          write(lp,*) 'ncountd= ',ncountd
          write(lp,*)
          endif !1st tile
          call xcsync(flush_lp)
        endif
        return
      endif
!
! --- read the target fields, and initialize the "inc" arrays.
!
      if     (incflg.gt.0) then
        call incupd_read_full(dtime0)
      else  !(incflg.le.0)
        call incupd_read_diff(dtime0)
      endif
!
      return
      end subroutine incupd_rd

      subroutine incupd(n, restrt)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
!
      logical restrt
      integer n
!
!**********
!*
! 1)  update hycom variables with increments.
!
! 2)  parameters:
!
!     output:
!      incremental updated model variables
!
! 4)  Ole Martin Smedstad (PSI), December 2005
!
!**********
!
      logical, parameter :: lpipe_incupd=.false.  !extra checking
!
      character utxt*12,vtxt*12
      integer   i,j,k
      real      q,utotij,vtotij
!
# include "stmt_fns.h"
!
      if     (ncountd.gt.incupf) then
!-------if (mnproc.eq.1) then
!-------write(lp,*) '... ended updating fields with increments .....'
!-------write(lp,*) 'ncountd= ',ncountd
!-------write(lp,*)
!-------endif !1st tile
!-------call xcsync(flush_lp)
        return
      endif
!
! --- update counter
!
      if     (incstp.ne.1) then
        ncount=ncount+1
        if     (ncount.eq.1 .or. ncount.eq.incstp+1) then
          q = 0.5  !half increment on 1st and incstp+1th step
        else
          q = 1.0
        endif !ncount
      else
        q = 1.0
      endif
!
      if     (ncount.gt.incstp+1) then
        if     (ncount.eq.incstp+2) then
          if (mnproc.eq.1) then
          write(lp,*) '... ended updating fields with increments ...'
          write(lp,*) 'ncount= ',ncount
          write(lp,*)
          endif !1st tile
          call xcsync(flush_lp)
        endif !ncount==incstp+2
        return
      endif !ncount>incstp+1
!
! --- ncount <= incstp+1
!
      if (mnproc.eq.1) then
      write(lp,*)
      if     (abs(incflg).eq.1) then
        write(lp,'(2a)') 'update fields with increments, ', &
                         'but not ubavg and vbavg'
      else   !abs(incflg).eq.2
        write(lp,'(2a)') 'update fields with increments, ', &
                         'including ubavg and vbavg'
      endif !incflg
      write(lp,*) '..........ncount,q= ',ncount,q
      endif !1st tile
      call xcsync(flush_lp)
!
! --- incremental update of dp (dpu, dpv).
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
            do k=1,kk-1
! ---         dp must be non-negative
              dp(i,j,k,n) = max( dp(i,j,k,n) + q*dpinc(i,j,k), 0.0 )
! ---          p must be at or above the bottom
               p(i,j,k+1) = min( p(i,j,k) + dp(i,j,k,n), &
                                 p(i,j,kk+1) )
              dp(i,j,k,n) =      p(i,j,k+1) - p(i,j,k)
            enddo !k
! ---       layer kk always touches the bottom
            dp(i,j,kk,n) = p(i,j,kk+1) - p(i,j,kk)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                  dpv(1-nbdy,1-nbdy,1,n), &
                  p,depthu,depthv, 0,0)
!
      if     (lpipe_incupd) then
        do k= 1,kk
          write (utxt,'(a9,i3)') 'up dpu k=',k
          write (vtxt,'(a9,i3)') 'up dpv k=',k
          call pipe_compare_sym2(dpu(1-nbdy,1-nbdy,1,n),iu,utxt, &
                                 dpv(1-nbdy,1-nbdy,1,n),iv,vtxt)
        enddo !k
      endif !lpipe_incupd
!
! --- incremental update of the other fields.
! --- salinity from updated th&S.
! --- rebalance u and v via utotij and vtotij.
!
!$OMP PARALLEL DO PRIVATE(j,i,k,utotij,vtotij) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
            do k=1,kk
              if     (tinc(i,j,k).ne.0.0 .or. &
                      sinc(i,j,k).ne.0.0     ) then
                temp(i,j,k,n) = temp(i,j,k,n) + q*tinc(i,j,k)
                saln(i,j,k,n) = saln(i,j,k,n) + q*sinc(i,j,k)
                th3d(i,j,k,n) = sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
              endif !non-zero increment
            enddo ! k
          endif !ip
          if (SEA_U) then
            utotij = 0.0
            do k=1,kk
              u(i,j,k,n) = u(i,j,k,n) + q*uinc(i,j,k)
              utotij = utotij + u(i,j,k,n)*dpu(i,j,k,n)
            enddo ! k
            utotij=utotij/depthu(i,j)
            do k=1,kk
              u(i,j,k,n) = u(i,j,k,n) - utotij
            enddo ! k
            if     (abs(incflg).eq.2) then !update ubavg
              ubavg(i,j,n) = ubavg(i,j,n) + q*ubinc(i,j)
!             ubavg(i,j,n) = ubavg(i,j,n) + q*ubinc(i,j) + utotij
            endif !incflg==2
          endif !iu
          if (SEA_V) then
            vtotij = 0.0
            do k=1,kk
              v(i,j,k,n) = v(i,j,k,n) + q*vinc(i,j,k)
              vtotij = vtotij + v(i,j,k,n)*dpv(i,j,k,n)
            enddo ! k
            vtotij=vtotij/depthv(i,j)
            do k=1,kk
              v(i,j,k,n) = v(i,j,k,n) - vtotij
            enddo ! k
            if     (abs(incflg).eq.2) then !update vbavg
              vbavg(i,j,n) = vbavg(i,j,n) + q*vbinc(i,j)
!             vbavg(i,j,n) = vbavg(i,j,n) + q*vbinc(i,j) + vtotij
            endif !incflg==2
          endif !iv
        enddo !i
      enddo ! j
!$OMP END PARALLEL DO
!
      if     (iniflg.lt.0 .and. incice.eq.1 .and. ncount.eq.incstp) then
        if (mnproc.eq.1) then
        write(lp,*) '... direct insertion of covice and thkice'
        write(lp,*)
        endif !1st tile
        call xcsync(flush_lp)
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              covice(i,j)=max(0.0,min(1.0,covice(i,j)+covinc(i,j)))
              thkice(i,j)=max(0.0,        thkice(i,j)+thkinc(i,j) )
            endif !ip
          enddo !li
        enddo !j
!$OMP   END PARALLEL DO
      endif !covice and thkice
!
      if (mnproc.eq.1) then
       write(lp,*) 'finished incupdate',ncount
       write(lp,*)
      endif !1st tile
      call xcsync(flush_lp)
!
      if     (lrestart .and. ncount.eq.incstp) then
        restrt = .true.
!
        if (mnproc.eq.1) then
         write(lp,*) 'incupd: set the flag to write a restart'
         write(lp,*)
        endif !1st tile
        call xcsync(flush_lp)
      endif
      return
      end subroutine incupd

      subroutine incupd_read_full(dtime)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
      use mod_pipe       ! HYCOM debugging interface
!
      real*8    dtime
!
! --- input 3-d HYCOM fields (from an archive file) on model day dtime.
! --- directly insert the input covice and thkice (if they exist).
! --- calculate the increment between the input and the initial state.
!
! --- filenames incup/incupd.iyear_iday_ihour.[ab].
! --- I/O and array I/O unit 925 used here, but not reserved.
!
      logical, parameter :: ldebug_incupd_read=.false. !usually .false.
!
      character flnm*24, cline*80, cvarin*6, cfield*8
      character ptxt*12,utxt*12,vtxt*12
      integer   i,idmtst,ios,j,jdmtst,k,l,layer,nskip
      integer   iyear,iday,ihour
      logical   nodens
      real      tincstp
!     real      sumdpi
!
      integer   nstep0
      real*8    dtime0
!
# include "stmt_fns.h"
!
      call forday(dtime, yrflag, iyear,iday,ihour)
!
      write(flnm,'("incup/incupd.",i4.4,"_",i3.3,"_",i2.2)') &
                                 iyear,iday,ihour
!
      if(dtime.ge.dtimeu) then
!
      ncountd=ncountd+1
      ncount=0
!
      if     (ncountd.gt.incupf) then
        if     (ncountd.eq.incupf+1) then
          if (mnproc.eq.1) then
          write(lp,*) '... ended updating fields with increments ...'
          write(lp,*) 'ncountd= ',ncountd
          write(lp,*)
          endif !1st tile
          call xcsync(flush_lp)
          ncountd=incupf+99 !turn off "ended" printout
        endif
        return
      endif !ncountd>incupf
!
      if (mnproc.eq.1) then
      write(lp,*) 'read (full) incremental updating ...'
      write(lp,*) 'ncountd ...',ncountd
      write (lp,*) 'incupd_read: ',flnm
      write (lp,*) '       time: ',dtime
      write (lp,*) 'iyear,iday,ihour: ',iyear,iday,ihour
      endif !1st tile
      call xcsync(flush_lp)
!
      call zaiopf(flnm//'.a','old', 925)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+925,file=flnm//'.b',form='formatted', &
              status='old',action='read')
!
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
!
        read(uoff+925,'(a)') cline  !'iversn'
        read(uoff+925,'(a)') cline  !'iexpt '
        read(uoff+925,'(a)') cline  !'yrflag'
      endif !1st tile
!
      call zagetc(cline,ios, uoff+925)
      read(cline,*) idmtst,cvarin
!     if     (mnproc.eq.1) then
!     write(lp,*) cvarin,' = ',idmtst
!     endif !1st tile
      if (cvarin.ne.'idm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_full - input ',cvarin, &
                              ' but should be idm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_full)')
               stop '(incupd_read_full)'
      endif
      call zagetc(cline,ios, uoff+925)
      read(cline,*) jdmtst,cvarin
!     if     (mnproc.eq.1) then
!     write(lp,*) cvarin,' = ',jdmtst
!     endif !1st tile
      if (cvarin.ne.'jdm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_full - input ',cvarin, &
                              ' but should be jdm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_full)')
               stop '(incupd_read_full)'
      endif
!
      if (idmtst.ne.itdm .or. jdmtst.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_full - input idm,jdm', &
                              ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',itdm,  jtdm,  '  (dimensions.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_full)')
               stop '(incupd_read_full)'
      endif
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        read (uoff+925,*)
      endif
!
! --- skip (most) surface fields.
!
      call zaiosk(925)
      call zagetc(cline,ios, uoff+925)
      i = index(cline,'=')
      read(cline(i+1:),*) nstep0,dtime0,layer
      if     (mnproc.eq.1) then
        write(lp,*) 'dtime0= ',dtime0
      endif
      if (dtime0.ne.dtime) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_full - input ',dtime0, &
                            ' but dtime should be ',dtime
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_full)')
               stop '(incupd_read_full)'
      endif
      nodens = layer.ne.0  !new or original archive type
      if     (nodens .and. layer.ne.sigver) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_full - input ',layer, &
                           ' sigver but should be ',sigver
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_full)')
               stop '(incupd_read_full)'
      endif
!
! assumes that there is a new incremental updating file once a day
! for "incupf" days, see blkdat.input
!
      dtimeu=dtime0+1.d0
!
      if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'dtime, dtime0, dtimeu = ',dtime, &
                     dtime0, dtimeu
        write(lp,*)
      endif !1st tile
      call xcsync(flush_lp)
!
      if     (nodens) then
        do i= 2,6
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo
      else
        do i= 2,11
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo
      endif
!
      call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice or mix_dpth or kemix
      if     (cfield.eq.'mix_dpth' .or. cfield.eq.'kemix') then
! ---   archive contains 'steric  '
        call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice
      endif
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
      if     (cfield.eq.'covice  ') then
!
! ---   directly insert covice and thkice.
!
        call rd_archive(util5, cfield,layer, 925)  !thkice
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              covice(i,j)=ubinc(i,j)
              thkice(i,j)=util5(i,j)
            endif !ip
          enddo !li
        enddo !j
!$OMP   END PARALLEL DO
        call zaiosk(925)  !temice
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call rd_archive(ubinc, cfield,layer, 925)
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
        incice =  1  !have     input covice, don't direct insert si_c
      else
        incice = -1  !have not input covice, might direct insert si_c
      endif
      call rd_archive(vbinc, cfield,layer, 925)
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
!
           if     (mnproc.eq.1) then
           write (lp,*) 'start 3-D archive file read'
           endif
           call xcsync(flush_lp)
!
! --- 3-d fields.
!
      nskip = 0
      do k=1,kk
        call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'u-vel.  ' .and. k.ne.2) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_full - expected ','u-vel.  '
          endif !1st tile
          call xcstop('(incupd_read_full)')
                 stop '(incupd_read_full)'
        elseif (cfield.ne.'u-vel.  ') then !k==2
!
! ---     count "tracer" fields (to be skipped)
!
          if     (mnproc.eq.1) then
          write(lp,'(2a)') "counting tracers: ",cfield
          endif
          do nskip= 2,99
            call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
            if     (mnproc.eq.1) then
            write(lp,'(2a)') "counting tracers: ",cfield
            endif
            if     (cfield.eq.'u-vel.  ') then
              exit
            endif
          enddo !nskip
          nskip = nskip - 1
          write(lp,'(a,i3)') "nskip =",nskip
        endif
        call rd_archive(vinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'v-vel.  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_full - expected ','v-vel.  '
          endif !1st tile
          call xcstop('(incupd_read_full)')
                 stop '(incupd_read_full)'
        endif
!          if     (mnproc.eq.1) then
!            write (lp,*) 'read v-vel archive file'
!          endif
        call rd_archive(dpinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'thknss  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                   'error in incupd_read_full - expected ','thknss  '
          endif !1st tile
          call xcstop('(incupd_read_full)')
                 stop '(incupd_read_full)'
        endif
!          if     (mnproc.eq.1) then
!            write (lp,*) 'read dpinc archive file'
!          endif
        call rd_archive(tinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'temp    ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_full - expected ','temp    '
          endif !1st tile
          call xcstop('(incupd_read_full)')
                 stop '(incupd_read_full)'
        endif
        call rd_archive(sinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'salin   ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_full - expected ','salin   '
          endif !1st tile
          call xcstop('(incupd_read_full)')
                 stop '(incupd_read_full)'
        endif
        if     (.not. nodens) then
! ---     skip density
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        endif !nodens:else
!
! ---   skip (nskip) tracers
!
        do l= 1,nskip
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo !l
      enddo !k
!
      call xctilr(dpinc,1,kk, 1,1, halo_ps)
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      close( unit=uoff+925)
      endif
      call zaiocl(925)
!
! --- calculate increments
! --- the "inc" reads, above, are full HYCOM fields (not increments yet).
!
      if(incstp.eq.1) then
        tincstp=1.0
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'tincstp=1.0        ',tincstp,incstp
        endif
      else
! ---   allow for LeapFrog and Asselin filter
        if     (mod(incstp,2).eq.1) then
          incstp = incstp-1  !must be even
        endif
        tincstp=(2.0/sqrt(1.0-ra2fac))/real(incstp)
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,'(a,f6.4,a,f12.7,i5)') &
          'tincstp=',2.0/sqrt(1.0-ra2fac),'/incstp ', &
           tincstp,incstp
        endif
      endif !incstp
!
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,*) 'calculate t,s,u,v and dp increments'
      endif !1st tile
      call xcsync(flush_lp)
!
!     if     (iutest.gt.0 .and. jutest.gt.0) then
!       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
!       write(lp,*) '*********** dpinc input ************'
!               sumdpi=0.0
!             write(lp,'(a)')
!    &                'k,dp1,dp2,dpinc='
!               do k= 1,kk
!                sumdpi=sumdpi+dpinc(iutest,jutest,k)
!                   write(lp,'(a,i3,3f20.5)')
!    &                'k= ',
!    &                 k,dp(iutest,jutest,k,1)*qonem,
!    &                 dp(iutest,jutest,k,2)*qonem,
!    &                 dpinc(iutest,jutest,k)*qonem
!                   call flush(lp)
!               enddo !k
!            write(lp,*) 'sumdpi= ', sumdpi*qonem
!            call flush(lp)
!     endif
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
          do i=1,ii
            if (SEA_U) then
            ubinc(i,j)=(ubinc(i,j) - ubavg(i,j,1))*tincstp
            do k=1,kk
!             use an approximate 2*dpu
              if     (dpinc(i,j,k)+dpinc(i-1,j,k).gt.2.0*onemm) then
                uinc(i,j,k)=(uinc(i,j,k) - u(i,j,k,1))*tincstp
              else
                uinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          endif !iu
          if (SEA_V) then
            vbinc(i,j)=(vbinc(i,j) - vbavg(i,j,1))*tincstp
            do k=1,kk
!             use an approximate 2*dpv
              if     (dpinc(i,j,k)+dpinc(i,j-1,k).gt.2.0*onemm) then
                vinc(i,j,k)=(vinc(i,j,k) - v(i,j,k,1))*tincstp
              else  
                vinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          endif !iv
          if (SEA_P) then
            do k=1,kk
              if     (dpinc(i,j,k).gt.onemm) then
                sinc(i,j,k)=(sinc(i,j,k) - saln(i,j,k,1))*tincstp
                tinc(i,j,k)=(tinc(i,j,k) - temp(i,j,k,1))*tincstp
              else
                tinc(i,j,k)=0.0  !thin target layer
                sinc(i,j,k)=0.0  !thin target layer
              endif
              dpinc(i,j,k)=(dpinc(i,j,k) - dp(i,j,k,1))*tincstp
            enddo !k
          endif !ip
        enddo !li
      enddo !j
!$OMP END PARALLEL DO
!
      call xctilr(dpinc,1,kk, 1,1, halo_ps)
!
      if     (ldebug_incupd_read) then
         call pipe_compare_sym2(ubinc,iu,'incupd:ubinc', &
                                vbinc,iv,'incupd:vbinc')
         do k= 1,kk
           write (utxt,'(a9,i3)') '  uinc k=',k
           write (vtxt,'(a9,i3)') '  vinc k=',k
           call pipe_compare_sym2(uinc(1-nbdy,1-nbdy,k),iu,utxt, &
                                  vinc(1-nbdy,1-nbdy,k),iv,vtxt)
           write (ptxt,'(a9,i3)') ' dpinc k=',k
           call pipe_compare_sym1(dpinc(1-nbdy,1-nbdy,k),ip,ptxt)
           write (ptxt,'(a9,i3)') '  tinc k=',k
           call pipe_compare_sym1( tinc(1-nbdy,1-nbdy,k),ip,ptxt)
           write (ptxt,'(a9,i3)') '  sinc k=',k
           call pipe_compare_sym1( sinc(1-nbdy,1-nbdy,k),ip,ptxt)
         enddo !k
       endif !ldebug_incupd_read
!
!     if     (iutest.gt.0 .and. jutest.gt.0) then
!       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
!       write(lp,*) '*********** dpinc out ************'
!             write(lp,'(a)')
!    &                'k,dp1,dp2,dpinc='
!               sumdpi=0.0
!               do k= 1,kk
!                sumdpi=sumdpi+dpinc(iutest,jutest,k)
!                   write(lp,'(a,i3,3f20.5)')
!    &                'k= ',
!    &                 k,dp(iutest,jutest,k,1)*qonem,
!    &                 dp(iutest,jutest,k,2)*qonem,
!    &                 dpinc(iutest,jutest,k)*qonem
!                   call flush(lp)
!               enddo !k
!            write(lp,*) 'inc sumdpi= ', sumdpi*qonem
!            call flush(lp)
!     endif
!
      if (mnproc.eq.1) then
       write(lp,*) '... finished reading incupd',dtime,dtime0
      endif !1st tile
      call xcsync(flush_lp)
!
      endif ! dtime
!
      return
      end subroutine incupd_read_full

      subroutine incupd_read_diff(dtime)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
      use mod_pipe       ! HYCOM debugging interface
!
      real*8    dtime
!
! --- input 3-d HYCOM increment fields (from a difference archive file)
! --- on model day dtime or dtime + insertion period.
! --- read in covice and thkice increments, for later direct insertion.
!
! --- filenames incup/incupd.iyear_iday_ihour.[ab].
! --- I/O and array I/O unit 925 used here, but not reserved.
!
      logical, parameter :: ldebug_incupd_read=.false. !usually .false.
!
      character flnm*24, cline*80, cvarin*6, cfield*8
      character ptxt*12,utxt*12,vtxt*12
      integer   i,idmtst,ios,j,jdmtst,k,l,layer,nskip
      integer   iyear,iday,ihour
      real      tincstp
!     real      sumdpi
!
      logical   lexist
      integer   nstep0
      real*8    dtimei,dtime0
!
# include "stmt_fns.h"
!
      dtimei = dtime
!
      call forday(dtimei, yrflag, iyear,iday,ihour)
!
      write(flnm,'("incup/incupd.",i4.4,"_",i3.3,"_",i2.2)') &
                                 iyear,iday,ihour
!
      inquire(file=flnm//'.a',exist=lexist)
!
! --- is input at end of insertion period
!
      if     (.not. lexist) then
        dtimei = dtime + incstp*(baclin/86400.d0)
!
        call forday(dtimei, yrflag, iyear,iday,ihour)
!
        write(flnm,'("incup/incupd.",i4.4,"_",i3.3,"_",i2.2)') &
                                   iyear,iday,ihour
      endif
!
      if(dtime.ge.dtimeu) then
!
      ncountd=ncountd+1
      ncount=0
!
      if     (ncountd.gt.incupf) then
        if     (ncountd.eq.incupf+1) then
          if (mnproc.eq.1) then
          write(lp,*) '... ended updating fields with increments ...'
          write(lp,*) 'ncountd= ',ncountd
          write(lp,*)
          endif !1st tile
          call xcsync(flush_lp)
          ncountd=incupf+99 !turn off "ended" printout
        endif
        return
      endif !ncountd>incupf
!
      if (mnproc.eq.1) then
      write(lp,*) 'read (diff) incremental updating ...'
      write(lp,*) 'ncountd ...',ncountd
      write (lp,*) 'incupd_read: ',flnm
      if     (dtime.eq.dtimei) then
        write (lp,*) '       time: ',dtime
      else
        write (lp,*) '       time: ',dtime
        write (lp,*) '    time_in: ',dtimei
      endif
      write (lp,*) 'iyear,iday,ihour: ',iyear,iday,ihour
      endif !1st tile
      call xcsync(flush_lp)
!
      call zaiopf(flnm//'.a','old', 925)
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        open (unit=uoff+925,file=flnm//'.b',form='formatted', &
              status='old',action='read')
!
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
        read(uoff+925,'(a)') cline
!
        read(uoff+925,'(a)') cline  !'iversn'
        read(uoff+925,'(a)') cline  !'jexpt '
        read(uoff+925,'(a)') cline  !'iexpt '
        read(uoff+925,'(a)') cline  !'yrflag'
      endif !1st tile
!
      call zagetc(cline,ios, uoff+925)
      read(cline,*) idmtst,cvarin
!     if     (mnproc.eq.1) then
!     write(lp,*) cvarin,' = ',idmtst
!     endif !1st tile
      if (cvarin.ne.'idm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_diff - input ',cvarin, &
                              ' but should be idm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_diff)')
               stop '(incupd_read_diff)'
      endif
      call zagetc(cline,ios, uoff+925)
      read(cline,*) jdmtst,cvarin
!     if     (mnproc.eq.1) then
!     write(lp,*) cvarin,' = ',jdmtst
!     endif !1st tile
      if (cvarin.ne.'jdm   ') then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_diff - input ',cvarin, &
                              ' but should be jdm   '
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_diff)')
               stop '(incupd_read_diff)'
      endif
!
      if (idmtst.ne.itdm .or. jdmtst.ne.jtdm) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_diff - input idm,jdm', &
                              ' not consistent with parameters'
        write(lp,*) 'idm,jdm = ',itdm,  jtdm,  '  (dimensions.h)'
        write(lp,*) 'idm,jdm = ',idmtst,jdmtst,'  (input)'
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_diff)')
               stop '(incupd_read_diff)'
      endif
!
! --- check for a difference archive
!
      call zagetc(cline,ios, uoff+925)
      i = index(cline,'diff day')
      if     (i.eq.0) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_diff - ', &
                              'not a difference archive'
        write(lp,*) 'header = '
        write(lp,*)  trim(cline)
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_diff)')
               stop '(incupd_read_diff)'
      endif
!
! --- skip (most) surface fields.
!
      call zaiosk(925)
      call zagetc(cline,ios, uoff+925)
      i = index(cline,'=')
      read(cline(i+1:),*) nstep0,dtime0,layer
      if     (mnproc.eq.1) then
        write(lp,*) 'dtime0= ',dtime0
      endif
      if (dtime0.ne.dtimei) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in incupd_read_diff - input ',dtime0, &
                                ' but dtimei should be ',dtimei
        write(lp,*)
        endif !1st tile
        call xcstop('(incupd_read_diff)')
               stop '(incupd_read_diff)'
      endif
!
! assumes that there is a new incremental updating file once a day
! for "incupf" days, see blkdat.input
!
      dtimeu=dtime0+1.d0
!
      if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'dtime, dtime0, dtimeu = ',dtime, &
                     dtime0, dtimeu
        write(lp,*)
      endif !1st tile
      call xcsync(flush_lp)
!
      do i= 2,12
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call zaiosk(925)
      enddo
!
      call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice or kemix
      if     (cfield.eq.'kemix') then
! ---   archive contains 'steric  '
        call rd_archive(ubinc, cfield,layer, 925)  !u_btrop or covice
      endif
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
      if     (cfield.eq.'covice  ') then
        covinc(:,:) = ubinc(:,:)
        call rd_archive(thkinc, cfield,layer, 925)  !thkice
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
        call zaiosk(925)  !temice
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call rd_archive(ubinc, cfield,layer, 925)
        if     (mnproc.eq.1) then
        write(lp,'(2a)') "surface: ",cfield
        endif
        call xcsync(flush_lp)
        incice =  1  !have     input covice, don't direct insert si_c
      else
        incice = -1  !have not input covice, might direct insert si_c
      endif
      call rd_archive(vbinc, cfield,layer, 925)
      if     (mnproc.eq.1) then
      write(lp,'(2a)') "surface: ",cfield
      endif
      call xcsync(flush_lp)
! --- kebtrop
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
        read (uoff+925,*)
      endif
      call zaiosk(925)
!
           if     (mnproc.eq.1) then
           write (lp,*) 'start 3-D archive file read'
           endif
           call xcsync(flush_lp)
!
! --- 3-d fields.
!
      nskip = 0
      do k=1,kk
        call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'u-vel.  ' .and. k.ne.2) then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_diff - expected ','u-vel.  '
          endif !1st tile
          call xcstop('(incupd_read_diff)')
                 stop '(incupd_read_diff)'
        elseif (cfield.ne.'u-vel.  ') then !k==2
!
! ---     count "tracer" fields (to be skipped)
!
          if     (mnproc.eq.1) then
          write(lp,'(2a)') "counting tracers: ",cfield
          endif
          do nskip= 2,99
            call rd_archive(uinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
            if     (mnproc.eq.1) then
            write(lp,'(2a)') "counting tracers: ",cfield
            endif
            if     (cfield.eq.'u-vel.  ') then
              exit
            endif
          enddo !nskip
          nskip = nskip - 1
          write(lp,'(a,i3)') "nskip =",nskip
        endif
        call rd_archive(vinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'v-vel.  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_diff - expected ','v-vel.  '
          endif !1st tile
          call xcstop('(incupd_read_diff)')
                 stop '(incupd_read_diff)'
        endif
!          if     (mnproc.eq.1) then
!            write (lp,*) 'read v-vel archive file'
!          endif
! ---   skip k.e.
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call zaiosk(925)
! ---   skip mnthknss
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call zaiosk(925)
!
        call rd_archive(dpinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'thknss  ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                   'error in incupd_read_diff - expected ','thknss  '
          endif !1st tile
          call xcstop('(incupd_read_diff)')
                 stop '(incupd_read_diff)'
        endif
!          if     (mnproc.eq.1) then
!            write (lp,*) 'read dpinc archive file'
!          endif
        call rd_archive(tinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'temp    ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_diff - expected ','temp    '
          endif !1st tile
          call xcstop('(incupd_read_diff)')
                 stop '(incupd_read_diff)'
        endif
        call rd_archive(sinc(1-nbdy,1-nbdy,k), cfield,layer, 925)
        if     (cfield.ne.'salin   ') then
          if     (mnproc.eq.1) then
          write(lp,'(/ a / a,a /)') cfield, &
                 'error in incupd_read_diff - expected ','salin   '
          endif !1st tile
          call xcstop('(incupd_read_diff)')
                 stop '(incupd_read_diff)'
        endif
! ---   skip density
        if     (mnproc.eq.1) then  ! .b file from 1st tile only
          read (uoff+925,*)
        endif
        call zaiosk(925)
!
! ---   skip (nskip) tracers
!
        do l= 1,nskip
          if     (mnproc.eq.1) then  ! .b file from 1st tile only
            read (uoff+925,*)
          endif
          call zaiosk(925)
        enddo !l
      enddo !k
!
      call xctilr(dpinc,1,kk, 1,1, halo_ps)
!
      if     (mnproc.eq.1) then  ! .b file from 1st tile only
      close( unit=uoff+925)
      endif
      call zaiocl(925)
!
! --- calculate increments
! --- the "inc" reads, above, are HYCOM diff fields.
!
      if(incstp.eq.1) then
        tincstp=1.0
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'tincstp=1.0        ',tincstp,incstp
        endif
      else
! ---   allow for LeapFrog and Asselin filter
        if     (mod(incstp,2).eq.1) then
          incstp = incstp-1  !must be even
        endif
        tincstp=(2.0/sqrt(1.0-ra2fac))/real(incstp)
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,'(a,f6.4,a,f12.7,i5)') &
          'tincstp=',2.0/sqrt(1.0-ra2fac),'/incstp ', &
           tincstp,incstp
        endif
      endif !incstp
!
      if (mnproc.eq.1) then
      write(lp,*)
      write(lp,*) 'calculate t,s,u,v and dp increments'
      endif !1st tile
      call xcsync(flush_lp)
!
!     if     (iutest.gt.0 .and. jutest.gt.0) then
!       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
!       write(lp,*) '*********** dpinc input ************'
!               sumdpi=0.0
!             write(lp,'(a)')
!    &                'k,dp1,dp2,dpinc='
!               do k= 1,kk
!                sumdpi=sumdpi+dpinc(iutest,jutest,k)
!                   write(lp,'(a,i3,3f20.5)')
!    &                'k= ',
!    &                 k,dp(iutest,jutest,k,1)*qonem,
!    &                 dp(iutest,jutest,k,2)*qonem,
!    &                 dpinc(iutest,jutest,k)*qonem
!                   call flush(lp)
!               enddo !k
!            write(lp,*) 'sumdpi= ', sumdpi*qonem
!            call flush(lp)
!     endif
!
! --- use background thickness fields for thinness detection
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
          do i=1,ii
            if (SEA_U) then
            ubinc(i,j)=ubinc(i,j)*tincstp
            do k=1,kk
              if     (dpu(i,j,k,1).gt.onemm) then
                uinc(i,j,k)=uinc(i,j,k)*tincstp
              else
                uinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          endif !iu
          if (SEA_V) then
            vbinc(i,j)=vbinc(i,j)*tincstp
            do k=1,kk
              if     (dpv(i,j,k,1).gt.onemm) then
                vinc(i,j,k)=vinc(i,j,k)*tincstp
              else  
                vinc(i,j,k)=0.0  !thin target layer
              endif
            enddo !k
          endif !iv
          if (SEA_P) then
            do k=1,kk
              if     (dp(i,j,k,1).gt.onemm) then
                sinc(i,j,k)=sinc(i,j,k)*tincstp
                tinc(i,j,k)=tinc(i,j,k)*tincstp
              else
                tinc(i,j,k)=0.0  !thin target layer
                sinc(i,j,k)=0.0  !thin target layer
              endif
              dpinc(i,j,k)=dpinc(i,j,k)*tincstp
            enddo !k
          endif !ip
        enddo !li
      enddo !j
!$OMP END PARALLEL DO
!
      call xctilr(dpinc,1,kk, 1,1, halo_ps)
!
      if     (ldebug_incupd_read) then
         call pipe_compare_sym2(ubinc,iu,'incupd:ubinc', &
                                vbinc,iv,'incupd:vbinc')
         do k= 1,kk
           write (utxt,'(a9,i3)') '  uinc k=',k
           write (vtxt,'(a9,i3)') '  vinc k=',k
           call pipe_compare_sym2(uinc(1-nbdy,1-nbdy,k),iu,utxt, &
                                  vinc(1-nbdy,1-nbdy,k),iv,vtxt)
           write (ptxt,'(a9,i3)') ' dpinc k=',k
           call pipe_compare_sym1(dpinc(1-nbdy,1-nbdy,k),ip,ptxt)
           write (ptxt,'(a9,i3)') '  tinc k=',k
           call pipe_compare_sym1( tinc(1-nbdy,1-nbdy,k),ip,ptxt)
           write (ptxt,'(a9,i3)') '  sinc k=',k
           call pipe_compare_sym1( sinc(1-nbdy,1-nbdy,k),ip,ptxt)
         enddo !k
       endif !ldebug_incupd_read
!
!     if     (iutest.gt.0 .and. jutest.gt.0) then
!       write(lp,*) '*',' iutest= ',iutest+i0,' jutest= ',jutest+j0,' *'
!       write(lp,*) '*********** dpinc out ************'
!             write(lp,'(a)')
!    &                'k,dp1,dp2,dpinc='
!               sumdpi=0.0
!               do k= 1,kk
!                sumdpi=sumdpi+dpinc(iutest,jutest,k)
!                   write(lp,'(a,i3,3f20.5)')
!    &                'k= ',
!    &                 k,dp(iutest,jutest,k,1)*qonem,
!    &                 dp(iutest,jutest,k,2)*qonem,
!    &                 dpinc(iutest,jutest,k)*qonem
!                   call flush(lp)
!               enddo !k
!            write(lp,*) 'inc sumdpi= ', sumdpi*qonem
!            call flush(lp)
!     endif
!
      if (mnproc.eq.1) then
       write(lp,*) '... finished reading incupd',dtime,dtime0
      endif !1st tile
      call xcsync(flush_lp)
!
      endif ! dtime
!
      return
      end subroutine incupd_read_diff

      subroutine incupd_si_c(dtime)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_za         ! HYCOM I/O interface
      use mod_pipe       ! HYCOM debugging interface
!
      real*8    dtime
!
! --- directly insert si_c into covice and thkice.
!
      integer   i,j,l
!
      if     (incice.eq.-1) then
        incice = 0  !will have directly inserted si_c
!
! ---   directly insert covice and thkice.
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-nbdy,jj+nbdy
          do i=1,ii
            if (SEA_P) then
              covice(i,j)=si_c(i,j)
              thkice(i,j)=covice(i,j)*hicemn
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        if (mnproc.eq.1) then
         write(lp,*) '... finished inserting si_c',dtime
        endif !1st tile
        call xcsync(flush_lp)
      endif !incice
!
      return
      end subroutine incupd_si_c
 
      subroutine stfupd(n)
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
!
      integer n
!
!**********
!*
! 1)  update hycom variables with stochastic anomalies.
!
! 2)  parameters:
!
!     output:
!      updated model variables
!
! 4)  Alan Wallcaft, NRL, September 2017.
!     Based on incupd.
!
!**********
!
      logical, parameter :: lpipe_stfupd=.false.  !extra checking
!
      integer   i,j,k
      real      q,q_s,q_t,zij,utotij,vtotij
!
# include "stmt_fns.h"
!
      if     (stfflg.eq.0) then
        return
      endif
!
! --- incremental update of T&S&vel.
! --- rebalance u and v via utotij and vtotij.
!
!$OMP PARALLEL DO PRIVATE(j,i,k,utotij,vtotij) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
            zij = -0.5*dp(i,j,1,n)*qonem  !layer 1 correction
            do k=1,kk
              zij = zij + 0.5*qonem*(dp(i,j,max(k-1,1),n)+dp(i,j,k,n))  !layer k center
              q_t = (delt1/3600.0)*exp(-zij/stfrdt)
              q_s = (delt1/3600.0)*exp(-zij/stfrds)
              temp(i,j,k,n) = temp(i,j,k,n) + q_t*(stoc_t(i,j,l0)*w0+ &
                                                   stoc_t(i,j,l1)*w1 )
              saln(i,j,k,n) = saln(i,j,k,n) + q_s*(stoc_s(i,j,l0)*w0+ &
                                                   stoc_s(i,j,l1)*w1 )
              th3d(i,j,k,n) = sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
            enddo ! k
          endif !ip
          if     (stfflg.eq.2) then
          if (SEA_U) then
            utotij = 0.0
            zij = -0.5*dpu(i,j,1,n)*qonem  !layer 1 correction
            do k=1,kk
              zij = zij + 0.5*qonem*(dpu(i,j,max(k-1,1),n)+dpu(i,j,k,n))  !layer k center
              q   = (delt1/3600.0)*exp(-zij/stfrdv)
              u(i,j,k,n) = u(i,j,k,n) + q*(stoc_u(i,j,l0)*w0+ &
                                           stoc_u(i,j,l1)*w1 )
            enddo ! k
            utotij=utotij/depthu(i,j)
            do k=1,kk
              u(i,j,k,n) =   u(i,j,k,n) - utotij
            enddo ! k
            ubavg(i,j,n) = ubavg(i,j,n) + utotij
          endif !iu
          if (SEA_V) then
            vtotij = 0.0
            zij = -0.5*dpv(i,j,1,n)*qonem  !layer 1 correction
            do k=1,kk
              zij = zij + 0.5*qonem*(dpv(i,j,max(k-1,1),n)+dpv(i,j,k,n))  !layer k center
              q   = (delt1/3600.0)*exp(-zij/stfrdv)
              v(i,j,k,n) = v(i,j,k,n) + q*(stoc_v(i,j,l0)*w0+ &
                                           stoc_v(i,j,l1)*w1 )
              vtotij = vtotij + v(i,j,k,n)*dpv(i,j,k,n)
            enddo ! k
            vtotij=vtotij/depthv(i,j)
            do k=1,kk
              v(i,j,k,n) =   v(i,j,k,n) - vtotij
            enddo ! k
            vbavg(i,j,n) = vbavg(i,j,n) + vtotij
          endif !iv
          endif !stfflg==2
        enddo !i
      enddo ! j
!$OMP END PARALLEL DO
!
      return
      end subroutine stfupd
!
      end module mod_incupd
!
!
!> Revision history:
!>
!> Feb  2006 - 1st module version
!> May  2006 - changed to read multiple increment files
!> Jul  2011 - thin layer is now 1mm (no longer 1m)
!> Jul  2011 - replace thinc with sinc
!> Nov  2012 - bugfix: added xctilr(dpinc to update halo
!> Dec  2012 - cleaned up printing to .log file
!> Apr  2012 - added incice and incupd_si_c
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Nov  2015 - added iniflg<0 for reading a difference archive
!> Dec  2015 - allow for Asselin filter when calculating tincstp
!> Dec  2015 - half increment for 1st and incstp+1th time step
!> Dec  2015 - difference archive can be at end of insertion period
!> Jul  2017 - bugfix for steric in archive
!> Aug  2017 - incupf -ve to write a restart after each update completes
!> Sep  2017 - added subroutine stfupd
