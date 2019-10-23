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
      module mod_asselin
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
!
      implicit none
!
! --- module for HYCOM Robert-Asselin filter of scalar fields
!
      private !! default is private
      public  :: asselin_save,asselin_filter

      contains

      subroutine asselin_save(m,n)
!
      integer m,n
!
! --- save time level t-1 for Robert-Asselin filter
! --- dpo is initialized in cnuity
!
! --- on exit:
! ---    otemp(:,:,:)   = time step t-1
! ---    osaln(:,:,:)   = time step t-1
! ---    oth3d(:,:,:)   = time step t-1
! ---  otracer(:,:,:,:) = time step t-1
!
! ---   onetao(:,:,  n) = time step t-1
! ---   oneta( :,:,  n) = time step t-1
! ---   onetao(:,:,  m) = time step t
! ---   oneta( :,:,  m) = time step t
!
      integer i,j,k,ktr
!
!$OMP PARALLEL DO PRIVATE(j,i,k,ktr)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
          oneta( i,j,n) = max( oneta0, 1.0 + pbavg(i,j,n)/pbot(i,j) )  !t-1
          oneta( i,j,m) = max( oneta0, 1.0 + pbavg(i,j,m)/pbot(i,j) )  !t
          onetao(i,j,n) = oneta( i,j,n)
          onetao(i,j,m) = oneta( i,j,m)
          endif
        enddo !i
        do k= 1,kk
          do i=1,ii
            otemp(i,j,k)   = temp(i,j,k,n)
            osaln(i,j,k)   = saln(i,j,k,n)
            oth3d(i,j,k)   = th3d(i,j,k,n)
            do ktr= 1,ntracr
              otracer(i,j,k,ktr) = tracer(i,j,k,n,ktr)
            enddo !ktr
          enddo !i
        enddo !k
        if     (mxlmy) then
          do k= 1,kk
            do i=1,ii
               oq2(i,j,k) =  q2(i,j,k,n)
              oq2l(i,j,k) = q2l(i,j,k,n)
            enddo !i
          enddo !k
        endif !mxkmy
      enddo !j
!
      call xctilr(oneta( 1-nbdy,1-nbdy,1),1,2, 6,6, halo_ps)
      call xctilr(onetao(1-nbdy,1-nbdy,1),1,2, 6,6, halo_ps)
!
      return
      end subroutine asselin_save

      subroutine asselin_filter(m,n)
!
      integer m,n
!
! --- Apply Robert-Asselin filter
!
      logical, parameter ::  lpipe_asselin=.false.  !extra checking (when pipe on)
      logical, parameter :: ldebug_asselin=.false.  !debugging, usually false
!
      real,    parameter :: onezm=9806.e-20 !zepto-m, insignificant thickness
!
      logical          latemp,lath3d
      integer          i,j,k,ktr,margin
      real             q,dpsold,dpsmid,dpsnew, &
                       dpold,dpmid,dpmidn,dpnew,qdpmidn,smin
      double precision ssum(4),s1(4)
!
      character*12 text
!
# include "stmt_fns.h"
!
      margin = 0  !at end of time step
!
      if     (ldebug_asselin) then
        ssum(1:4) = 0.0d0
      endif !debug
!
!$OMP PARALLEL DO PRIVATE(j,k,i,ktr,dpsold,dpsmid,dpsnew,q, &
!$OMP                     dpold,dpmid,dpmidn,dpnew,qdpmidn,smin) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            oneta(i,j,n) = max( oneta0, 1.0 + pbavg(i,j,n)/pbot(i,j) )  !t+1
            oneta(i,j,m) = max( oneta0, 1.0 + pbavg(i,j,m)/pbot(i,j) )  !t with RA
          endif !ip
        enddo !i
        do k= 1,kk
          latemp =  k.le.nhybrd .and. advflg.eq.0       ! advect temp
          lath3d = (k.le.nhybrd .and. advflg.eq.1) .or. &
                   (k.eq.1      .and. isopyc     )      ! advect th3d
          do i=1-margin,ii+margin
            if (SEA_P) then
!
! ---         Robert-Asselin time filter of scalar fields
! ---         Note that this is smoothing oneta * dp * scalar
!
              dpold  = dpo(i,j,k,n)*onetao(i,j,n)  !t-1
              dpmid  = dpo(i,j,k,m)*onetao(i,j,m)  !t
              dpnew  = dp( i,j,k,n)*oneta( i,j,n)  !t+1
              q      = 0.5*ra2fac*(dpold+dpnew-2.0*dpmid)
              dpmidn = dpmid + q
              dp(i,j,k,m) = dpmidn/oneta(i,j,m)  !t with RA
              if     (dpmidn.gt.onezm) then !effectively non-zero
                if     (ldebug_asselin .and. i.eq.itest &
                                       .and. j.eq.jtest) then
                  ssum(1) = ssum(1) + dpold*osaln(i,j,k)
                  ssum(2) = ssum(2) + dpmid* saln(i,j,k,m)
                  ssum(3) = ssum(3) + dpnew* saln(i,j,k,n)
                endif !debug
!
                qdpmidn = 1.0/dpmidn
!
! ---           version that exactly conserves constant salinity
                smin   =   min( osaln(i,j,k), &
                                 saln(i,j,k,m), &
                                 saln(i,j,k,n) )
                dpsold = dpold*(osaln(i,j,k)  -smin)  !t-1
                dpsmid = dpmid*( saln(i,j,k,m)-smin)  !t
                dpsnew = dpnew*( saln(i,j,k,n)-smin)  !t+1
                q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
                saln(i,j,k,m)= smin + (dpsmid + q) * qdpmidn  !t with RA
!
                if     (ldebug_asselin .and. i.eq.itest &
                                       .and. j.eq.jtest) then
                  ssum(4) = ssum(4) + dpmidn*saln(i,j,k,m)
!                 write(6,*) k,ssum(4)*qonem
                  if     (k.eq.1) then
! ---               normalize by onem, layer 1 is often 1 m thick
                    s1(1) = ssum(1)/onem
                    s1(2) = ssum(2)/onem
                    s1(3) = ssum(3)/onem
                    s1(4) = ssum(4)/onem
                    write(lp,'(a,i9,1p4e16.8)') &
                     'asseln1:',nstep,s1(1),s1(2),s1(4),s1(3)
                  endif
                endif !debug
!
                if     (latemp) then
                  smin   =   min( otemp(i,j,k), &
                                   temp(i,j,k,m), &
                                   temp(i,j,k,n) )
                  dpsold = dpold*(otemp(i,j,k)  -smin)  !t-1
                  dpsmid = dpmid*( temp(i,j,k,m)-smin)  !t
                  dpsnew = dpnew*( temp(i,j,k,n)-smin)  !t+1
                  q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
                  temp(i,j,k,m) = smin + (dpsmid + q) * qdpmidn
! ---             update dependent thermodynamic variable
                  th3d(i,j,k,m) = sig(temp(i,j,k,m), &
                                      saln(i,j,k,m))-thbase
                elseif (lath3d) then
                  smin   =   min( oth3d(i,j,k), &
                                   th3d(i,j,k,m), &
                                   th3d(i,j,k,n) )
                  dpsold = dpold*(oth3d(i,j,k)  -smin)  !t-1
                  dpsmid = dpmid*( th3d(i,j,k,m)-smin)  !t
                  dpsnew = dpnew*( th3d(i,j,k,n)-smin)  !t+1
                  q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
                  th3d(i,j,k,m) = smin + (dpsmid + q) * qdpmidn
! ---             update dependent thermodynamic variable
                  temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase, &
                                         saln(i,j,k,m))
                else   ! exactly isopycnal layer
                  th3d(i,j,k,m) = theta(i,j,k)
! ---             update dependent thermodynamic variable
                  temp(i,j,k,m) = tofsig(th3d(i,j,k,m)+thbase, &
                                         saln(i,j,k,m))
                endif
                do ktr= 1,ntracr
! ---             non-negative version that exactly conserves constant tracers
                  smin   = min( otracer(i,j,k,  ktr), &
                                 tracer(i,j,k,m,ktr), &
                                 tracer(i,j,k,n,ktr) )
                  dpsold = dpold*(otracer(i,j,k,  ktr) - smin)  !>=0
                  dpsmid = dpmid*( tracer(i,j,k,m,ktr) - smin)  !>=0
                  dpsnew = dpnew*( tracer(i,j,k,n,ktr) - smin)  !>=0
                  q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
!diag               if (i.eq.itest.and.j.eq.jtest.and.ktr.eq.1) then
!diag                   write(lp,'(a,i3,1p3g15.6)') &
!diag                     'RA dpo:',k,dpold*qonem,dpmid*qonem,dpnew*qonem
!diag                   write(lp,'(a,i3,3f15.6)') &
!diag                     'RA tro:',k,otracer(i,j,k,  ktr), &
!diag                                  tracer(i,j,k,m,ktr), &
!diag                                  tracer(i,j,k,n,ktr)
!diag                   write(lp,'(a,i3,1p2g15.6)') &
!diag                     'RA dpm:',k,dpmid*qonem, &
!diag                                 dp(i,j,k,m)*qonem
!diag                   write(lp,'(a,i3,2f15.6)') &
!diag                     'RA trm:',k,tracer(i,j,k,m,ktr), &
!diag                                 smin + (dpsmid + q) * qdpmidn
!diag                   call flush(lp)
!diag               endif !debug
                  tracer(i,j,k,m,ktr) = smin + (dpsmid + q) * qdpmidn
                enddo !ktr
                if     (mxlmy) then
                  dpsold = dpold*oq2(i,j,k)
                  dpsmid = dpmid* q2(i,j,k,m)
                  dpsnew = dpnew* q2(i,j,k,n)
                  q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
                  q2(i,j,k,m) = (dpsmid + q) * qdpmidn
                  dpsold = dpold*oq2l(i,j,k)
                  dpsmid = dpmid* q2l(i,j,k,m)
                  dpsnew = dpnew* q2l(i,j,k,n)
                  q      = 0.5*ra2fac*(dpsold+dpsnew-2.0*dpsmid)
                  q2l(i,j,k,m) = (dpsmid + q) * qdpmidn
                endif !mxlmy
              endif !effectively non-zero 
!
              if     (ldebug_asselin .and. k.eq.kk .and. &
                          i.eq.itest .and. j.eq.jtest   ) then
! ---           normalize by initial depth
                ssum(1) = ssum(1)/pbot(i,j)
                ssum(2) = ssum(2)/pbot(i,j)
                ssum(3) = ssum(3)/pbot(i,j)
                ssum(4) = ssum(4)/pbot(i,j)
                write(lp,'(a,i9,1p4e16.8)') &
                 'asselin:',nstep,ssum(1),ssum(2),ssum(4),ssum(3)
              endif !debug
!
            endif !ip
          enddo !i
        enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_asselin) then
! ---   compare two model runs.
        do k= 1,kk
          write (text,'(a9,i3)') 'otemp  k=',k
          call pipe_compare_sym1(otemp(1-nbdy,1-nbdy,k),  ip,text)
          write (text,'(a9,i3)') 'temp.m k=',k
          call pipe_compare_sym1( temp(1-nbdy,1-nbdy,k,m),ip,text)
          write (text,'(a9,i3)') 'temp.n k=',k
          call pipe_compare_sym1( temp(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'osaln  k=',k
          call pipe_compare_sym1(osaln(1-nbdy,1-nbdy,k),  ip,text)
          write (text,'(a9,i3)') 'saln.m k=',k
          call pipe_compare_sym1( saln(1-nbdy,1-nbdy,k,m),ip,text)
          write (text,'(a9,i3)') 'saln.n k=',k
          call pipe_compare_sym1( saln(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'oth3d  k=',k
          call pipe_compare_sym1(oth3d(1-nbdy,1-nbdy,k),  ip,text)
          write (text,'(a9,i3)') 'th3d.m k=',k
          call pipe_compare_sym1( th3d(1-nbdy,1-nbdy,k,m),ip,text)
          write (text,'(a9,i3)') 'th3d.n k=',k
          call pipe_compare_sym1( th3d(1-nbdy,1-nbdy,k,n),ip,text)
        enddo !k
      endif
!
!     call pipe_comparall(m,n, 'asseln, step')
!
      return
      end subroutine asselin_filter

      end module mod_asselin
!
!
!> Revision history:
!>
!> Aug. 2018 - 1st version
!> Feb. 2019 - replaced onetai by 1.0
!> Sep. 2019 - added oneta0
