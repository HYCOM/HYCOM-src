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
      subroutine thermf_oi(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n
!
! --- ----------------------------------------------------------
! --- thermal forcing - combine ocean and sea ice surface fluxes
! ---                 - complete surface salinity forcing
! ---                 - allow for ice shelves (no surface flux)
! --- ----------------------------------------------------------
!
      logical, parameter ::  ldebug_empbal=.false.   !usually .false.
!
      integer i,j,k
      real    emnp,dpemnp,dplay1,onetanew,onetaold, &
              pbanew,pbaold,q,salt1n,salt1o
      real*8  d1,d2,d3,d4,ssum(2),s1(2)
!
      if     (ishelf.ne.0) then  !sea ice and an ice shelf
!$OMP PARALLEL DO PRIVATE(j,i)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              if     (ishlf(i,j).eq.1) then  !standard ocean point
                if ( cpl_swflx  .and. cpl_lwmdnflx .and. cpl_lwmupflx &
                                .and. cpl_precip    ) then

                  sstflx(i,j) = (1.0-covice(i,j))*sstflx(i,j)   !relax over ocean
                  surflx(i,j) =     surflx(i,j) + flxice(i,j)   !ocn/ice frac. in coupler
                  salflx(i,j) =                   sflice(i,j) +  & !ocn/ice frac. in coupler
                                                  sssflx(i,j)   !relax  evrywhere
                  wtrflx(i,j) =     wtrflx(i,j) + wflice(i,j) +  & !ocn/ice frac. in coupler
                                                  rivflx(i,j)   !rivers evrywhere 
                else

                  sswflx(i,j) = (1.0-covice(i,j))*sswflx(i,j) +  & !ocean
                                                  fswice(i,j)   !ice cell average
                  surflx(i,j) = (1.0-covice(i,j))*surflx(i,j) +  & !ocean
                                                  flxice(i,j)   !ice cell average
                  sstflx(i,j) = (1.0-covice(i,j))*sstflx(i,j)   !relax over ocean
                  salflx(i,j) =                   sflice(i,j) +  & !ice cell average
                                                  sssflx(i,j)   !relax everywhere
                  wtrflx(i,j) = (1.0-covice(i,j))*wtrflx(i,j) +  & !ocean
                                                  wflice(i,j) +  & !ice cell average
                                                  rivflx(i,j)   !rivers everywhere
                endif ! .not.cpl
              else  !under an ice shelf
                sswflx(i,j) = 0.0
                surflx(i,j) = 0.0
                sstflx(i,j) = 0.0
                salflx(i,j) = 0.0
                wtrflx(i,j) = 0.0
              endif  !ishlf
              util1(i,j) = surflx(i,j)*scp2(i,j)
              util2(i,j) = wtrflx(i,j)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP END PARALLEL DO
      elseif (iceflg.ne.0) then
          if ( cpl_swflx  .and. cpl_lwmdnflx .and. cpl_lwmupflx &
                          .and. cpl_precip    ) then
! ---     allow for sea ice
!$OMP PARALLEL DO PRIVATE(j,i)
            do j=1,jj
               do i=1,ii
                  if (SEA_P) then
                     sstflx(i,j) = (1.0-covice(i,j))*sstflx(i,j)   !relax over ocean
                     surflx(i,j) =     surflx(i,j) + flxice(i,j)   !ocn/ice frac. in coupler
                     salflx(i,j) =                   sflice(i,j) +  & !ocn/ice frac. in coupler
                                                     sssflx(i,j)   !relax  evrywhere
                     wtrflx(i,j) =     wtrflx(i,j) + wflice(i,j) +  & !ocn/ice frac. in coupler
                                                     rivflx(i,j)   !rivers evrywhere
                     util1(i,j) = surflx(i,j)*scp2(i,j)
                     util2(i,j) = wtrflx(i,j)*scp2(i,j)
                  endif
               enddo !i
            enddo
!$OMP END PARALLEL DO
          else ! cpl_
!$OMP PARALLEL DO PRIVATE(j,i)
            do j=1,jj
               do i=1,ii
                 if (SEA_P) then
                  sswflx(i,j) = (1.0-covice(i,j))*sswflx(i,j) + &
                                                  fswice(i,j)   !ice cell average
                  surflx(i,j) = (1.0-covice(i,j))*surflx(i,j) + &
                                                  flxice(i,j)   !ice cell average
                  sstflx(i,j) = (1.0-covice(i,j))*sstflx(i,j)   !relax over ocean
                  salflx(i,j) =                   sflice(i,j) +  & !ice cell average
                                                  sssflx(i,j)   !relax everywhere
                  wtrflx(i,j) = (1.0-covice(i,j))*wtrflx(i,j) + &
                                                  wflice(i,j) +  & !ice cell average
                                                  rivflx(i,j)   !rivers everywhere
                  util1(i,j) = surflx(i,j)*scp2(i,j)
                  util2(i,j) = wtrflx(i,j)*scp2(i,j)
                  endif
             enddo !i
           enddo !j
!$OMP END PARALLEL DO
          endif ! .not. cpl
      else !no sea ice
!$OMP PARALLEL DO PRIVATE(j,i)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then  !sswflx, surflx and sstflx unchanged
              salflx(i,j) = sssflx(i,j)
              wtrflx(i,j) = wtrflx(i,j) + rivflx(i,j)
               util1(i,j) = surflx(i,j)*scp2(i,j)
               util2(i,j) = wtrflx(i,j)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP END PARALLEL DO
      endif !ishlf:iceflg
!
      if     (epmass .and. empbal.eq.1) then
!
! ---   balance wtrflx to emptgt, via a region-wide offset:
!
        call xcsum(d2, util2,ipa)  !total wtrflx
        d2 = d2/area               !basin-wide average
        if     (d2.ne.emptgt) then  !not balanced at emptgt
          q = emptgt - d2
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                wtrflx(i,j) = wtrflx(i,j) + q
                 util2(i,j) = wtrflx(i,j)*scp2(i,j)
              endif
            enddo !i
          enddo !j
          if     (ldebug_empbal .and. &
                  mnproc.eq.1 .and. mod(nstep,100).le.1) then
            write (lp,'(i9,a,1pe12.5)')  &
             nstep,' offset wtrflx by ',q
          endif !master
        endif !not already balanced
      elseif (.not.epmass .and. empbal.eq.1) then
!
! ---   balance wtrflx to emptgt, via a region-wide offset:
! ---   virtual salt flux case, balance sss*wtrflx
!
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              util3(i,j) = wtrflx(i,j)*saln(i,j,1,n)*scp2(i,j)
              util4(i,j) =             saln(i,j,1,n)*scp2(i,j)
            endif
          enddo !i
        enddo !j
        call xcsum(d3, util3,ipa)  !total sss*wtrflx
        call xcsum(d4, util4,ipa)  !total sss
        d2 = d3/d4                 !basin-wide weighted average
        if     (d2.ne.emptgt) then  !not balanced at emptgt
          q = emptgt - d2
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                wtrflx(i,j) = wtrflx(i,j) + q
                 util2(i,j) = wtrflx(i,j)*scp2(i,j)
              endif
            enddo !i
          enddo !j
          if     (ldebug_empbal .and. &
                  mnproc.eq.1 .and. mod(nstep,100).le.1) then
            write (lp,'(i9,a,1pe12.5)')  &
             nstep,' offset wtrflx by ',q
          endif !master
        endif !not already balanced
      elseif (epmass .and. empbal.eq.2) then
!
! ---   balance wtrflx to emptgt, by scaling down +ve or -ve anomalies
!
        call xcsum(d2, util2,ipa)  !total
        d2 = d2/area               !basin-wide average
        if     (d2.ne.emptgt) then  !not balanced at emptgt
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                util3(i,j) = max(wtrflx(i,j)-emptgt,0.0)*scp2(i,j)
              endif
            enddo !i
          enddo !j
          call xcsum(d3, util3,ipa)  !positive anomaly only, >emptgt
          d3 = d3/area               !basin-wide average positive anomaly 
          d1 = d2-emptgt - d3        !basin-wide average negative anomaly
          if     (-d1.gt.d3) then
! ---       scale down the negative values
            q = -d3/d1  !<1
!$OMP       PARALLEL DO PRIVATE(j,i) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if     (wtrflx(i,j).lt.emptgt) then
                    wtrflx(i,j) = emptgt + q*(wtrflx(i,j)-emptgt)
                     util2(i,j) = wtrflx(i,j)*scp2(i,j)
                  endif
                endif
              enddo !i
            enddo !j
            if     (ldebug_empbal .and. &
                    mnproc.eq.1 .and. mod(nstep,100).le.1) then
              write (lp,'(i9,a,f8.5)')  &
               nstep,' scale -ve wtrflx anomaly by ',q
            endif !master
          else !-d1.lt.d3
! ---       scale down the positive values
            q = -d1/d3  !<1
!$OMP       PARALLEL DO PRIVATE(j,i) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if     (wtrflx(i,j).gt.emptgt) then
                    wtrflx(i,j) = emptgt + q*(wtrflx(i,j)-emptgt)
                     util2(i,j) = wtrflx(i,j)*scp2(i,j)
                  endif
                endif
              enddo !i
            enddo !j
            if     (ldebug_empbal .and. &
                    mnproc.eq.1 .and. mod(nstep,100).le.1) then
              write (lp,'(i9,a,f8.5)')  &
               nstep,' scale +ve wtrflx anomaly by ',q
            endif !master
          endif !reduce -ve:reduce +ve
        endif !not already balanced
      elseif (.not.epmass .and. empbal.eq.2) then
!
! ---   balance wtrflx to emptgt, by scaling down +ve or -ve anomalies
! ---   virtual salt flux case, balance sss*wtrflx
!
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              util2(i,j) =     wtrflx(i,j)*saln(i,j,1,n)*scp2(i,j)
              util3(i,j) = max(wtrflx(i,j)-emptgt,0.0)* &
                                           saln(i,j,1,n)*scp2(i,j)
              util4(i,j) =                 saln(i,j,1,n)*scp2(i,j)
            endif
          enddo !i
        enddo !j
        call xcsum(d2, util2,ipa)  !total sss*wtrflx
        call xcsum(d3, util3,ipa)  !anom. sss*wtrflx, wtrflx >emptgt
        call xcsum(d4, util4,ipa)  !total sss
        d2 = d2/d4                 !basin-wide weighted average
        d3 = d3/d4                 !basin-wide weighted average positive anomaly
        d1 = d2-emptgt - d3        !basin-wide weighted average negative anomaly
        if     (-d1.gt.d3) then
! ---     scale down the negative values
          q = -d3/d1  !<1
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                if     (wtrflx(i,j).lt.emptgt) then
                  wtrflx(i,j) = emptgt + q*(wtrflx(i,j)-emptgt)
                endif
                util2(i,j) = wtrflx(i,j)*scp2(i,j)
              endif
            enddo !i
          enddo !j
          if     (ldebug_empbal .and. &
                  mnproc.eq.1 .and. mod(nstep,100).le.1) then
            write (lp,'(i9,a,f8.5)')  &
             nstep,' scale -ve wtrflx anomaly by ',q
          endif !master
        else !-d1.lt.d3
! ---     scale down the positive values
          q = -d1/d3  !<1
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                if     (wtrflx(i,j).gt.emptgt) then
                  wtrflx(i,j) = emptgt + q*(wtrflx(i,j)-emptgt)
                endif
                util2(i,j) = wtrflx(i,j)*scp2(i,j)
              endif
            enddo !i
          enddo !j
          if     (ldebug_empbal .and. &
                  mnproc.eq.1 .and. mod(nstep,100).le.1) then
            write (lp,'(i9,a,f8.5)')  &
             nstep,' scale +ve wtrflx anomaly by ',q
          endif !master
        endif !reduce -ve:reduce +ve
      endif !empbal

!!Alex New calculation of epmass, with E-P applied to the top layer
!!AJW  Modified new epmass calculation to use actual dp rather than h
!!AJW  updates pbavg, dp and S.1
!
      if     (epmass) then  !requires btrlfr=.true., see blkdat.F
!$OMP   PARALLEL DO PRIVATE(j,i,k,emnp,dpemnp,dplay1,onetanew, &
!$OMP                       onetaold,pbanew,pbaold,q)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
! ---         "dp" is dp', use actual dp = dp'/oneta, this is exact for
! ---         btrmas=true and approximately correct for btrmas=false
! ---         This only works if pbavg and dp have just been integrated from
! ---         the same time, hence btrlfr=false is not allowed
!
              emnp     = -wtrflx(i,j)*svref  !m/s
              pbaold   = pbavg(i,j,n)
              pbanew   = pbavg(i,j,n) - delt1*emnp*onem  !emnp mass = rho0*vol
              pbanew   = max(pbanew, -pbot(i,j))
              onetaold = max( oneta0, 1.0 + pbaold/pbot(i,j) )
              onetanew = max( oneta0, 1.0 + pbanew/pbot(i,j) )
!
!             if     (i.eq.itest.and.j.eq.jtest) then
!               s1(1) = dp(i,j,1,n)*onetaold*saln(i,j,1,n)
!               ssum(1) = 0.0d0
!               do k= 1,kk
!                 ssum(1) = ssum(1) + dp(i,j,k,n)*onetaold*saln(i,j,k,n)
!               enddo !k
!             endif !test
!
              pbavg(i,j,n) = pbanew
              oneta(i,j,n) = max(oneta0, 1.0 + pbavg(i,j,n)/pbot(i,j))
              if     (delt1.ne.baclin) then
! ---           Robert-Asselin time filter correction
                pbavg(i,j,m) = pbavg(i,j,m)+0.5*ra2fac*(pbanew-pbaold)
                oneta(i,j,m) = max(oneta0,1.0 + pbavg(i,j,m)/pbot(i,j))
              endif
! ---         treat E-P as a layer above layer 1, merged into layer 1
! ---         E-P salinity is 0 psu, E-P temperature is SST (T.1 is unchanged)
                          q = onetaold/onetanew
                     dplay1 = dp(i,j,1,n)*q
                     dpemnp = (pbanew - pbaold)/onetanew
                     salt1o = saln(i,j,1,n)*onetaold*dp(i,j,1,n)*qonem  !diagnostic
              saln(i,j,1,n) = saln(i,j,1,n) + &
                               (0.0-saln(i,j,1,n))* &
                               dpemnp/(dplay1 + dpemnp)
                dp(i,j,1,n) = dplay1 + dpemnp
                 p(i,j,2  ) = dp(i,j,1,n)
                     salt1n = saln(i,j,1,n)*onetanew*dp(i,j,1,n)*qonem  !diagnostic
              do k= 2,kk
                dp(i,j,k,n) = dp(i,j,k,n)*q
                 p(i,j,k+1) = dp(i,j,k,n) + p(i,j,k)
              enddo !k
!
!             if     (i.eq.itest.and.j.eq.jtest) then
!                 s1(2) = dp(i,j,1,n)*onetanew*saln(i,j,1,n)
!               ssum(2) = 0.0d0
!               do k= 1,kk
!                 ssum(2) = ssum(2) + dp(i,j,k,n)*onetanew*saln(i,j,k,n)
!               enddo !k
! ---           normalize by onem, layer 1 is often 1 m thick
!               s1(1) = s1(1)/onem
!               s1(2) = s1(2)/onem
!               write(lp,'(a,i9,1p3e16.8)')
!    &            'epmass1:',nstep,s1(1),s1(2),s1(2)-s1(1)
! ---           normalize by initial depth
!               ssum(1) = ssum(1)/pbot(i,j)
!               ssum(2) = ssum(2)/pbot(i,j)
!               write(lp,'(a,i9,1p3e16.8)')
!    &           'epmassd:',nstep,ssum(1),ssum(2),ssum(2)-ssum(1)
!               write(lp,'(i9,a,3f12.6)')
!    .            nstep,',oneta   =',onetaold,onetanew,onetaold/onetanew
!               call flush(lp)
!             endif !test
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
      endif !epmass
!
! --- regiona-wide statistics
!
      call xcsum(d1, util1,ipa)
      call xcsum(d2, util2,ipa)
      watcum=watcum+d1
      empcum=empcum+d2
!
!diag if     (itest.gt.0 .and. jtest.gt.0) then
!diag   write(lp,'(i9,2i5,a/19x,4f10.4)') &
!diag     nstep,i0+itest,j0+jtest, &
!diag     '    sswflx    surflx    sstflx    wtrflx', &
!diag     sswflx(itest,jtest), &
!diag     surflx(itest,jtest), &
!diag     sstflx(itest,jtest), &
!diag     wtrflx(itest,jtest)
!diag endif !test
      return
      end subroutine thermf_oi

      subroutine thermf(m,n, dtime)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n
      real*8  dtime
!
! --- ---------------
! --- thermal forcing
! --- note: on exit flux is for ocean fraction of each grid cell
! --- ---------------
!
      logical, parameter ::  ldebug_sssbal=.false.   !usually .false.
!
      integer i,j,k,ktr,nm,l, iyear,iday,ihour
      real    day365,pwl,q,utotij,vtotij
      real*8  t1mean,s1mean,tmean,smean,pmean,rmean, &
              rareac,runsec,secpyr
      real*8  d1,d2,d3,d4
!
      real    pwij(kk+1),trwij(kk,ntracr), &
              prij(kk+1),trcij(kk,ntracr)
!
      real*8  tmean0,smean0,rmean0
      save    tmean0,smean0,rmean0
!
      double precision dtime_diurnl
      save             dtime_diurnl
      data             dtime_diurnl / -99.d0 /
!
# include "stmt_fns.h"
!
!$OMP PARALLEL DO PRIVATE(j,k,i,ktr) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do k=1,kk
          do i=1,ii
            if (SEA_P) then
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            endif !ip
          enddo !i
        enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
! --- ----------------------------
! --- thermal forcing at nestwalls
! --- ----------------------------
!
      if (nestfq.ne.0.0 .and. delt1.ne.baclin) then  !not on very 1st time step
!
!$OMP PARALLEL DO PRIVATE(j,i,k,pwl,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (ip(i,j).eq.1 .and. rmunp(i,j).ne.0.0) then
! ---       Newtonian relaxation with implict time step,
! ---         result is positive if source and nest are positive
! ---       Added by Remy Baraille, SHOM
            k=1
            saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmunp(i,j)* &
                 (snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1) )/ &
                          (1.0+delt1*rmunp(i,j))
            temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmunp(i,j)* &
                 (tnest(i,j,k,ln0)*wn0+tnest(i,j,k,ln1)*wn1) )/ &
                          (1.0+delt1*rmunp(i,j))
            th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
!
            if     (hybrid) then
              do k=kk,2,-1
                pwl=pnest(i,j,k,ln0)*wn0+pnest(i,j,k,ln1)*wn1
                if     (pwl.gt.p(i,j,kk+1)-tencm) then
                  pwl=p(i,j,kk+1)
                endif
                p(i,j,k)=min(p(i,j,k+1), &
                             ( p(i,j,k)+delt1*rmunp(i,j)*pwl )/ &
                             (1.0+delt1*rmunp(i,j)))
                dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
!
                if     (pwl.lt.p(i,j,kk+1)) then
                  saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmunp(i,j)* &
                       (snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1) )/ &
                                (1.0+delt1*rmunp(i,j))
                  if     (k.le.nhybrd) then
                    temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmunp(i,j)* &
                         (tnest(i,j,k,ln0)*wn0+tnest(i,j,k,ln1)*wn1) )/ &
                                  (1.0+delt1*rmunp(i,j))
                    th3d(i,j,k,n)=sig(temp(i,j,k,n), &
                                      saln(i,j,k,n))-thbase
                  else
                    th3d(i,j,k,n)=       theta(i,j,k)
                    temp(i,j,k,n)=tofsig(theta(i,j,k)+thbase, &
                                         saln(i,j,k,n))
                  endif
                endif
              enddo  !k
              dp(i,j,1,n)=p(i,j,2)-p(i,j,1)
            else  ! isopyc
              do k=kk,2,-1
                saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmunp(i,j)* &
                     (snest(i,j,k,ln0)*wn0+snest(i,j,k,ln1)*wn1) )/ &
                              (1.0+delt1*rmunp(i,j))
                temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
                if (k.ge.3) then
                  pwl=pnest(i,j,k,ln0)*wn0+pnest(i,j,k,ln1)*wn1
                  pwl=max(p(i,j,2),pwl)
                  if     (pwl.gt.p(i,j,kk+1)-tencm) then
                    pwl=p(i,j,kk+1)
                  endif
                  p(i,j,k)=min(p(i,j,k+1), &
                               ( p(i,j,k)+delt1*rmunp(i,j)*pwl )/ &
                               (1.0+delt1*rmunp(i,j)))
                endif
                dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
              enddo  !k
            endif  ! hybrid:isopyc
!
! ---       minimal tracer support (non-negative in buffer zone).
            do ktr= 1,ntracr
              tracer(i,j,k,n,ktr)=max(tracer(i,j,k,n,ktr),0.0)
            enddo
          endif  !ip.eq.1 .and. rmunp.ne.0.0
!
          if (iu(i,j).eq.1) then
            q  =rmunvu(i,j)
            if     (q.ne.0.0) then
              do k= 1,kk
                pwl=u(i,j,k,n)
                u(i,j,k,n)=( u(i,j,k,n)+delt1*q* &
                   (unest(i,j,k,ln0)*wn0+unest(i,j,k,ln1)*wn1) )/ &
                           (1.0+delt1*q)
              enddo  !k
            endif !rmunvu.ne.0.0
          endif  !iu.eq.1
!
          if (iv(i,j).eq.1) then
            q  =rmunvv(i,j)
            if     (q.ne.0.0) then
              do k= 1,kk
                pwl=v(i,j,k,n)
                v(i,j,k,n)=( v(i,j,k,n)+delt1*q* &
                   (vnest(i,j,k,ln0)*wn0+vnest(i,j,k,ln1)*wn1) )/ &
                           (1.0+delt1*q)
              enddo  !k
            endif  !rmunvv.ne.0.0
          endif  !iv.eq.1
        enddo  !i
      enddo  !j
!$OMP END PARALLEL DO
!
      endif  !  nestfq.ne.0.0
!
! --- ----------------------------
! --- thermal forcing at sidewalls
! --- ----------------------------
!
      if (relax .and. delt1.ne.baclin) then  !not on very 1st time step
!
!$OMP PARALLEL DO PRIVATE(j,i,k,pwl) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
      do i=1,ii
      if (SEA_P) then
        if (rmu(i,j).ne.0.0) then
! ---     Newtonian relaxation with implict time step,
! ---       result is positive if source and wall are positive
          k=1
          saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmu(i,j)* &
             ( swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1 &
              +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3) )/ &
                        (1.0+delt1*rmu(i,j))
          if     (lwflag.eq.2 .or. sstflg.gt.2   .or. &
                  icmflg.eq.2 .or. ticegr.eq.0.0     ) then
! ---       use seatmp, since it is the best available SST
            if(cpl_seatmp) then
               temp(i,j,k,n)=temp(i,j,k,n)+delt1*rmu(i,j)* &
               imp_seatmp(i,j,1)/(1.0+delt1*rmu(i,j))
            elseif (natm.eq.2) then
               temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmu(i,j)* &
                ( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1) )/ &
                                  (1.0+delt1*rmu(i,j))
            else
               temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmu(i,j)* &
                ( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1 &
                 +seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3) )/ &
                                 (1.0+delt1*rmu(i,j))
            endif
          else
            temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmu(i,j)* &
               ( twall(i,j,k,lc0)*wc0+twall(i,j,k,lc1)*wc1 &
                +twall(i,j,k,lc2)*wc2+twall(i,j,k,lc3)*wc3) )/ &
                          (1.0+delt1*rmu(i,j))
          endif
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
!
          if     (hybrid) then
            do k=kk,2,-1
              pwl=pwall(i,j,k,lc0)*wc0+pwall(i,j,k,lc1)*wc1 &
                 +pwall(i,j,k,lc2)*wc2+pwall(i,j,k,lc3)*wc3
              if     (pwl.gt.p(i,j,kk+1)-tencm) then
                pwl=p(i,j,kk+1)
              endif
              p(i,j,k)=min( p(i,j,k+1), &
                            ( p(i,j,k)+delt1*rmu(i,j)*pwl )/ &
                            (1.0+delt1*rmu(i,j)) )
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
!
              if     (pwl.lt.p(i,j,kk+1)) then
                saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmu(i,j)* &
                   ( swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1 &
                    +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3) )/ &
                              (1.0+delt1*rmu(i,j))
                if     (k.le.nhybrd) then
                  temp(i,j,k,n)=( temp(i,j,k,n)+delt1*rmu(i,j)* &
                     ( twall(i,j,k,lc0)*wc0+twall(i,j,k,lc1)*wc1 &
                      +twall(i,j,k,lc2)*wc2+twall(i,j,k,lc3)*wc3) )/ &
                                (1.0+delt1*rmu(i,j))
                  th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                else
                  th3d(i,j,k,n)=       theta(i,j,k)
                  temp(i,j,k,n)=tofsig(theta(i,j,k)+thbase, &
                                       saln(i,j,k,n))
                endif !hybrid:else
              endif !pwl.lt.p(i,j,kk+1)
            enddo !k
            dp(i,j,1,n)=p(i,j,2)-p(i,j,1)
          else  ! isopyc
            do k=kk,2,-1
              saln(i,j,k,n)=( saln(i,j,k,n)+delt1*rmu(i,j)* &
                 ( swall(i,j,k,lc0)*wc0+swall(i,j,k,lc1)*wc1 &
                  +swall(i,j,k,lc2)*wc2+swall(i,j,k,lc3)*wc3) )/ &
                            (1.0+delt1*rmu(i,j))
              temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
              if (k.ge.3) then
                pwl=pwall(i,j,k,lc0)*wc0+pwall(i,j,k,lc1)*wc1 &
                   +pwall(i,j,k,lc2)*wc2+pwall(i,j,k,lc3)*wc3
                pwl=max(p(i,j,2),pwl)
                if     (pwl.gt.p(i,j,kk+1)-tencm) then
                  pwl=p(i,j,kk+1)
                endif
                p(i,j,k)=min(p(i,j,k+1), &
                             p(i,j,k)+delt1*rmu(i,j)*(pwl-p(i,j,k)))
              endif !k.ge.3
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
            enddo !k
          endif !hybrid:isopyc
        endif !rmu(i,j).ne.0.0
      endif !ip
      enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      endif  !  relax = .true.
!
! --- ----------------------------
! --- tracer forcing at sidewalls
! --- ----------------------------
!
      if (trcrlx .and. delt1.ne.baclin) then  !not on very 1st time step
!
!$OMP   PARALLEL DO PRIVATE(j,i,k,ktr,pwij,trwij,prij,trcij) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              if     (rmutra(i,j).ne.0.0) then !at least one mask is non-zero
                prij(1)=0.0
                do k=1,kk
                  prij(k+1) =  prij(k)+dp(i,j,k,n)
                  pwij(k)   =  pwall(i,j,k,lc0)*wc0 &
                              +pwall(i,j,k,lc1)*wc1 &
                              +pwall(i,j,k,lc2)*wc2 &
                              +pwall(i,j,k,lc3)*wc3
                  do ktr= 1,ntracr
                    trwij(k,ktr) =  trwall(i,j,k,lc0,ktr)*wc0 &
                                   +trwall(i,j,k,lc1,ktr)*wc1 &
                                   +trwall(i,j,k,lc2,ktr)*wc2 &
                                   +trwall(i,j,k,lc3,ktr)*wc3
                  enddo !ktr
                enddo !k
                pwij(kk+1)=prij(kk+1)
!               call plctrc(trwij,pwij,kk,ntracr,
!    &                      trcij,prij,kk        )
                call plmtrc(trwij,pwij,kk,ntracr, &
                            trcij,prij,kk        )
                do ktr= 1,ntracr
                  if     (rmutr(i,j,ktr).ne.0.0) then
                    do k=1,kk
                      tracer(i,j,k,n,ktr) = ( tracer(i,j,k,n,ktr)+ &
                                delt1*rmutr(i,j,ktr)*trcij(k,ktr) )/ &
                                            (1.0+delt1*rmutr(i,j,ktr))
                    enddo !k
                  endif !rmutr.ktr.ne.0.0
                enddo !ktr
              endif !rmutra.ne.0.0
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
      endif  !  trcrlx = .true.
!
! --- ---------------------------------------------------------
! --- Update dpu,dpv, and rebalance velocity, if dp has changed
! --- ---------------------------------------------------------
!
      if ((nestfq.ne.0.0 .and. delt1.ne.baclin) .or. &
          (relax         .and. delt1.ne.baclin)     ) then
        call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                    dpv(1-nbdy,1-nbdy,1,n), &
                    p,depthu,depthv, 0,0)
!
!$OMP   PARALLEL DO PRIVATE(j,i,k,utotij,vtotij) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (iu(i,j).eq.1 .and. &
                max(rmunvu(i,j), &
                    rmu(   i,j),rmu(i-1,j) ).ne.0.0) then
              utotij = 0.0                                     
              do k=1,kk                                        
                utotij = utotij + u(i,j,k,n)*dpu(i,j,k,n)
              enddo ! k
              utotij=utotij/depthu(i,j)
              do k=1,kk
                u(i,j,k,n) = u(i,j,k,n) - utotij
              enddo ! k
            endif  !rebalance u
!
            if (iv(i,j).eq.1 .and. &
                max(rmunvv(i,j), &
                    rmu(   i,j),rmu(i,j-1) ).ne.0.0) then
              vtotij = 0.0
              do k=1,kk
                vtotij = vtotij + v(i,j,k,n)*dpv(i,j,k,n)
              enddo ! k
              vtotij=vtotij/depthv(i,j)
              do k=1,kk
                v(i,j,k,n) = v(i,j,k,n) - vtotij
              enddo ! k
            endif  !rebalance v
          enddo  !i
        enddo  !j
!$OMP   END PARALLEL DO
      endif !update dpu,dpv and rebalance u,v
!
! --- --------------------------------
! --- thermal forcing of ocean surface
! --- --------------------------------
!
!
      if (thermo .or. sstflg.gt.0 .or. srelax) then
!
      if     (dswflg.eq.1 .and. dtime-dtime_diurnl.gt.1.0) then
! ---   update diurnal factor table
        call forday(dtime,yrflag, iyear,iday,ihour)
        day365 = mod(iday+364,365)
        call thermf_diurnal(diurnl, day365)
        dtime_diurnl = dtime
!diag       if     (mnproc.eq.1) then
!diag       write (lp,'(a)') 'diurnl updated'
!diag       endif !1st tile
      endif
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call thermfj(m,n,dtime, j)
      enddo
!$OMP END PARALLEL DO
!
!diag if     (itest.gt.0 .and. jtest.gt.0) then
!diag   write(lp,'(i9,2i5,a/19x,4f10.4)') &
!diag     nstep,i0+itest,j0+jtest, &
!diag     '    sstflx     ustar    hekman    surflx', &
!diag     sstflx(itest,jtest), &
!diag     ustar( itest,jtest), &
!diag     hekman(itest,jtest), &
!diag     surflx(itest,jtest)
!diag   write(lp,'(i9,2i5,a/19x,4f10.4)') &
!diag     nstep,i0+itest,j0+jtest, &
!diag     '    sswflx     wtrflx   rivflx    sssflx', &
!diag     sswflx(itest,jtest), &
!diag     wtrflx(itest,jtest), &
!diag     rivflx(itest,jtest), &
!diag     sssflx(itest,jtest)
!diag endif !test
!
! --- smooth surface fluxes?
!
      if     (flxsmo) then
        call psmooth_ice(surflx, 0,0, ishlf, util1)  !uses covice
        call psmooth_ice(wtrflx, 0,0, ishlf, util1)  !uses covice
      endif
!
      if     (sssbal.eq.1) then
!
! ---   balance sssflx, via a region-wide offset:
! ---     util2 from thermfj holds sssflx*scp2
!
        call xcsum(d2, util2,ipa)  !total
        if     (d2.ne.0.0d0) then  !not balanced
          q = -d2/area
          do j=1,jj
            do i=1,ii
              if (SEA_P) then
                sssflx(i,j) = sssflx(i,j) + q
              endif
            enddo !i
          enddo !j
          if     (ldebug_sssbal .and. &
                  mnproc.eq.1 .and. mod(nstep,100).le.1) then
            write (lp,'(i9,a,1pe12.5)')  &
             nstep,' offset sssflx by ',q
          endif !master
        endif !not already balanced
      elseif (sssbal.eq.2) then
!
! ---   balance sssflx, by scaling up relaxation:
! ---     util2 from thermfj holds     sssflx*scp2
! ---     util1 from thermfj holds max(sssflx*scp2), i.e. positive only
!
        call xcsum(d2, util2,ipa)  !total
        if     (d2.ne.0.0d0) then  !not balanced
          call xcsum(d1, util1,ipa)  !positive only
          d3 = d2 - d1               !negative only
          if     (d1.gt.-d3) then
! ---       scale up the negative values
            q = min(-d1/d3, 5.0 )  ! >1
!$OMP       PARALLEL DO PRIVATE(j,k,l,i) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if     (sssflx(i,j).lt.0.0) then
                    sssflx(i,j) = q*sssflx(i,j)
                  endif
                endif
              enddo !i
            enddo !j
            if     (ldebug_sssbal .and. &
                    mnproc.eq.1 .and. mod(nstep,100).le.1) then
              write (lp,'(i9,a,f6.3)')  &
               nstep,' scale -ve sssflx by ',q
            endif !master
          else !d1.lt.-d3
! ---       scale up the positive values
            q = min(-d3/d1, 5.0 )  ! >1
!$OMP       PARALLEL DO PRIVATE(j,k,l,i) &
!$OMP                SCHEDULE(STATIC,jblk)
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if     (sssflx(i,j).gt.0.0) then
                    sssflx(i,j) = q*sssflx(i,j)
                  endif
                endif
              enddo !i
            enddo !j
            if     (ldebug_sssbal .and. &
                    mnproc.eq.1 .and. mod(nstep,100).le.1) then
              write (lp,'(i9,a,f6.3)')  &
               nstep,' scale +ve sssflx by ',q
            endif !master
          endif !reduce +ve:reduce -ve
        endif !not already balanced
      endif !sssbal
!###
!
      if (nstep.eq.nstep1+1 .or. diagno) then
        if (nstep.eq.nstep1+1) then
          nm=m
        else
          nm=n
        endif
!$OMP   PARALLEL DO PRIVATE(j,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              util1(i,j)=temp(i,j,1,nm)*scp2(i,j)
              util2(i,j)=saln(i,j,1,nm)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
        call xcsum(d1, util1,ipa)
        call xcsum(d2, util2,ipa)
        t1mean=d1
        s1mean=d2
!
!$OMP   PARALLEL DO PRIVATE(j,k,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          k=1
            do i=1,ii
              if (SEA_P) then
! ---           always use mass-conserving diagnostics
                oneta(i,j,nm)= max(oneta0,1.0 + pbavg(i,j,nm)/pbot(i,j))
                q=oneta(i,j,nm)*dp(i,j,k,nm)*scp2(i,j)
                util1(i,j)=q
                util2(i,j)=q*temp(i,j,k,nm)
                util3(i,j)=q*saln(i,j,k,nm)
                util4(i,j)=q*th3d(i,j,k,nm)
              endif !ip
            enddo !i
          do k=2,kk
            do i=1,ii
              if (SEA_P) then
                q=oneta(i,j,nm)*dp(i,j,k,nm)*scp2(i,j)
                util1(i,j)=util1(i,j)+q
                util2(i,j)=util2(i,j)+q*temp(i,j,k,nm)
                util3(i,j)=util3(i,j)+q*saln(i,j,k,nm)
                util4(i,j)=util4(i,j)+q*th3d(i,j,k,nm)
              endif !ip
            enddo !i
          enddo !k
        enddo !j
!$OMP   END PARALLEL DO
        call xcsum(d1, util1,ipa)
        call xcsum(d2, util2,ipa)
        call xcsum(d3, util3,ipa)
        call xcsum(d4, util4,ipa)
        pmean=d1
        tmean=d2/pmean
        smean=d3/pmean
        rmean=d4/pmean
        if     (mnproc.eq.1) then
        write (lp,'(i9,a,3f10.3)')  &
          nstep,' mean basin temp, saln, dens ', &
          tmean,smean,rmean+thbase
        endif !1st tile
        if     (nstep.eq.nstep1+1) then
!
! ---     save initial basin means.
          tmean0=tmean
          smean0=smean
          rmean0=rmean
        else
!
! ---     diagnostic printout of fluxes.
          rareac=1.0/(area*(nstep-nstep1))
          runsec=   baclin*(nstep-nstep1)
          if      (yrflag.eq.0) then
            secpyr=360.00d0*86400.0d0
          elseif (yrflag.lt.3) then
            secpyr=366.00d0*86400.0d0
          elseif (yrflag.eq.3) then
            secpyr=365.25d0*86400.0d0
          elseif (yrflag.eq.4) then
            secpyr=365.00d0*86400.0d0
          endif
          if     (mnproc.eq.1) then
          write (lp,'(i9,a,2f10.3)')  &
           nstep,' mean surface temp and saln  ', &
           t1mean/area,s1mean/area
          write (lp,'(i9,a,2f10.3,a)')  &
           nstep,' energy residual (atmos,tot) ', &
           watcum*rareac, &
           (tmean-tmean0)*(spcifh*avgbot*rhoref)/runsec, &
          ' (W/m^2)'
          write (lp,'(i9,a,2f10.3,a)') &
           nstep,'  e - p residual (atmos,tot) ', &
           empcum*svref*rareac*100.0*secpyr, &
           (smean-smean0)/(saln0*runsec)*avgbot*100.0*secpyr, &
          ' (cm/year)'
          write (lp,'(i9,a,2f10.3)')  &
           nstep,' temp drift per century      ', &
           (watcum*rareac/(spcifh*avgbot*rhoref))*(secpyr*100.0d0), &
           (tmean-tmean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,f10.3)')  &
           nstep,' saln drift per century      ', &
           (smean-smean0)*(secpyr*100.0d0)/runsec
          write (lp,'(i9,a,9x,f10.3)')  &
           nstep,' dens drift per century      ', &
           (rmean-rmean0)*(secpyr*100.0d0)/runsec
          endif !1st tile
          call xcsync(flush_lp)
        endif !master
      endif !diagno
!c
      endif   !  thermo .or.  sstflg.gt.0 .or. srelax
!
      return
      end subroutine thermf
!
      subroutine thermfj(m,n,dtime, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
#if defined(STOKES)
      use mod_stokes  !    HYCOM Stokes Drift
#endif
      implicit none
!
      integer m,n, j
      real*8  dtime
!
! --- thermal forcing of ocean surface, for row j.
!
      integer i,ihr,it_a,ilat
      real    radfl,swfl,swflc,sstrlx,wind,airt,vpmx,prcp,xtau,ytau, &
              evap,evape,emnp,esst,sssf, &
              snsibl,dsgdt,sssc,sssdif,sstdif,rmut
      real    cd0,clh,cl0,cl1,csh, &
              pair,rair,slat,ssen,tdif,tsur,wsph, &
              tamts,q,qva,va
      real    swscl,xhr,xlat
      real    u10,v10,uw10,uw
      real    cd_n10,cd_n10_rt,ce_n10,ch_n10,cd_rt,stab, &
              tv,tstar,qstar,bstar,zeta,x2,x,xx, &
              psi_m,psi_h,z0,qrair,zi
      real    cd10,ce10,ch10,ustar1
      real*8  dloc
      integer k
!
! --- 'ustrmn' = minimum ustar
! --- 'cormn4' = 4 times minimum coriolis magnitude
! --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
! --- 'evaplh' = latent heat of evaporation (j/kg)
! --- 'csice'  = ice-air sensible exchange coefficient
!
      real       ustrmn,cormn4,csubp,evaplh,csice
      parameter (ustrmn=1.0e-5,  &
                 cormn4=4.0e-5,   & ! corio(4N) is about 1.e-5
                 csubp =1005.7, &
                 evaplh=2.47e6, &
                 csice =0.0006)
!
! --- parameter for lwflag=-1
! --- 'sb_cst' = Stefan-Boltzman constant
      real       sb_cst
      parameter (sb_cst=5.67e-8)
!
! --- parameters primarily for flxflg=1 (ustflg=1)
! --- 'airdns' = air density at sea level (kg/m**3)
! --- 'cd'     = drag coefficient
! --- 'ctl'    = thermal transfer coefficient (latent)
! --- 'cts1'   = thermal transfer coefficient (sensible, stable)
! --- 'cts2'   = thermal transfer coefficient (sensible, unstable)
!
      real       airdns,cd,ctl,cts1,cts2
      parameter (airdns=1.2)
      parameter (cd  =0.0013, ctl =0.0012, &
                 cts1=0.0012, cts2=0.0012)
!
! --- parameters primarily for flxflg=2 (ustflg=2)
! --- 'pairc'  = air pressure (Pa) or (mb) * 100
! --- 'rgas'   = gas constant (j/kg/k)
! --- 'tzero'  = celsius to kelvin temperature offset
! --- 'clmin'  = minimum allowed cl
! --- 'clmax'  = maximum allowed cl
! --- 'wsmin'  = minimum allowed wind speed (for cl and cd)
! --- 'wsmax'  = maximum allowed wind speed (for cl and cd)
!
      real       pairc,rgas,tzero,clmin,clmax,wsmin,wsmax
      parameter (pairc=1013.0*100.0, &
                 rgas =287.1,   tzero=273.16,  &
                 clmin=0.0003,  clmax=0.002, &
                 wsmin=3.5,     wsmax=27.5)
!
! --- parameters primarily for flxflg=4
! --- 'lvtc'   = include a virtual temperature correction
! --- 'vamin'  = minimum allowed wind speed (for cl)
! --- 'vamax'  = maximum allowed wind speed (for cl)
! --- 'tdmin'  = minimum allowed Ta-Ts      (for cl)
! --- 'tdmax'  = maximum allowed Ta-Ts      (for cl)
!
! --- 'as0_??' =  stable Ta-Ts  polynominal coefficients, va<=5m/s
! --- 'as5_??' =  stable Ta-Ts  polynominal coefficients, va>=5m/s
! --- 'au0_??' = unstable Ta-Ts polynominal coefficients, va<=5m/s
! --- 'au5_??' = unstable Ta-Ts polynominal coefficients, va>=5m/s
! --- 'an0_??' =  neutral Ta-Ts polynominal coefficients, va<=5m/s
! --- 'an5_??' =  neutral Ta-Ts polynominal coefficients, va>=5m/s
! --- 'ap0_??' =    +0.75 Ta-Ts polynominal coefficients, va<=5m/s
! --- 'ap5_??' =    +0.75 Ta-Ts polynominal coefficients, va>=5m/s
! --- 'am0_??' =    -0.75 Ta-Ts polynominal coefficients, va<=5m/s
! --- 'am5_??' =    -0.75 Ta-Ts polynominal coefficients, va>=5m/s
!
      logical, parameter :: lvtc =.true.
      real,    parameter :: vamin= 1.2, vamax=34.0
      real,    parameter :: tdmin=-8.0, tdmax= 7.0
!
      real, parameter :: &
        as0_00=-2.925e-4,   as0_10= 7.272e-5,  as0_20=-6.948e-6, &
        as0_01= 5.498e-4,   as0_11=-1.740e-4,  as0_21= 1.637e-5, &
        as0_02=-5.544e-5,   as0_12= 2.489e-5,  as0_22=-2.618e-6
      real, parameter :: &
        as5_00= 1.023e-3,   as5_10=-2.672e-6,  as5_20= 1.546e-6, &
        as5_01= 9.657e-6,   as5_11= 2.103e-4,  as5_21=-6.228e-5, &
        as5_02=-2.281e-8,   as5_12=-5.329e-3,  as5_22= 5.094e-4
      real, parameter :: &
        au0_00= 2.077e-3,   au0_10=-2.899e-4,  au0_20=-1.954e-5, &
        au0_01=-3.933e-4,   au0_11= 7.350e-5,  au0_21= 5.483e-6, &
        au0_02= 3.971e-5,   au0_12=-6.267e-6,  au0_22=-4.867e-7
      real, parameter :: &
        au5_00= 1.074e-3,   au5_10= 6.912e-6,  au5_20= 1.849e-7, &
        au5_01= 5.579e-6,   au5_11=-2.244e-4,  au5_21=-2.167e-6, &
        au5_02= 5.263e-8,   au5_12=-1.027e-3,  au5_22=-1.010e-4
      real, parameter :: &
        an0_00= 1.14086e-3, an5_00= 1.073e-3, &
        an0_01=-3.120e-6,   an5_01= 5.531e-6, &
        an0_02=-9.300e-7,   an5_02= 5.433e-8
      real, parameter :: &
        ap0_00= as0_00 + as0_10*0.75 + as0_20*0.75**2, &
        ap0_01= as0_01 + as0_11*0.75 + as0_21*0.75**2, &
        ap0_02= as0_02 + as0_12*0.75 + as0_22*0.75**2
      real, parameter :: &
        ap5_00= as5_00 + as5_10*0.75 + as5_20*0.75**2, &
        ap5_01= as5_01, &
        ap5_02= as5_02, &
        ap5_11=          as5_11*0.75 + as5_21*0.75**2, &
        ap5_12=          as5_12*0.75 + as5_22*0.75**2
      real, parameter :: &
        am0_00= au0_00 - au0_10*0.75 + au0_20*0.75**2, &
        am0_01= au0_01 - au0_11*0.75 + au0_21*0.75**2, &
        am0_02= au0_02 - au0_12*0.75 + au0_22*0.75**2
      real, parameter :: &
        am5_00= au5_00 - au5_10*0.75 + au5_20*0.75**2, &
        am5_01= au5_01, &
        am5_02= au5_02, &
        am5_11=        - au5_11*0.75 + au5_21*0.75**2, &
        am5_12=        - au5_12*0.75 + au5_22*0.75**2
!
! --- parameters primarily for flxflg=4
      real, parameter :: vonkar=0.4         !Von Karmann constant
      real, parameter :: cpcore=1000.5      !specific heat of air (j/kg/deg)
!
      real satvpr,qsatur6,qsatur,qsatur5,t6,p6,f6,qra !t declared in stmt_fns.h
# include "stmt_fns.h"
!
! --- saturation vapor pressure (Pa),
! --- from a polynominal approximation (lowe, j.appl.met., 16, 100-103, 1976)
      satvpr(t)=  100.0*(6.107799961e+00+t*(4.436518521e-01 &
                     +t*(1.428945805e-02+t*(2.650648471e-04 &
                     +t*(3.031240396e-06+t*(2.034080948e-08 &
                     +t* 6.136820929e-11))))))
!
! --- pressure dependent saturation mixing ratio (kg/kg)
! --- p6 is pressure in Pa, f6 is fractional depression from SSS
      qsatur6(t6,p6,f6)=0.622*(f6*satvpr(t6)/(p6-f6*satvpr(t6)))
!
! --- saturation mixing ratio (kg/kg), from a polynominal approximation
! --- for saturation vapor pressure (lowe, j.appl.met., 16, 100-103, 1976)
! --- assumes that (mslprs-satvpr(t)) is 1.e5 Pa
      qsatur(t)=.622e-3*(6.107799961e+00+t*(4.436518521e-01 &
                     +t*(1.428945805e-02+t*(2.650648471e-04 &
                     +t*(3.031240396e-06+t*(2.034080948e-08 &
                     +t* 6.136820929e-11))))))
!
! --- saturation specific humidity (flxflg=5)
! --- qra is 1/rair
      qsatur5(t,qra)= 0.98*qra*6.40380e5*exp(-5107.4/(t+tzero))
!
! --- temperature relaxation coefficient
      rmut=1./(30.0*86400.0)  !1/30 days
!
! --- ------------------------------------------------------
! --- thermal forcing of ocean surface (positive into ocean)
! --- ------------------------------------------------------
!
      do i=1,ii
      if (SEA_P) then
      if     (flxflg.gt.0) then
! ---   wind = wind, or wind-ocean, speed (m/s)
        if     (flxflg.eq.6 .and. amoflg.ne.0) then
          wind=wndocn(i,j)  !magnitude of wind minus ocean current
        elseif(cpl_wndspd) then
          wind=imp_wndspd(i,j,1)
        elseif (natm.eq.2) then
          wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
        else
          wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1 &
              +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
        endif !natm
! ---   swfl = shortwave radiative thermal flux (W/m^2) +ve into ocean/ice
! ---          Qsw includes the atmos. model's surface albedo,
! ---          i.e. it already allows for sea-ice&snow where it is observed.
        if(cpl_swflx) then
          swfl=imp_swflx (i,j,1)
        elseif (natm.eq.2) then
          swfl=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
        else
          swfl=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1 &
              +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
        endif !natm
        if     (dswflg.eq.1) then
! ---     daily to diurnal shortwave correction to swfl and radfl.
          dloc  = dtime + plon(i,j)/360.0
          xhr   = (dloc - int(dloc))*24.0  !local time of day
          ihr   = int(xhr)
          xhr   =     xhr - ihr
          if     (plat(i,j).ge.0.0) then
            ilat  = int(plat(i,j))
            xlat  =     plat(i,j) - ilat
          else
            ilat  = int(plat(i,j)) - 1
            xlat  =     plat(i,j) - ilat
          endif
          swscl = (1.0-xhr)*(1.0-xlat)*diurnl(ihr,  ilat  ) + &
                  (1.0-xhr)*     xlat *diurnl(ihr,  ilat+1) + &
                       xhr *(1.0-xlat)*diurnl(ihr+1,ilat  ) + &
                       xhr *     xlat *diurnl(ihr+1,ilat+1)
          if(cpl_swflx) then
            swflc = (swscl-1.0)*imp_swflx (i,j,1)
            swfl  =  swscl     *imp_swflx (i,j,1)
          else
            swflc = (swscl-1.0)*swfl  !diurnal correction only
            swfl  =  swscl     *swfl
          endif
!diag         if     (i.eq.itest.and.j.eq.jtest) then
!diag           write(lp,'(i9,a,2i5,2f8.5)') &
!diag             nstep,', hr,lat =',ihr,ilat,xhr,xlat
!diag           write(lp,'(i9,a,5f8.5)') &
!diag             nstep,', swscl  =',swscl,diurnl(ihr,  ilat  ), &
!diag                                      diurnl(ihr,  ilat+1), &
!diag                                      diurnl(ihr+1,ilat  ), &
!diag                                      diurnl(ihr+1,ilat+1)
!diag           call flush(lp)
!diag         endif !test
        else
          swflc = 0.0 !no diurnal correction
        endif !dswflg
! ---   radfl= net       radiative thermal flux (W/m^2) +ve into ocean/ice
! ---        = Qsw+Qlw across the atmosphere to ocean or sea-ice interface
        if(cpl_swflx .and. cpl_lwmdnflx .and. cpl_lwmupflx) then
           radfl= imp_swflx (i,j,1) &
                 +imp_lwdflx(i,j,1) &
                 +imp_lwuflx(i,j,1)
        elseif (natm.eq.2) then

          radfl=(radflx(i,j,l0)*w0+radflx(i,j,l1)*w1)
        else
          radfl=(radflx(i,j,l0)*w0+radflx(i,j,l1)*w1 &
                +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3)
        endif !natm
        if     (lwflag.eq.-1) then
! ---     input radflx is Qlwdn, convert to Qlw + Qsw
          if(cpl_swflx .and. cpl_lwmdnflx .and. cpl_lwmupflx ) then
             radfl= imp_swflx (i,j,1) &
                   +imp_lwdflx(i,j,1) &
!                   - sb_cst*(temp(i,j,1,n)+tzero)**4
                   +imp_lwuflx(i,j,1)
          else
             radfl = radfl - sb_cst*(temp(i,j,1,n)+tzero)**4 + swfl
          endif

          sstflx(i,j) = 0.0
        elseif (lwflag.gt.0) then
! ---     over-ocean longwave correction to radfl (Qsw+Qlw).
          tsur = temp(i,j,1,n)
          if     (lwflag.eq.1) then !from climatology
            tdif = tsur - &
                   ( twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1 &
                    +twall(i,j,1,lc2)*wc2+twall(i,j,1,lc3)*wc3)
          else !w.r.t. atmospheric model's sst
            if(cpl_surtmp) then
                  tdif = tsur - imp_surtmp(i,j,1)
            elseif (natm.eq.2) then
              tdif = tsur - &
                     ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1)
            else
              tdif = tsur - &
                     ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1 &
                      +surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3)
            endif !natm
          endif
          !correction is blackbody radiation from tdif at tsur
          radfl = radfl - (4.506+0.0554*tsur) * tdif
          !count the correction as a relaxation term
          sstflx(i,j) = - (4.506+0.0554*tsur) * tdif
          !allow for any diurnal correction
          radfl = radfl + swflc
        else
          sstflx(i,j) = 0.0
          !allow for any diurnal correction
          radfl = radfl + swflc
        endif
!diag         if     (i.eq.itest.and.j.eq.jtest) then
!diag           write(lp,'(i9,a,4f8.2)') &
!diag             nstep,', radfl  =',radfl,swflc,swfl,swfl-swflc
!diag           call flush(lp)
!diag         endif !test
        if     (pcipf) then
! ---     prcp = precipitation (m/sec; positive into ocean)
! ---     note that if empflg==3, this is actually P-E
          if(cpl_precip) then
            prcp=imp_precip(i,j,1)
          elseif (natm.eq.2) then
            prcp=precip(i,j,l0)*w0+precip(i,j,l1)*w1
          else
            prcp=precip(i,j,l0)*w0+precip(i,j,l1)*w1 &
                +precip(i,j,l2)*w2+precip(i,j,l3)*w3
          endif !natm
        endif
        if     (empflg.lt.0) then  !observed (or NWP) SST
          if (cpl_seatmp) then
            esst = imp_seatmp(i,j,1)
          elseif (natm.eq.2) then
            esst = seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1
          else
            esst = seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1+ &
                   seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3
          endif !natm
        endif
        if     (flxflg.ne.3) then
          if(cpl_airtmp .and. cpl_vapmix) then
             airt=imp_airtmp(i,j,1)
             vpmx=imp_vapmix(i,j,1)
          elseif (natm.eq.2) then
            airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
            vpmx=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1
          else
! ---       airt = air temperature (C)
            airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1 &
                +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
! ---       vpmx = water vapor mixing ratio (kg/kg)
! ---       vpmx = specific humidity (kg/kg) when flxflg==5
            vpmx=vapmix(i,j,l0)*w0+vapmix(i,j,l1)*w1 &
                +vapmix(i,j,l2)*w2+vapmix(i,j,l3)*w3
          endif !natm
        endif
! ---   pair = msl pressure (Pa)
        if     (mslprf .or. flxflg.eq.6) then
          if     (natm.eq.2) then
            pair=mslprs(i,j,l0)*w0+mslprs(i,j,l1)*w1 &
                +prsbas
          else
            pair=mslprs(i,j,l0)*w0+mslprs(i,j,l1)*w1 &
                +mslprs(i,j,l2)*w2+mslprs(i,j,l3)*w3 &
                +prsbas
          endif !natm
        else
          pair=pairc
        endif
! ---   ustar = U* (sqrt(N.m/kg))                 
        if     (ustflg.eq.3) then !ustar from input
          if(cpl_ustara) then
            ustar(i,j)=imp_ustara(i,j,1)
          elseif (natm.eq.2) then
            ustar(i,j)=ustara(i,j,l0)*w0+ustara(i,j,l1)*w1
          else
            ustar(i,j)=ustara(i,j,l0)*w0+ustara(i,j,l1)*w1 &
                      +ustara(i,j,l2)*w2+ustara(i,j,l3)*w3
          endif !natm
        elseif (ustflg.eq.1) then !ustar from wndspd, constant cd
          ustar(i,j)=sqrt(svref*cd*airdns)*wind
        elseif (ustflg.eq.2) then !ustar from wndspd, variable cd
          wsph = min( wsmax, max( wsmin, wind ) )
          cd0  = 0.862e-3 + 0.088e-3 * wsph - 0.00089e-3 * wsph**2
          rair = pair / (rgas * ( tzero + airt ))
          ustar(i,j)=sqrt(svref*cd0*rair)*wind
        elseif (ustflg.eq.4) then !ustar from surface stress, see montum_hs
          ustar(i,j)=sqrt(svref*sqrt(surtx(i,j)**2+surty(i,j)**2))
        endif !ustflg
        ustar( i,j)=max(ustrmn,ustar(i,j))
        hekman(i,j)=ustar(i,j)*(cekman*4.0)/ &
                     max( cormn4, &
                          abs(corio(i,j  ))+abs(corio(i+1,j  ))+ &
                          abs(corio(i,j+1))+abs(corio(i+1,j+1)) )
      else !flxlfg==0, i.e. no flux
        swfl=0.0
        sstflx(i,j)=0.0
        mixflx(i,j)=0.0
        buoflx(i,j)=0.0
        bhtflx(i,j)=0.0
        ustar( i,j)=0.0
        hekman(i,j)=0.0
      endif !flxflg
!
      if     (flxflg.eq.1) then
!
! ---   MICOM bulk air-sea flux parameterization
! ---   (constant Cl and constant stable/unstable Cs)
!
        if (temp(i,j,1,n).lt.airt) then
          csh=cts1  !stable
        else
          csh=cts2  !unstable
        endif
! ---   evap   = evaporation (W/m^2) into atmos from ocean.
! ---   snsibl = sensible heat flux  into atmos from ocean.
        if     (empflg.lt.0) then
          evape = ctl*airdns*evaplh*wind* &
                  max(0.,0.97*qsatur(esst)-vpmx)
        endif
! ---   Latent Heat flux (W/m2)
        if(cpl_latflx) then
            evap=imp_latflx(i,j,1)
        else
            evap=ctl*airdns*evaplh*wind* &
                 max(0.,0.97*qsatur(temp(i,j,1,n))-vpmx)
        endif
! ---   Sensible Heat flux (W/m2)
        if(cpl_sensflx) then
            snsibl=imp_sensflx(i,j,1)
        else
            snsibl=csh*airdns*csubp*wind*(temp(i,j,1,n)-airt)
        endif
! ---   surflx = thermal energy flux (W/m^2) into ocean
        surflx(i,j)=radfl - snsibl - evap
      elseif (flxflg.eq.2) then
!
! ---    Cl (and Cs) depend on wind speed and Ta-Ts.
! ---    Kara, A. B., P. A. Rochford, and H. E. Hurlburt, 2002:
! ---    Air-sea flux estimates and the 1997-1998 ENSO event.
! ---    Bound.-Layer Meteor., 103, 439-458.
! ---    http://www7320.nrlssc.navy.mil/pubs.php
!
        rair = pair / (rgas * ( tzero + airt ))
        slat = evaplh*rair
        ssen = csubp *rair
!
        tdif = temp(i,j,1,n) - airt
        wsph = min( wsmax, max( wsmin, wind ) )
        cl0  =  0.885e-3 + 0.0748e-3 * wsph - 0.00143e-3 * wsph**2
        cl1  = -0.113e-4 + 4.89e-4   / wsph
        clh  = min( clmax, max( clmin, cl0 + cl1 * tdif ) )
        csh  = 0.9554*clh
!
! ---   evap   = evaporation         (W/m^2) into atmos from ocean.
! ---   snsibl = sensible heat flux  (W/m^2) into atmos from ocean.
! ---   surflx = thermal energy flux (W/m^2) into ocean
        if     (empflg.lt.0) then
          evape = slat*clh*wind*(0.97*qsatur(esst)-vpmx)
        endif

! ---   Latent Heat flux (W/m2)
        if(cpl_latflx) then
           evap=imp_latflx(i,j,1)
        else
           evap   = slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
        endif
! ---   Sensible Heat flux (W/m2)
        if(cpl_sensflx) then
           snsibl=imp_sensflx(i,j,1)
        else
           snsibl = ssen*csh*wind* tdif
        endif
        surflx(i,j) = radfl - snsibl - evap
!
!diag   if     (i.eq.itest.and.j.eq.jtest) then
!diag     write(lp,'(i9,2i5,a,4f8.5)') &
!diag     nstep,i0+i,j0+j,' cl0,cl,cs,cd    = ',cl0,clh,csh,cd0
!diag     write(lp,'(i9,2i5,a,2f8.2,f8.5)') &
!diag     nstep,i0+i,j0+j,' wsph,tdif,ustar = ',wsph,tdif,ustar(i,j)
!diag     call flush(lp)
!diag   endif
      elseif (flxflg.eq.4) then
!
! ---   Similar to flxflg.eq.2, but with Cl based on an approximation
! ---   to values from the COARE 3.0 algorithm (Fairall et al., 2003), 
! ---   for Cl over the global ocean in the range 1m/s <= Va <= 40m/s
! ---   and -8degC <= Ta-Ts <= 7degC, that is quadratic in Ta-Ts and
! ---   quadratic in either Va or 1/Va (Kara et al.,  2005).
!
! ---   Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B.
! ---   Edson, 2003:  Bulk parameterization of air-sea fluxes:  Updates 
! ---   and verification for the COARE algorithm.  J. Climate, 16, 571-591.
!
! ---   Kara, A. B., H. E. Hurlburt, and A. J. Wallcraft, 2005:
! ---   Stability-dependent exchange coefficients for air-sea fluxes.
! ---   J. Atmos. Oceanic. Technol., 22, 1080-1094. 
! ---   http://www7320.nrlssc.navy.mil/pubs.php
!
        rair = pair / (rgas * ( tzero + airt ))
        slat = evaplh*rair
        ssen = csubp *rair
!
        tdif  = temp(i,j,1,n) - airt
        if     (lvtc) then !correct tamts for 100% humidity
          tamts = -tdif - 0.61*(airt+tzero)*(qsatur(airt)-vpmx)
          tamts = min( tdmax, max( tdmin, tamts ) )
        else
          tamts = min( tdmax, max( tdmin, -tdif ) )
        endif !lvtc:else
        va    = min( vamax, max( vamin,  wind ) )
        if     (va.le.5.0) then
          if     (tamts.gt. 0.75) then !stable
            clh =   (as0_00 + as0_01* va + as0_02* va**2) &
                  + (as0_10 + as0_11* va + as0_12* va**2)*tamts &
                  + (as0_20 + as0_21* va + as0_22* va**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au0_00 + au0_01* va + au0_02* va**2) &
                  + (au0_10 + au0_11* va + au0_12* va**2)*tamts &
                  + (au0_20 + au0_21* va + au0_22* va**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap0_00 + ap0_01* va + ap0_02* va**2) &
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am0_00 + am0_01* va + am0_02* va**2) &
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          endif !tamts
        else !va>5
          qva = 1.0/va
          if     (tamts.gt. 0.75) then !stable
            clh =   (as5_00 + as5_01* va + as5_02* va**2) &
                  + (as5_10 + as5_11*qva + as5_12*qva**2)*tamts &
                  + (as5_20 + as5_21*qva + as5_22*qva**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au5_00 + au5_01* va + au5_02* va**2) &
                  + (au5_10 + au5_11*qva + au5_12*qva**2)*tamts &
                  + (au5_20 + au5_21*qva + au5_22*qva**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap5_00 + ap5_01* va + ap5_02* va**2 &
                                  + ap5_11*qva + ap5_12*qva**2) &
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am5_00 + am5_01* va + am5_02* va**2 &
                                  + am5_11*qva + am5_12*qva**2) &
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          endif !tamts
        endif !va
        csh  = 0.9554*clh
!
! ---   evap   = evaporation         (W/m^2) into atmos from ocean.
! ---   snsibl = sensible heat flux  (W/m^2) into atmos from ocean.
! ---   surflx = thermal energy flux (W/m^2) into ocean
        if     (empflg.lt.0) then
          evape = slat*clh*wind*(0.97*qsatur(esst)-vpmx)
        endif
! ---   Latent Heat flux (W/m2)
        if(cpl_latflx) then
           evap=imp_latflx(i,j,1)
        else
           evap   = slat*clh*wind*(0.97*qsatur(temp(i,j,1,n))-vpmx)
        endif
! ---   Sensible Heat flux (W/m2)
        if(cpl_sensflx) then
           snsibl=imp_sensflx(i,j,1)
        else
           snsibl = ssen*csh*wind* tdif
        endif
        surflx(i,j) = radfl - snsibl - evap
!
!diag   if     (i.eq.itest.and.j.eq.jtest) then
!diag     write(lp,'(i9,2i5,a,3f8.5)') &
!diag     nstep,i0+i,j0+j,' cl,cs,cd    = ',clh,csh,cd0
!diag     write(lp,'(i9,2i5,a,2f8.2,f8.5)') &
!diag     nstep,i0+i,j0+j,' va,tamst,ustar = ',va,tamts,ustar(i,j)
!diag     call flush(lp)
!diag   endif
      elseif (flxflg.eq.6) then
!
! ---   Similar to flxflg.eq.4, but with more pressure dependance,
! ---   wind-ocean speed, and a SSS dependent depression of satvpr
!
! ---   Cl based on an approximation
! ---   to values from the COARE 3.0 algorithm (Fairall et al., 2003), 
! ---   for Cl over the global ocean in the range 1m/s <= Va <= 40m/s
! ---   and -8degC <= Ta-Ts <= 7degC, that is quadratic in Ta-Ts and
! ---   quadratic in either Va or 1/Va (Kara et al.,  2005).
!
! ---   Fairall, C. W., E. F. Bradley, J. E. Hare, A. A. Grachev, and J. B.
! ---   Edson, 2003:  Bulk parameterization of air-sea fluxes:  Updates 
! ---   and verification for the COARE algorithm.  J. Climate, 16, 571-591.
!
! ---   Kara, A. B., H. E. Hurlburt, and A. J. Wallcraft, 2005:
! ---   Stability-dependent exchange coefficients for air-sea fluxes.
! ---   J. Atmos. Oceanic. Technol., 22, 1080-1094. 
! ---   http://www7320.nrlssc.navy.mil/pubs.php
!
! ---   use virtual temperature for density
        rair = pair / (rgas * ( tzero + airt ) * (1.0+0.608*vpmx) )
        slat = evaplh*rair
        ssen = csubp *rair
!
        tdif  = temp(i,j,1,n) - airt
! ---   correct tamts for 100% humidity
        sssf  = 1.0
        tamts = -tdif - 0.608*(airt+tzero)* &
                              (qsatur6(airt,pair,sssf)-vpmx)
        tamts = min( tdmax, max( tdmin, tamts ) )
        va    = min( vamax, max( vamin,  wind ) )  !wind=samo
        if     (va.le.5.0) then
          if     (tamts.gt. 0.75) then !stable
            clh =   (as0_00 + as0_01* va + as0_02* va**2) &
                  + (as0_10 + as0_11* va + as0_12* va**2)*tamts &
                  + (as0_20 + as0_21* va + as0_22* va**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au0_00 + au0_01* va + au0_02* va**2) &
                  + (au0_10 + au0_11* va + au0_12* va**2)*tamts &
                  + (au0_20 + au0_21* va + au0_22* va**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap0_00 + ap0_01* va + ap0_02* va**2) &
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am0_00 + am0_01* va + am0_02* va**2) &
                     + q *(an0_00 + an0_01* va + an0_02* va**2)
          endif !tamts
        else !va>5
          qva = 1.0/va
          if     (tamts.gt. 0.75) then !stable
            clh =   (as5_00 + as5_01* va + as5_02* va**2) &
                  + (as5_10 + as5_11*qva + as5_12*qva**2)*tamts &
                  + (as5_20 + as5_21*qva + as5_22*qva**2)*tamts**2
          elseif (tamts.lt.-0.75) then !unstable
            clh =   (au5_00 + au5_01* va + au5_02* va**2) &
                  + (au5_10 + au5_11*qva + au5_12*qva**2)*tamts &
                  + (au5_20 + au5_21*qva + au5_22*qva**2)*tamts**2
          elseif (tamts.ge.-0.098)  then
            q = (tamts-0.75)/0.848  !linear between  0.75 and -0.098
            q = q**2  !favor  0.75
            clh = (1.0-q)*(ap5_00 + ap5_01* va + ap5_02* va**2 &
                                  + ap5_11*qva + ap5_12*qva**2) &
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          else
            q = (tamts+0.75)/0.652  !linear between -0.75 and -0.098
            q = q**2  !favor -0.75
            clh = (1.0-q)*(am5_00 + am5_01* va + am5_02* va**2 &
                                  + am5_11*qva + am5_12*qva**2) &
                     + q *(an5_00 + an5_01* va + an5_02* va**2)
          endif !tamts
        endif !va
        csh  = 0.9554*clh
!
! ---   sssf = fractional depression of satvpr due to SSS
! ---            Sud and Walker (1997), from Witting (1908)
! ---   evap = evaporation (W/m^2) into atmos from ocean.
        sssf = 1.0 - (0.001/1.805)*max(saln(i,j,1,n)-0.03,0.0)
        if     (empflg.lt.0) then
          evape = slat*clh*wind*(qsatur6(esst,pair,sssf)-vpmx)
        endif
! ---   Latent Heat flux (W/m2)
        if(cpl_latflx) then
           evap=imp_latflx(i,j,1)
        else
           evap=slat*clh*wind*(qsatur6(temp(i,j,1,n),pair,sssf)-vpmx)
        endif
! ---   snsibl = sensible heat flux  (W/m^2) into atmos from ocean.
        if(cpl_sensflx) then
           snsibl=imp_sensflx(i,j,1)
        else
           snsibl = ssen*csh*wind* tdif
        endif
! ---   surflx = thermal energy flux (W/m^2) into ocean
        surflx(i,j) = radfl - snsibl - evap
!
!diag   if     (i.eq.itest.and.j.eq.jtest) then
!diag     write(lp,'(i9,2i5,a,3f8.5)') &
!diag     nstep,i0+i,j0+j,' cl,cs,cd    = ',clh,csh,cd0
!diag     write(lp,'(i9,2i5,a,2f8.2,f8.5)') &
!diag     nstep,i0+i,j0+j,' va,tamst,ustar = ',va,tamts,ustar(i,j)
!diag     call flush(lp)
!diag   endif
      elseif (flxflg.eq.5) THEN
! ---   CORE v2 Large and Yeager 2009 CLym. Dyn.: The global climatology 
! ---    of an interannually varying air-sea flux dataset.
! ---   The bulk formulae effectively transform the problem of specifying
! ---   the turbulent surface fluxes (at 10m) into one of describing the
! ---   near surface atmospheric state (wind, temperature and humidity).
! ---   Note that vpmx actually contains specific humidity
!
        rair = pairc / (rgas * ( tzero + airt ))  !always uses pairc
        qrair=1.0/rair
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Over-ocean fluxes following Large and Yeager (used in NCAR models)
! !
! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
!
! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
! Stephen.Griffies@noaa.gov updated the code with the bug fix.
! Script to find the cd,ch,ce to transform fluxes at 10m to surface
! fluxes
! See Large and Yeager 2004 for equations :
! "Diurnal to Decadal Global Forcing For Ocean and Sea-Ice Models:The
! Data Sets
!  and Flux Climatologies", NCAR Technical report.
! The  code below is for values at zi = 10m high
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!
        zi=10.0  ! height in m
        tv = (airt+tzero)*(1.0+0.608*vpmx) !in Kelvin
        uw = max(wind, 0.5)      !0.5 m/s floor on wind (undocumented NCAR)
        uw10 = uw                !first guess 10m wind
   
        cd_n10 = (2.7/uw10+0.142+0.0764*uw10)*1.e-3         !L-Y eqn. 6a
        cd_n10_rt = sqrt(cd_n10)
        ce_n10 = 34.6 *cd_n10_rt*1.e-3                      !L-Y eqn. 6b
        stab   = 0.5 + sign(0.5,airt-temp(i,j,1,n))
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt*1.e-3  !L-Y eqn. 6c

        cd10 = cd_n10  !first guess for exchange coeff's at z
        ch10 = ch_n10
        ce10 = ce_n10
        do it_a= 1,3  !Monin-Obukhov iteration
          cd_rt = sqrt(cd10)
          ustar1 = cd_rt*uw                              !L-Y eqn. 7a
          tstar = (ch10/cd_rt)*(airt-temp(i,j,1,n))      !L-Y eqn. 7b
          qstar = (ce10/cd_rt)*(vpmx-qsatur5(temp(i,j,1,n),qrair))  !L-Y eqn. 7c
          bstar = g*(tstar/tv+qstar/(vpmx+1.0/0.608))
          zeta  = vonkar*bstar*zi/(ustar1*ustar1)        !L-Y eqn. 8a
          zeta  = sign( min(abs(zeta),10.0), zeta )      !undocumented NCAR
          x2 = sqrt(abs(1.0-16.0*zeta))                  !L-Y eqn. 8b
          x2 = max(x2, 1.0)                              !undocumented NCAR
          x  = sqrt(x2)
          if (zeta > 0.0) then
              psi_m = -5.0*zeta  !L-Y eqn. 8c
              psi_h = -5.0*zeta  !L-Y eqn. 8c
          else
              psi_m = log((1.0+2.0*x+x2)*(1.0+x2)/8.0) &
                       -2.0*(atan(x)-atan(1.0))  !L-Y eqn. 8d
              psi_h = 2.0*log((1.0+x2)/2.0)      !L-Y eqn. 8e
          end if
          uw10 = uw/(1.0+cd_n10_rt*(log(zi/10)-psi_m)/vonkar)  !L-Y eqn. 9
          cd_n10 = (2.7/uw10+0.142+0.0764*uw10)*1.e-3          !L-Y eqn. 6a again
          cd_n10_rt = sqrt(cd_n10)
          ce_n10 = 34.6*cd_n10_rt*1.e-3                        !L-Y eqn. 6b again
          stab   = 0.5 + sign(0.5,zeta)
          ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt*1.e-3   !L-Y eqn. 6c again
          z0 = 10*exp(-vonkar/cd_n10_rt) ! diagnostic
          xx   = (log(zi/10.0)-psi_m)/vonkar
          cd10 = cd_n10/(1.0+cd_n10_rt*xx)**2         !L-Y 10a
          xx   = (log(zi/10.0)-psi_h)/vonkar
          ch10 = ch_n10/(1.0+ch_n10*xx/cd_n10_rt) * &
                              sqrt(cd10/cd_n10)       !L-Y 10b
          ce10 = ce_n10/(1.0+ce_n10*xx/cd_n10_rt) * &
                              sqrt(cd10/cd_n10)       !L-Y 10c
        end do !it_a

! ---  Latent Heat flux
        slat = evaplh*rair
        if     (empflg.lt.0) then
          evape = slat*ce10*wind*(qsatur5(esst,qrair)-vpmx)
        endif
        if(cpl_latflx) then
          evap=imp_latflx(i,j,1)
        else
          evap = slat*ce10*wind*(qsatur5(temp(i,j,1,n),qrair)-vpmx)
        endif

! --- Sensible Heat flux
        ssen   = cpcore*rair
        if(cpl_sensflx) then
          snsibl=imp_sensflx(i,j,1)
        else
          snsibl = ssen*ch10*wind*(temp(i,j,1,n)-airt)
        endif

! --- Total surface fluxes
        surflx(i,j) = radfl - snsibl - evap
      elseif (flxflg.eq.3) then
!
! ---   input radiation flux is the net flux.
!
        evap=0.0
        surflx(i,j)=radfl
      else  ! no flux
        evap=0.0
        surflx(i,j)=0.0  
      endif  ! flxflg
!
! --- add a time-invarient net heat flux offset
      if     (flxoff) then
        surflx(i,j)=surflx(i,j)+offlux(i,j)
      endif
!
! --- relax to surface temperature
      if     (sstflg.ge.1) then 
! ---   use a reference relaxation thickness (min. mixed layer depth)
! ---   in shallow water, thkmlt is replaced by the total depth
! ---   actual e-folding time is (dpmixl(i,j,n)/(thkmlt*onem))/rmut
! ---   in shallow water this is (dpmixl(i,j,n)/p(i,j,kk+1)  )/rmut
        if     (sstflg.eq.1) then !climatological sst
        if     (natm.eq.2) then
          sstdif = ( twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1) - &
                   temp(i,j,1,n)
          else
          sstdif = ( twall(i,j,1,lc0)*wc0+twall(i,j,1,lc1)*wc1 &
                    +twall(i,j,1,lc2)*wc2+twall(i,j,1,lc3)*wc3) - &
                   temp(i,j,1,n)
          endif !natm
        else  !synoptic sst
          if(cpl_seatmp) then
            sstdif = imp_seatmp(i,j,1) - temp(i,j,1,n)
          elseif (natm.eq.2) then
            sstdif = ( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1) - &
                     temp(i,j,1,n)
          else
            sstdif = ( seatmp(i,j,l0)*w0+seatmp(i,j,l1)*w1 &
                      +seatmp(i,j,l2)*w2+seatmp(i,j,l3)*w3) - &
                     temp(i,j,1,n)
          endif !natm
        endif
        sstrlx=(rmut*spcifh*min(p(i,j,kk+1),thkmlt*onem)/g)*sstdif
        surflx(i,j)=surflx(i,j)+sstrlx
        sstflx(i,j)=sstflx(i,j)+sstrlx
      endif !sstflg
! --- sswflx = shortwave radiative energy flux (W/m^2) into ocean
      sswflx(i,j)=swfl
! --- emnp = evaporation minus precipitation   (m/sec) into atmos.
      if     (.not.pcipf) then
        prcp = 0.0
        emnp = 0.0   !no E-P
      elseif (empflg.eq.3) then
        emnp = -prcp !input prcp is P-E
      elseif (empflg.gt.0) then !E based on model SST
        emnp = evap *svref/evaplh - prcp  !input prcp is P
      else  !E based on observed SST
        emnp = evape*svref/evaplh - prcp  !input prcp is P
      endif
! --- wtrflx = water flux (m/s kg/m**3) into ocean
      wtrflx(i,j)=-emnp*rhoref
! --- allow for rivers as a precipitation bogas (m/s kg/m**3)
      if     (priver) then
        if(cpl_orivers.and.cpl_irivers) then
            rivflx(i,j) = (imp_orivers(i,j,1)+imp_irivers(i,j,1)) &
                        * rhoref
        else
            rivflx(i,j) = ( rivers(i,j,lr0)*wr0+rivers(i,j,lr1)*wr1    &
                        +   rivers(i,j,lr2)*wr2+rivers(i,j,lr3)*wr3)   &
                        * rhoref
        endif
!       wtrflx(i,j) = wtrflx(i,j)+rivflx(i,j) !update wtrflx in thermf_oi
      else
        rivflx(i,j) = 0.0
      endif
! --- relax to surface salinity
      if     (srelax) then
! ---   use a reference relaxation thickness (min. mixed layer depth)
! ---   in shallow water, thkmls is replaced by the total depth
! ---   actual e-folding time is (dpmixl(i,j,n)/(thkmls*onem))/rmus
! ---   in shallow water this is (dpmixl(i,j,n)/p(i,j,kk+1)  )/rmus
!
! ---   sssrmx is the maximum SSS difference for relaxation (psu)
! ---          set negative to stop relaxing entirely at -sssrlx
! ---          the default (sssflg=1) is 99.9, i.e. no limit
! ---          always fully relax when ice covered or when
! ---          fresher than half climatological sss.
!
        sssc   =  swall(i,j,1,lc0)*wc0+swall(i,j,1,lc1)*wc1 &
                 +swall(i,j,1,lc2)*wc2+swall(i,j,1,lc3)*wc3
        sssdif = sssc - saln(i,j,1,n)
        if     (saln(i,j,1,n).gt.0.5*sssc .and. &
                abs(sssdif).gt.abs(sssrmx(i,j))) then  !large sss anomaly
          if     (sssrmx(i,j).lt.0.0) then
            sssdif = covice(i,j)*sssdif !turn off relaxation except under ice
          elseif (sssdif.lt.0.0) then !sssdif < -sssrmx < 0
            sssdif =      covice(i,j)*   sssdif + &
                     (1.0-covice(i,j))*(-sssrmx(i,j))  !limit relaxation
          else !sssdif > sssrmx > 0
            sssdif =      covice(i,j)*    sssdif + &
                     (1.0-covice(i,j))*   sssrmx(i,j)  !limit relaxation
          endif
        endif
        sssflx(i,j)=(rmus(i,j)*min(p(i,j,kk+1),thkmls*onem)/g)* &
                    sssdif
         util2(i,j)=sssflx(i,j)*scp2(i,j)
         util1(i,j)=max(util2(i,j),0.0)
!       salflx(i,j)=salflx(i,j)+sssflx(i,j) !update salflx in thermf_oi
      else
        sssflx(i,j)=0.0
         util2(i,j)=0.0
         util1(i,j)=0.0
      endif !srelax
      endif !ip
      enddo !i
      return
      end subroutine thermfj

      subroutine thermf_diurnal(diurnal, date)
      implicit none
!
      real        diurnal(0:24,-91:91),date
!
! --- Calculate a table of latitude vs hourly scale factors
! --- for the distribution of daily averaged solar radiation
! --- the clear sky insolation formula of Lumb (1964) is used with 
! --- correction for the seasonally varying earth-sun distance.
! --- According to reed (1977) the lumb formula gives values in close
! --- agreement with the daily mean values of the seckel and beaudry 
! --- (1973) formulae derived from data in the smithsonian
! --- meteorological tables --- (list, 1958).
!
! --- Lumb, F. E., 1964: The influence of cloud on hourly amounts of
! --- total solar radiation at sea surface.Quart. J. Roy. Meteor. Soc.
! --- 90, pp43-56.
!
! ---   date = julian type real date - 1.0 (range 0. to 365.), 
! ---          where 00z jan 1 = 0.0.
!
! --- Base on "QRLUMB" created 2-4-81 by Paul J Martin. NORDA Code 322.
!
      real, parameter ::     pi = 3.14159265
      real, parameter :: raddeg = pi/180.0
!
      integer lat,ihr
      real    sindec,cosdec,alatrd,fd,ourang,sinalt,ri,qsum
      real*8  sum
!
!     calc sin and cosin of the declination angle of the sun.
      call declin(date,sindec,cosdec)
!
!     loop through latitudes
      do lat= -90,90
!       calc latitude of site in radians.
        alatrd = lat*raddeg
!
!       loop through hours
        sum = 0.0
        do ihr= 0,23
!         calc hour angle of the sun (the angular distance of the sun
!         from the site, measured to the west) in radians.
          fd     = real(ihr)/24.0
          ourang = (fd-0.5)*2.0*pi
!         calc sine of solar altitude.
          sinalt = sin(alatrd)*sindec+cos(alatrd)*cosdec*cos(ourang)
!
!         calc clear-sky solar insolation from lumb formula.
          if     (sinalt.le.0.0) then
            diurnal(ihr,lat) = 0.0
          else
            ri=1.00002+.01671*cos(0.01720242*(date-2.1))
            diurnal(ihr,lat) = 2793.0*ri*ri*sinalt*(.61+.20*sinalt)
          endif
          sum = sum + diurnal(ihr,lat)
        enddo !ihr
        if     (sum.gt.0.0) then
!         rescale so that sum is 24.0 (daily average to diurnal factor)
          qsum = 24.0/sum
          do ihr= 0,23
            diurnal(ihr,lat) = diurnal(ihr,lat)*qsum
          enddo !ihr
        endif
        diurnal(24,lat) = diurnal(0,lat) !copy for table lookup
      enddo !lat
      do ihr= 0,24
        diurnal(ihr,-91) = diurnal(ihr,-90) !copy for table lookup
        diurnal(ihr, 91) = diurnal(ihr, 90) !copy for table lookup
      enddo !ihr
      return
!
      contains
        subroutine declin(date,sindec,cosdec)
        implicit none
!
        real date,sindec,cosdec
!
!  subroutine to calc the sin and cosin of the solar declination angle
!  as a function of the date.
!       date = julian type real date - 1.0 (range 0. to 365.), where 00z
!              jan 1 = 0.0.
!       sindec = returned sin of the declination angle.
!       cosdec = returned cosin of the declination angle.
!  formula is from fnoc pe model.
!  created 10-7-81.   paul j martin.   norda code 322.
!
        real a
!
        a=date
        sindec=.39785*sin(4.88578+.0172*a+.03342*sin(.0172*a)- &
        .001388*cos(.0172*a)+.000348*sin(.0344*a)-.000028*cos(.0344*a))
        cosdec=sqrt(1.-sindec*sindec)
        return
        end subroutine declin
      end subroutine thermf_diurnal

      subroutine swfrac_ij(akpar,zz,kz,zzscl,jerlov,swfrac)
      implicit none
!
      integer kz,jerlov
      real    akpar,zz(kz),zzscl,swfrac(kz)
!
! --- calculate fraction of shortwave flux remaining at depths zz
!
! --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
! --- zzscl  = scale factor to convert zz to m
! --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
!
! --- zz(kz) must be the bottom, so swfrac(kz)=0.0 and
! --- any residual which would otherwise be at the bottom is 
! --- uniformly distrubuted across the water column
!
! --- betard is jerlov water types 1 to 5 red  extinction coefficient
! --- betabl is jerlov water types 1 to 5 blue extinction coefficient
! --- redfac is jerlov water types 1 to 5 fract. of penetr. red light
      real, parameter, dimension(5) :: &
        betard = (/ 1.0/0.35, 1.0/0.6, 1.0,     1.0/1.5, 1.0/1.4  /), &
        betabl = (/ 1.0/23.0, 1.0/20.0,1.0/17.0,1.0/14.0,1.0/ 7.9 /), &
        redfac = (/ 0.58,     0.62,    0.67,    0.77,    0.78     /)
!
! --- parameters for ZP LEE et al., 2005 SW attenuation scheme
      real, parameter ::  jjx0 = -0.057
      real, parameter ::  jjx1 =  0.482
      real, parameter ::  jjx2 =  4.221
      real, parameter ::  jjc0 =  0.183
      real, parameter ::  jjc1 =  0.702
      real, parameter ::  jjc2 = -2.567
!
! --- local variables for ZP LEE et al., 2005 SW attenuation scheme
      real chl                  ! surface chlorophyll value (mg m-3)
      real clog                 ! log10 transformed chl
      real a490                 ! total absorption coefficient 490 nm
      real bp550                ! particle scattering coefficient 550 nm
      real v1                   ! scattering emprical constant
      real bbp490               ! particle backscattering 490 nm
      real bb490                ! total backscattering coefficient 490 nm
      real k1                   ! internal vis attenuation term
      real k2                   ! internal vis attenuation term
!
      integer k,knot0
      real    beta_b,beta_r,frac_r,frac_b,d,swfbot
!
      if     (jerlov.ge.0) then
        if     (jerlov.gt.0) then
! ---     standard Jerlov
          beta_r = betard(jerlov)
          beta_b = betabl(jerlov)
          frac_r = redfac(jerlov)
          frac_b = 1.0 - frac_r
        else
! ---     Jerlov-like scheme, from Kpar
! ---       A. B. Kara, A. B., A. J. Wallcraft and H. E. Hurlburt, 2005:
! ---       A New Solar Radiation Penetration Scheme for Use in Ocean 
! ---       Mixed Layer Studies: An Application to the Black Sea Using
! ---       a Fine-Resolution Hybrid Coordinate Ocean Model (HYCOM)
! ---       Journal of Physical Oceanography vol 35, 13-32
          beta_r = 1.0/0.5
          beta_b = akpar
          beta_b = max( betabl(1), beta_b)  !time interp. kpar might be -ve
          frac_b = max( 0.27, 0.695 - 5.7*beta_b )
          frac_r = 1.0 - frac_b
        endif
!
! ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
!
        if     (-d*beta_r.gt.-10.0) then
          swfbot=frac_r*exp(-d*beta_r)+ &
                 frac_b*exp(-d*beta_b)
        elseif (-d*beta_b.gt.-10.0) then
          swfbot=frac_b*exp(-d*beta_b)
        else
          swfbot=0.0
        endif
!
! ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
!
! ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k  !deepest level not on the bottom
            exit
          endif
        enddo !k
!
! ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
!
          if     (-d*beta_r.gt.-10.0) then
            swfrac(k)=frac_r*exp(-d*beta_r)+ &
                      frac_b*exp(-d*beta_b)-swfbot*d
          elseif (-d*beta_b.gt.-10.0) then
            swfrac(k)=frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k)=0.0-swfbot*d
          endif
        enddo !k
      else   !jerlov.eq.-1
!
! --- ---------------------------------------------------------------------
! ---   shortwave attneuation scheme from:
! ---    Lee, Z., K. Du, R. Arnone, S. Liew, and B. Penta (2005),
! ---     Penetration of solar radiation in the upper ocean:
! ---     A numerical model for oceanic and coastal waters,
! ---     J. Geophys. Res., 110, C09019, doi:10.1029/2004JC002780.
! ---   This is a 2-band scheme with "frac_r" fixed. However,
! ---    "beta_b" and "beta_r" are now depth dependent.
! ---   Required input to the scheme is the total absorption coefficient
! ---    at the surface for 490 nm waveband (a490, m-1) and the
! ---    total backscattering coefficient at the surface at the same
! ---    waveband (bb490, m-1). 
! ---   However, here simple "CASE 1" relationships between these surface
! ---    optical properties and surface chlorophyll-a (mg m-3) are assumed.
! ---   These assumptions are considered valid for global, basin-scale
! ---    oceanography. However, coastal and regional applications tend
! ---    to be more complex, and a490 and bb490 should be determined
! ---    directly from the satellite data.
! ---   Authored by Jason Jolliff, NRL; week of 14 November 2011
! --- ---------------------------------------------------------------------
!
        frac_r = 0.52
        frac_b = 1.0 - frac_r
!
! ---   a490 as a function of chl, adapted from Morel et al 2007 Kd(490)
! ---   valid range for chl is 0.01 to 100 mg m-3
        chl  = akpar
        chl  = max(chl,  0.01)
        chl  = min(chl,100.0)
        clog = LOG10(chl)
        a490 = 10.0**(clog*clog*clog*(-0.016993) + &
                      clog*clog*0.0756296 + &
                      clog*0.55420 - 1.14881)
!
! ---   bb490 as a function of chl, from Morel and Maritorania 2001;
! ---   0.0012 is the pure water backscatter
! ---   valid range is restricted to 0.02 - 3.0 mg m-3 chl
        chl   = akpar
        chl   = max(chl,0.02)
        chl   = min(chl,3.0)
        clog  = LOG10(chl)
        bp550 = 0.416*chl**0.766
        if (chl .lt. 2.0) then
          v1 = 0.5*(clog-0.3)
        else
          v1 = 0.0
        endif
        bbp490 = (0.002 + 0.01*(0.50 - 0.25*clog)) &
               * (490.0/550.0)**v1 * bp550
        bb490 = bbp490 + 0.0012
!
! ---   functions of a490 and bb490 for beta_b
        k1 = jjx0 + jjx1*sqrt(a490) + jjx2*bb490
        k2 = jjc0 + jjc1*     a490  + jjc2*bb490
!
! ---   how much SW nominally reaches the bottom
        d = zz(kz)*zzscl
!
        beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
        beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
        if     (-d*beta_b.gt.-10.0) then
          swfbot = frac_r*exp(-d*beta_r)+ &
                   frac_b*exp(-d*beta_b)
        else
          swfbot = 0.0
        endif
!
! ---   spread swfbot uniformly across the water column
        swfbot = swfbot/d
!
! ---   no SW actually left at the bottom
        knot0 = 0
        do k= kz,1,-1
          if     (zz(k).ge.zz(kz)) then
            swfrac(k) = 0.0
          else
            knot0 = k
            exit
          endif
        enddo !k
!
! ---   how much SW reaches zz
        do k= 1,knot0
          d = zz(k)*zzscl
!
          beta_r = 0.560 + 2.34/((0.001 + d)**0.65)
          beta_b = k1    + k2 * 1.0/sqrt(1.0 + d)
          if     (-d*beta_b.gt.-10.0) then
            swfrac(k) = frac_r*exp(-d*beta_r)+ &
                        frac_b*exp(-d*beta_b)-swfbot*d
          else
            swfrac(k) = 0.0-swfbot*d
          endif
        enddo !k
      endif
      end
      subroutine swfrml_ij(akpar,hbl,bot,zzscl,jerlov,swfrml)
      implicit none
!
      integer jerlov
      real    akpar,hbl,bot,zzscl,swfrml
!
! --- calculate fraction of shortwave flux remaining at depth hbl
!
! --- akpar  = CHL (jerlov=-1) or KPAR (jerlov=0)
! --- zzscl  = scale factor to convert hbl and bot to m
! --- jerlov = 1-5 for jerlov type, or 0 for KPAR or -1 for CHL
!
      real zz(2),swf(2)
!
      if     (hbl.ge.bot) then
        swfrml = 0.0
      else
        zz(1) = hbl
        zz(2) = bot
        call swfrac_ij(akpar,zz,2,zzscl,jerlov,swf)
        swfrml = swf(1)
      endif
      return
      end
!
!
!> Revision history:
!>
!> Oct. 1999 - surface flux calculations modified for kpp mixed layer model,
!>             including penetrating solar radiation based on jerlov water type
!> Apr. 2000 - conversion to SI units
!> Oct. 2000 - added thermfj to simplify OpenMP logic
!> Dec. 2000 - modified fluxes when ice is present
!> Dec. 2000 - added Kara bulk air-sea flux parameterization (flxflg=2)
!> May  2002 - buoyfl now calculated in mixed layer routine
!> Aug. 2002 - added nested velocity relaxation
!> Nov. 2002 - separate sss and sst relaxation time scales (thkml[st])
!> Nov. 2002 - save sssflx and sstflx for diagnostics
!> Mar. 2003 - longwave radiation correction for model vs "longwave" SST
!> May  2003 - use seatmp in place of twall.1, when available
!> Mar. 2003 - add option to smooth surface fluxes
!> Mar. 2004 - added epmass for treating E-P as a mass exchange
!> Mar. 2005 - limit thkml[st] to no more than the actual depth
!> Mar. 2005 - added empflg
!> Mar. 2005 - replaced qsatur with 97% of qsatur in evap calculation
!> Mar. 2005 - added ustflg
!> Mar. 2005 - added flxoff
!> Apr. 2005 - add a virtual temperature correction to Ta-Ts for flxflg=4.
!> June 2006 - explicit separation of ocean and sea ice surface fluxes
!> June 2007 - rebalance velocity after sidewall and nestwall relaxation
!> Oct. 2008 - add dswflg
!> June 2009 - add sssrmx
!> Apr. 2010 - change sssrmx to an array
!> Nov. 2010 - added empflg<0 for using observed SST in E
!> Nov. 2011 - ignore sssrmx, i.e. fully relax to sss, under ice
!> July 2013 - vamax set to 34 m/s, same as for Cd (momtum.f)
!> Oct. 2013 - added subroutine swfrac_ij, called in mixed layer routines
!> Oct. 2013 - added subroutine swfrml_ij, called in mixed layer routines
!> Nov. 2013 - added rivflx, so that rivers under sea ice are handled correctly
!> Nov. 2013 - added lwflag=-1 for input radflx=Qlwdn
!> Nov. 2013 - added flxflg=5 for the CORE v2 bulk parameterization
!> Jan. 2014 - tv in Kelvin (flxflg=5)
!> Jan. 2014 - added pair for time varying msl pressure (mslprf)
!> Jan. 2014 - added natm
!> Apr. 2014 - added ice shelf logic (ishelf)
!> Apr. 2014 - replace ip with ipa for mass sums
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> June 2014 - always call thermfj
!> Oct. 2014 - added flxflg=6, similar to flxflg=4
!> Oct. 2014 - flxflg=6 uses sea level pressure optimally
!> Oct. 2014 - flxflg=6 replaces wind speed with wind-ocean speed
!> Oct. 2014 - Newtonian relaxation now uses an implict time step
!> Oct. 2014 - always fully relax when fresher than half climatological sss
!> Aug. 2016 - bugfix for dswflg=1 when lwflag=-1
!> July 2017 - always apply salflx where there is a river
!> Apr. 2018 - for flxflg=5, vpmx contains specific humidity
!> Apr. 2018 - thkmls is -ve for region-wide balanced relaxation (now sssbal)
!> Aug. 2018 - epmass now applies E-P to top layer
!> Aug. 2018 - always use mass-conserving diagnostics
!> Nov. 2018 - only apply E-P (wtrflx) over the ocean
!> Nov. 2018 - virtual salt flux replaced with water and actual salt flux
!> Nov. 2018 - rmus now an array, based on sefold - see blkdat and forfunr
!> Nov. 2018 - added empbal and sssbal
!> Nov. 2018 - added emptgt
!> Dec. 2018 - add /* USE_NUOPC_CESMBETA */ macros for coupled simulation
!> Feb. 2019 - replaced onetai by 1.0
!> Sep. 2019 - added oneta0
!> Oct. 2019 - rmunv replaced with rmunvu and rmunvv
!> Nov. 2019 - added amoflg
