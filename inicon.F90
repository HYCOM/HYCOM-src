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
      subroutine inicon(mnth)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
      use mod_restart    ! HYCOM restart
!
! --- hycom version 1.0
      implicit none
!
      integer mnth
!
! --- ------------------------------------------------------
! --- initializatize all fields (except tracers, see initrc)
! --- ------------------------------------------------------
!
      logical    lpipe_inicon
      parameter (lpipe_inicon=.false.)
!
      real      pinit,pk1p5,pmin(0:kdm),realat,cenlat,tempk,qdep
      integer   i,j,k,k1,kkap,m,n
!diag character text*24
      character ptxt*12,utxt*12,vtxt*12
!
      real     poflat,roflat
      external poflat,roflat
!
# include "stmt_fns.h"
!
      if     (iniflg.eq.3) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in inicon - invalid iniflg value'
        write(lp,*) 'iniflg = ',iniflg
        write(lp,*) 'use restart/src/restart_archv to convert'
        write(lp,*) ' an archive to a restart file (off-line).'
        write(lp,*) 'then rerun with this as restart_in, and with'
        write(lp,*) ' a positive initial value in limits'
        write(lp,*)
        endif !1st tile
        call xcstop('(inicon)')
               stop '(inicon)'
      elseif (iniflg.lt.0 .or. iniflg.gt.3) then
        if     (mnproc.eq.1) then
        write(lp,*)
        write(lp,*) 'error in inicon - invalid iniflg value'
        write(lp,*) 'iniflg = ',iniflg
        write(lp,*)
        endif !1st tile
        call xcstop('(inicon)')
               stop '(inicon)'
      endif
!
! --- set all land to zero
!
      call restart_zero
!
      if     (iniflg.eq.2) then
        call rdrlax(mnth,1)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              do k=1,kk
                if (k.eq.1 .or. k.le.nhybrd) then
                  temp(i,j,k,1)=twall(i,j,k,1)
                  saln(i,j,k,1)=swall(i,j,k,1)
                  th3d(i,j,k,1)=sig(temp(i,j,k,1),saln(i,j,k,1))-thbase
                else  ! isopyc
                  temp(i,j,k,1)=tofsig(theta(i,j,k)+thbase, &
                                       swall(i,j,k,1))
                  saln(i,j,k,1)=swall(i,j,k,1)
                  th3d(i,j,k,1)=theta(i,j,k)
                endif
!
                temp(i,j,k,2)=temp(i,j,k,1)
                saln(i,j,k,2)=saln(i,j,k,1)
                th3d(i,j,k,2)=th3d(i,j,k,1)
              enddo !k
            endif !ip
          enddo !i
        enddo !j
      else
        do j=1,jj
          do k=1,kk
            do i=1,ii
              if (SEA_P) then
                tempk=tofsig(theta(i,j,k)+thbase,saln0)
!
                temp(i,j,k,1)=tempk
                saln(i,j,k,1)=saln0
                th3d(i,j,k,1)=theta(i,j,k)
!
                temp(i,j,k,2)=tempk
                saln(i,j,k,2)=saln0
                th3d(i,j,k,2)=theta(i,j,k)
              endif !ip
            enddo !i
          enddo !k
        enddo !j
      endif
!
      if     (lpipe .and. lpipe_inicon) then
         do k= 1,kk
           write (ptxt,'(a9,i3)') 'temp.1 k=',k
           call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'temp.2 k=',k
           call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') 'saln.1 k=',k
           call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'saln.2 k=',k
           call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') 'th3d.1 k=',k
           call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'th3d.2 k=',k
           call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,2),ip,ptxt)
         enddo
       endif
!
      if     (mnproc.eq.1) then
      write (lp,'('' sigma(k):'',9f7.2/(15x,9f7.2))') &
         (sigma(k),k=1,kk)
      endif !1st tile
      call xcsync(flush_lp)
!
      i = (itdm+1)/2
      j = (jtdm+1)/2
      call xceget(cenlat, plat, i,j)
!
      do j=1,jj
      do i=1,ii
      if (SEA_P) then
      p(i,j,   1)=0.0
      p(i,j,kk+1)=depths(i,j)*onem
!
       if     (dpns.eq.dsns) then
            qdep = 1.0  !not terrain following
       else
            qdep = max( 0.0, min( 1.0, &
                   (depths(i,j) - dsns)/ &
                   (dpns        - dsns)  ) )
      endif

      pmin(0)=0.0
      do k=1,kk
      if     (k.le.nhybrd) then
        pmin(k)=pmin(k-1)+qdep*dp0k(k)+(1.0-qdep)*ds0k(k)
      else  ! isopyc
        pmin(k)=pmin(k-1)
      endif
!
      if     (mxlmy) then
        q2(    i,j,k,1)=smll
        q2l(   i,j,k,2)=smll
        vctymy(i,j,k  )=diwm(i,j)
        difqmy(i,j,k  )=diws(i,j)
        diftmy(i,j,k  )=diws(i,j)
        if (k.eq.kk) then
          q2(    i,j,0  ,1)=smll
          q2l(   i,j,0  ,2)=smll
          q2(    i,j,k+1,1)=smll
          q2l(   i,j,k+1,2)=smll
          vctymy(i,j,0    )=diwm(i,j)
          difqmy(i,j,0    )=diws(i,j)
          diftmy(i,j,0    )=diws(i,j)
          vctymy(i,j,k+1  )=diwm(i,j)
          difqmy(i,j,k+1  )=diws(i,j)
          diftmy(i,j,k+1  )=diws(i,j)
        endif
      endif !mxlmy
!
      if     (iniflg.le.1) then
!
!       initial interfaces from zonal mean climatology.
!
        if (k.lt.kk) then
          if (iniflg.eq.0) then
!
! ---       initial interfaces are flat,
! ---       based on zonal mean climatology at center of the basin.
!
            realat=cenlat
          else  ! iniflg==1
            if (mapflg.ne.4) then
              realat=plat(i,j)
            else
              realat=cenlat
            endif
          endif
          pinit=poflat(.5*(sigma(k)+sigma(k+1)),realat)
!
          if     (i.eq.itest .and. j.eq.jtest) then
             write (lp,'(a,i3,2f12.3,2f10.3)') &
               'k,pmin,poflat,sigma,realat = ', &
               k,pmin(k)*qonem, &
                   pinit*qonem,.5*(sigma(k)+sigma(k+1)),realat
            call flush(lp)
          endif
!
        else  ! k==kk
          pinit=hugel
        endif
        p(i,j,k+1)=max(pmin(k),pinit)
        if     (k.gt.2                .and. &
                k.le.nhybrd+1         .and. &
                p(i,j,k).le.pmin(k-1) .and. &
                (k.eq.kk .or. p(i,j,k+1).gt.pmin(k))) then
          do k1=1,k
            pk1p5 = 0.5*(min(p(i,j,k1)  ,depths(i,j)*onem)+ &
                         min(p(i,j,k1+1),depths(i,j)*onem) )
            th3d(i,j,k1,1)=roflat(pk1p5,realat) -thbase
            temp(i,j,k1,1)=tofsig(th3d(i,j,k1,1)+thbase,saln0)
            saln(i,j,k1,1)=saln0
!
            th3d(i,j,k1,2)=th3d(i,j,k1,1)
            temp(i,j,k1,2)=temp(i,j,k1,1)
            saln(i,j,k1,2)=saln(i,j,k1,1)
!
            if     (kapref.eq.0) then !not thermobaric
              thstar(i,j,k1,1)=th3d(i,j,k1,1)
            elseif (kapref.gt.0) then
              thstar(i,j,k1,1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     kapref)
            else !variable kapref
              thstar(i,j,k1,1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     2)
              thstar(i,j,k1,2)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     kapi(i,j))
            endif
!
            if     (i.eq.itest .and. j.eq.jtest) then
               write (lp,'(a,i3,4f12.3)') &
                 'k,pk+.5,roflat,realat = ', &
                 k1,pk1p5*qonem, &
                 th3d(i,j,k1,1)+thbase,temp(i,j,k1,1),realat
              call flush(lp)
            endif
          end do
        endif
        if (k.eq.kk) then
          do k1=1,kk
            p( i,j,k1+1)=min(p(i,j,k1+1),depths(i,j)*onem)
            dp(i,j,k1,1)=    p(i,j,k1+1)-p(i,j,k1)
            dp(i,j,k1,2)=   dp(i,j,k1,1)
            if     (kapref.eq.0) then !not thermobaric
              thstar(i,j,k1,1)=th3d(i,j,k1,1)
            elseif (kapref.gt.0) then
              thstar(i,j,k1,1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     kapref)
            else !variable kapref
              thstar(i,j,k1,1)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     2)
              thstar(i,j,k1,2)=th3d(i,j,k1,1)+kappaf(temp(i,j,k1,1), &
                                                     saln(i,j,k1,1), &
                                              thbase+th3d(i,j,k1,1), &
#if defined(KAPPAF_CENTERED)
                                       0.5*(p(i,j,k1+1)+p(i,j,k1)), &
#else
                                                        p(i,j,k1), &
#endif
                                                     kapi(i,j))
            endif
          enddo
        endif
      elseif (iniflg.eq.2) then
!
!       initial interfaces from relaxation fields.
!
        if     (k.lt.kk) then
          p(i,j,k+1) = pwall(i,j,k+1,1)
        else
          p(i,j,k+1) = depths(i,j)*onem
        endif
        dp(i,j,k,1) = p(i,j,k+1)-p(i,j,k)
        dp(i,j,k,2) = dp(i,j,k,1)
        if     (kapref.eq.0) then !not thermobaric
          thstar(i,j,k,1)=th3d(i,j,k,1)
        elseif (kapref.gt.0) then
          thstar(i,j,k,1)=th3d(i,j,k,1)+kappaf(temp(i,j,k,1), &
                                               saln(i,j,k,1), &
                                        thbase+th3d(i,j,k,1), &
#if defined(KAPPAF_CENTERED)
                                  0.5*(p(i,j,k+1)+p(i,j,k)), &
#else
                                                  p(i,j,k), &
#endif
                                               kapref)
        else !variable kapref
          thstar(i,j,k,1)=th3d(i,j,k,1)+kappaf(temp(i,j,k,1), &
                                               saln(i,j,k,1), &
                                        thbase+th3d(i,j,k,1), &
#if defined(KAPPAF_CENTERED)
                                  0.5*(p(i,j,k+1)+p(i,j,k)), &
#else
                                                  p(i,j,k), &
#endif
                                               2)
          thstar(i,j,k,2)=th3d(i,j,k,1)+kappaf(temp(i,j,k,1), &
                                               saln(i,j,k,1), &
                                        thbase+th3d(i,j,k,1), &
#if defined(KAPPAF_CENTERED)
                                  0.5*(p(i,j,k+1)+p(i,j,k)), &
#else
                                                  p(i,j,k), &
#endif
                                               kapi(i,j))
        endif
      endif
!
!diag if (mod(k,3).ne.1) go to 55
!diag write (text,'(''intf.pressure (m), k='',i3)') k+1
!diag call prtmsk(ip,p(1-nbdy,1-nbdy,k+1),util1,idm,ii,jj,0.,1.*qonem,text)
!
      enddo !k
!
      if     (isopyc) then
!
! ---   MICOM-like mixed layer no thinner than thkmin.
!
        p( i,j,2)  =max(p(i,j,2),min(depths(i,j),thkmin)*onem)
        dp(i,j,1,1)=p(i,j,2)-p(i,j,1)
        dp(i,j,1,2)=dp(i,j,1,1)
        do k=2,kk
          p( i,j,k+1)=max(p(i,j,k+1),p(i,j,k))
          dp(i,j,k,1)=    p(i,j,k+1)-p(i,j,k)
          dp(i,j,k,2)=dp(i,j,k,1)
        enddo !k
      endif !isopyc
      endif !ip
      enddo !i
      enddo !j
!
      if     (iniflg.eq.0) then
        do k= 1,kk
          tempk = 0.0
          do j= 1,jj
            do i= 1,ii
              if (ip(i,j).eq.1 .and. &
                  abs(th3d(i,j,k,1)-th3d(ii/2,jj/2,k,1)).gt. &
                  abs(tempk)) then
                write(6,*) 'inicon: i,j,k,th3d = ', &
                  i,j,k,th3d(i,j,k,1),th3d(ii/2,jj/2,k,1), &
                        th3d(i,j,k,1)-th3d(ii/2,jj/2,k,1)
                tempk = th3d(i,j,k,1)-th3d(ii/2,jj/2,k,1)
              endif
            enddo
          enddo
          if     (tempk.eq.0.0) then
            write(6,*) 'inicon: constant layer k = ',k
          else
            write(6,*) 'inicon: variable layer k = ',k
          endif
        enddo
      endif
!
!
      do j=1,jj
      do i=1,ii
      if (SEA_P) then
      pbavg(i,j,1)=0.0 - montg_c(i,j)*rhoref ! montg_c non-zero for sshflg=2
      pbavg(i,j,2)=0.0 - montg_c(i,j)*rhoref ! montg_c non-zero for sshflg=2
      pbavg(i,j,3)=0.0 - montg_c(i,j)*rhoref ! montg_c non-zero for sshflg=2
      pbot(i,j)=p(i,j,kk+1)
!
      klist(i,j)=kk  !for MY2.5 mixed layer
!
      steric(i,j)=0.0
      srfhgt(i,j)=0.0
      montg1(i,j)=0.0
!
      do kkap= 1,kapnum
        montg(i,j,1,kkap)=0.0
        do k=1,kk-1
          montg(i,j,k+1,kkap)=montg(i,j,k,kkap)- &
          p(i,j,k+1)*(thstar(i,j,k+1,kkap)-thstar(i,j,k,kkap))*svref**2
        enddo
!
        thkk( i,j,kkap)=thstar(i,j,kk,kkap)
        psikk(i,j,kkap)=montg( i,j,kk,kkap)
      enddo !kkap
      montg1(i,j) = montg(i,j,1,1)
      srfhgt(i,j) = montg1(i,j)
!
! --- start with a thin mixed layer
      if     (hybrid) then
        dpmixl(i,j,1)=min(depths(i,j)*onem-onem, &
                          max(thkmin*onem,p(i,j,2)))
      else  ! isopyc
        dpmixl(i,j,1)=p(i,j,2)
      endif
      dpmixl(i,j,2)=dpmixl(i,j,1)
      dpbl(  i,j)  =dpmixl(i,j,1)
      dpbbl( i,j)  =thkbot*onem
!
      temice(i,j) = temp(i,j,1,1)
      covice(i,j) = 0.0
      thkice(i,j) = 0.0
      endif !ip
      enddo !i
      do i=1,ii
        do k= 1,3
          ubavg(i,j,k) = 0.0
          vbavg(i,j,k) = 0.0
        enddo !k
        do k= 1,kk
          u(i,j,k,1) = 0.0
          u(i,j,k,2) = 0.0
          v(i,j,k,1) = 0.0
          v(i,j,k,2) = 0.0
        enddo !k
      enddo !i
      enddo !j
!
      call xctilr( psikk,1,kapnum,nbdy,nbdy, halo_ps)
      call xctilr(  thkk,1,kapnum,nbdy,nbdy, halo_ps)
      call xctilr( dpmixl,1,2,    nbdy,nbdy, halo_ps)
!
      if     (itest.gt.0 .and. jtest.gt.0) then
         write (lp,103) nstep,i0+itest,j0+jtest, &
         '  istate:  temp    saln  thstar   thkns    dpth   montg', &
         dpmixl(itest,jtest,1)*qonem, &
         (k,temp(itest,jtest,k,1),saln(itest,jtest,k,1), &
          thstar(itest,jtest,k,1)+thbase,dp(itest,jtest,k,1)*qonem, &
          (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem, &
          montg(itest,jtest,k,1)/g,k=1,kk)
         write(lp,104) depths(itest,jtest)
      endif !test tile
      call xcsync(flush_lp)
!
      if     (lpipe .and. lpipe_inicon) then
         do k= 1,kk
           write (ptxt,'(a9,i3)') 'th3d.1 k=',k
           call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'th3d.2 k=',k
           call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') 'thstar k=',k
           call pipe_compare_sym1(thstar(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'saln.1 k=',k
           call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'saln.2 k=',k
           call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') 'temp.1 k=',k
           call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') 'temp.2 k=',k
           call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') '  dp.1 k=',k
           call pipe_compare_sym1(  dp(1-nbdy,1-nbdy,k,1),ip,ptxt)
           write (ptxt,'(a9,i3)') '  dp.2 k=',k
           call pipe_compare_sym1(  dp(1-nbdy,1-nbdy,k,2),ip,ptxt)
           write (ptxt,'(a9,i3)') 'montg  k=',k
           call pipe_compare_sym1(montg(1-nbdy,1-nbdy,k,1),ip,ptxt)
         enddo
         write (ptxt,'(a9,i3)') 'thkk   k=',kk
         call pipe_compare_sym1(thkk( 1-nbdy,1-nbdy,1),ip,ptxt)
         write (ptxt,'(a9,i3)') 'psikk  k=',kk
         call pipe_compare_sym1(psikk(1-nbdy,1-nbdy,1),ip,ptxt)
       endif
!
      if(mxlkrt) then
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              do k=1,kk
                if(dpmixl(i,j,1).gt.p(i,j,k  ) .and. &
                   dpmixl(i,j,1).le.p(i,j,k+1)) then
                  t1sav(i,j,1)=temp(i,j,k,1)
                  s1sav(i,j,1)=saln(i,j,k,1)
                  tmlb( i,j,1)=temp(i,j,k,1)
                  smlb( i,j,1)=saln(i,j,k,1)
                  nmlb( i,j,1)=k
                  t1sav(i,j,2)=t1sav(i,j,1)
                  s1sav(i,j,2)=s1sav(i,j,1)
                  tmlb( i,j,2)=tmlb(i,j,1)
                  smlb( i,j,2)=smlb(i,j,1)
                  nmlb( i,j,2)=k
                endif !dpmixl
              enddo !k
            endif !ip
          enddo !i
        enddo !j
      endif !mxlkrt
!
      if (hybrid) then
        m=2
        n=1
            call pipe_comparall(m,n, 'inicon, step')
        call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                    dpv(1-nbdy,1-nbdy,1,n), &
                    p,depthu,depthv, 0,0)
            if     (lpipe) then
              do k= 1,kk
                write (utxt,'(a9,i3)') 'dpu    k=',k
                write (vtxt,'(a9,i3)') 'dpv    k=',k
                call pipe_compare_sym2(dpu(1-nbdy,1-nbdy,1,n),iu,utxt, &
                                       dpv(1-nbdy,1-nbdy,1,n),iv,vtxt)
              enddo
            endif
        call hybgen(m,n, .false.)
            call pipe_comparall(m,n, 'inicn1, step')
        m=1
        n=2
        call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                    dpv(1-nbdy,1-nbdy,1,n), &
                    p,depthu,depthv, 0,0)
            if     (lpipe) then
              do k= 1,kk
                write (utxt,'(a9,i3)') 'dpu    k=',k
                write (vtxt,'(a9,i3)') 'dpv    k=',k
                call pipe_compare_sym2(dpu(1-nbdy,1-nbdy,1,n),iu,utxt, &
                                       dpv(1-nbdy,1-nbdy,1,n),iv,vtxt)
              enddo
            endif
        call hybgen(m,n, .false.)
            call pipe_comparall(m,n, 'inicn2, step')
      endif
!
      if     (itest.gt.0 .and. jtest.gt.0) then
         write (lp,103) nstep,i0+itest,j0+jtest, &
         '  istate:  temp    saln  thstar   thkns    dpth   montg', &
         dpmixl(itest,jtest,1)*qonem, &
         (k,temp(itest,jtest,k,1),saln(itest,jtest,k,1), &
          thstar(itest,jtest,k,1)+thbase,dp(itest,jtest,k,1)*qonem, &
          (p(itest,jtest,k+1)+p(itest,jtest,k))*0.5*qonem, &
          montg(itest,jtest,k,1)/g,k=1,kk)
         write(lp,104) depths(itest,jtest)
 103     format (i9,2i5,a/23x,'mxl',32x,     f8.1/ &
                         (23x,i3,2f8.2,f8.2,2f8.1,f8.3))
 104     format (         23x,'bot',32x,     f8.1)
      endif !test tile
      call xcsync(flush_lp)
!
      return
!
      contains
      include 'internal_kappaf.h'
      end !inicon
!
!
!> Revision history:
!>
!> Nov. 1999 - added code to initialize homogeneous values of thermodynamical
!>             variables near the surface
!> May  2000 - conversion to SI units
!> Aug. 2000 - added hybrid and isopycnic vertical coordinate options
!> Mar  2009 - more accurate kappaf, with potential density
!> Mar  2012 - replaced dssk with dpns and dsns
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Aug  2017 - kappaf via internal function
!> Aug  2018 - added onetai
!> Aug  2018 - updated the halo of onetai,psikk,thkk,dpmixl
!> Oct  2018 - centered pressure for kapref available via a compile-time macro
!> Feb  2019 - onetai is 1.0
!> Feb  2019 - montg_c correction to pbavg (see momtum for correction to psikk)
!> Feb  2019 - removed onetai
