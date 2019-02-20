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
      subroutine mxpwp(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- hycom version 2.1
!
      integer m,n
!
! -------------------------------------------------------------------
! --- price-weller-pinkel dynamical instability vertical mixing model
! -------------------------------------------------------------------
!
! --- background diapycnal mixing is provided by the explicit diapycnal
! --- mixing model, subroutine diapf2
!
      integer i,j,k
      real    delp,sigmlj
!
# include "stmt_fns.h"
!
      call xctilr(u(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
      call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
!
! --- diffisuvity/viscosity calculation
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxpwpaj(m,n, j)
      enddo
!$OMP END PARALLEL DO
!
! --- final velocity mixing at u,v points
!
      call xctilr(vcty(1-nbdy,1-nbdy,1),1,kk, 1,1, halo_ps)
      call xctilr(dpbl(1-nbdy,1-nbdy),  1, 1, 1,1, halo_ps)
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxpwpbj(m,n, j)
      enddo
!$OMP END PARALLEL DO
!
! --- mixed layer diagnostics
!
      if (diagno) then
!
! --- diagnose new mixed layer depth based on density jump criterion
!$OMP   PARALLEL DO PRIVATE(j,i,k,sigmlj) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
!
! --- depth of mixed layer base set to interpolated depth where
! --- the density jump is equivalent to a tmljmp temperature jump.
! --- this may not vectorize, but is used infrequently.
              sigmlj = -tmljmp*dsigdt(temp(i,j,1,n),saln(i,j,1,n))
              sigmlj = max(sigmlj,tmljmp*0.1)  !cold-water fix
              do k=2,kk
                if     (p(i,j,k+1).ge.p(i,j,kk+1)-onem) then
                  dpmixl(i,j,n) = p(i,j,k+1)
                  exit !k
                elseif ((th3d(i,j,k,n)-th3d(i,j,1,n)).ge.sigmlj) then
                  dpmixl(i,j,n)=max(dp(i,j,1,n), &
                                    p(i,j,k) + dp(i,j,k,n)* &
                     (th3d(i,j,1,n)+sigmlj-th3d(i,j,k-1,n))/ &
                     (th3d(i,j,k,n) +epsil-th3d(i,j,k-1,n)) )
                  exit
                endif
              enddo !k
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        call xctilr(p(     1-nbdy,1-nbdy,2),1,kk, 1,1, halo_ps)
        call xctilr(dpmixl(1-nbdy,1-nbdy,n),1, 1, 1,1, halo_ps)
!
! --- calculate bulk mixed layer t, s, theta
!
!$OMP   PARALLEL DO PRIVATE(j,i,k,delp) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              dpmixl(i,j,m)=dpmixl(i,j,n)
              tmix(i,j)=temp(i,j,1,n)*dp(i,j,1,n)
              smix(i,j)=saln(i,j,1,n)*dp(i,j,1,n)
              do k=2,kk
                delp=min(p(i,j,k+1),dpmixl(i,j,n)) &
                    -min(p(i,j,k  ),dpmixl(i,j,n))
                tmix(i,j)=tmix(i,j)+delp*temp(i,j,k,n)
                smix(i,j)=smix(i,j)+delp*saln(i,j,k,n)
              enddo !k
              tmix(i,j)=tmix(i,j)/dpmixl(i,j,n)
              smix(i,j)=smix(i,j)/dpmixl(i,j,n)
              thmix(i,j)=sig(tmix(i,j),smix(i,j))-thbase
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! --- calculate bulk mixed layer u
!
!$OMP   PARALLEL DO PRIVATE(j,i,k,delp) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_U) then
              umix(i,j)=u(i,j,1,n)*2.*dpu(i,j,1,n)
              do k=2,kk
                delp= &
                   (min(p(i,j,k+1)+p(i-1,j,k+1), &
                        dpmixl(i,j,n)+dpmixl(i-1,j,n)) &
                   -min(p(i,j,k  )+p(i-1,j,k  ), &
                        dpmixl(i,j,n)+dpmixl(i-1,j,n)))
                umix(i,j)=umix(i,j)+delp*u(i,j,k,n)
              enddo !k
              umix(i,j)=umix(i,j)/(dpmixl(i,j,n)+dpmixl(i-1,j,n))
            endif !iu
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! --- calculate bulk mixed layer v
!
!$OMP   PARALLEL DO PRIVATE(j,i,k,delp) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_V) then
              vmix(i,j)=v(i,j,1,n)*2.*dpv(i,j,1,n)
              do k=2,kk
                delp= &
                   (min(p(i,j,k+1)+p(i,j-1,k+1), &
                        dpmixl(i,j,n)+dpmixl(i,j-1,n)) &
                   -min(p(i,j,k  )+p(i,j-1,k  ), &
                        dpmixl(i,j,n)+dpmixl(i,j-1,n)))
                vmix(i,j)=vmix(i,j)+delp*v(i,j,k,n)
              enddo !k
              vmix(i,j)=vmix(i,j)/(dpmixl(i,j,n)+dpmixl(i,j-1,n))
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
      endif                                           ! diagno
!
      return
      end
      subroutine mxpwpaj(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
!
! --- calculate viscosity and diffusivity
!
      integer i
!
      do i=1,ii
        if (SEA_P) then
          call mxpwpaij(m,n, i,j)
        endif !ip
      enddo !i
!
      return
      end
!
      subroutine mxpwpbj(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
!
! --- final velocity mixing at u,v points
!
      integer i
!
      do i=1,ii
        if (SEA_U) then
          call mxpwpbiju(m,n, i,j)
        endif !iu
      enddo !i
!
      do i=1,ii
        if (SEA_V) then
          call mxpwpbijv(m,n, i,j)
        endif !iv
      enddo !i
!
      return
      end
!
      subroutine mxpwpaij(m,n, i,j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- hycom version 2.1
!
      integer m,n, i,j
!
! ----------------------------------------------
! --- pwp vertical mixing, single j-row (part A)
! ----------------------------------------------
!
! local variables for pwp mixing
      real swfrac(kdm+1)       ! fractional surface shortwave radiation flux
!
      real t1d(kdm),s1d(kdm),th1d(kdm),tr1d(kdm,mxtrcr), &
           dp1d(kdm),p1d(kdm+1),u1d(kdm),v1d(kdm),rig(kdm+1)
!
      real dtemp,dsaln,rib,rigf,rig1,rig2,told,sold,trold,uold,vold, &
           sflux1,tsum,ssum,trsum,usum,vsum,dpsum,tup,sup,thup, &
           alfadt,betads, &
           chl,qoneta
!
      integer k,k1,k2,k3,k10,kmax,kmlb,kmlb1,kintf,ktr,iter,jrlv
!
# include "stmt_fns.h"
!
! -----------------------------------------------------------
! --- set 1-d arrays and locate deepest mass-containing layer
! -----------------------------------------------------------
!
      p1d(1)=0.0
      do k=1,kk
        t1d (k)=temp(i,j,k,n)
        s1d (k)=saln(i,j,k,n)
        th1d(k)=sig(t1d(k),s1d(k))-thbase
        do ktr= 1,ntracr
          tr1d(k,ktr)=tracer(i,j,k,n,ktr)
        enddo
        dp1d(k)=dp(i,j,k,n)
        p1d(k+1)=p1d(k)+dp1d(k)
        u1d (k)=0.5*(u(i  ,j,  k,n)+u(i+1,j  ,k,n))
        v1d (k)=0.5*(v(i  ,j  ,k,n)+v(i  ,j+1,k,n))
      enddo !k
!
      do k=kk,1,-1
        if (dp1d(k).gt.tencm) then
          exit !k
        endif 
      enddo !k
      kmax=max(k,2)  !always consider at least 2 layers
!
      qoneta = 1.0/oneta(i,j,n)
!
! ---------------------------------
! --- distribute surface t,s fluxes
! ---------------------------------
!
! --- forcing of t,s by surface fluxes. flux positive into ocean.
      if (pensol) then
! ---   shortwave flux penetration depends on kpar or chl or jerlov water type.
        if     (jerlv0.le.0) then
          chl =  akpar(i,j,lk0)*wk0+akpar(i,j,lk1)*wk1 &
                +akpar(i,j,lk2)*wk2+akpar(i,j,lk3)*wk3
        endif
        call swfrac_ij(chl,p1d,kdm+1,qonem*oneta(i,j,n), &
                       jerlov(i,j),swfrac)
      endif !pensol
!
      do k=1,kk
        if (thermo .or. sstflg.gt.0 .or. srelax) then
          if (k.eq.1) then
            if (pensol) then
              sflux1=surflx(i,j)-sswflx(i,j)
              dtemp=(sflux1+(1.-swfrac(k+1))*sswflx(i,j))* &
                    delt1*g*qoneta/(spcifh*max(onemm,dp1d(k)))
              if (epmass) then  !only actual salt flux
                dsaln= salflx(i,j)* &
                      delt1*g*qoneta/max(onemm,dp1d(k))
              else  !water flux treated as a virtual salt flux
                dsaln=(salflx(i,j)-wtrflx(i,j)*saln(i,j,1,n))* &
                      delt1*g*qoneta/max(onemm,dp1d(k))
              endif
!diag         if (i.eq.itest.and.j.eq.jtest) then
!diag           write (lp,100) &
!diag             nstep,i+i0,j+j0,k,0.,1.-swfrac(k+1),dtemp,dsaln
!diag           call flush(lp)
!diag         endif
 100          format(i9,2i5,i3,'absorbup,dn,dtemp,dsaln ',2f6.3,2f10.6)
            else !.not.pensol
              dtemp=surflx(i,j)* &
                    delt1*g*qoneta/(spcifh*max(onemm,dp1d(k)))
              if (epmass) then  !only actual salt flux
                dsaln= salflx(i,j)* &
                      delt1*g*qoneta/max(onemm,dp1d(k))
              else  !water flux treated as a virtual salt flux
                dsaln=(salflx(i,j)-wtrflx(i,j)*saln(i,j,1,n))* &
                      delt1*g*qoneta/max(onemm,dp1d(k))
              endif
            endif
          elseif (k.le.kmax) then
            if (pensol) then
              dtemp=(swfrac(k)-swfrac(k+1))*sswflx(i,j)*delt1*g*qoneta/ &
                    (spcifh*max(onemm,dp1d(k)))
              dsaln=0.
!diag         if (i.eq.itest.and.j.eq.jtest) then
!diag           write (lp,100) &
!diag             nstep,i+i0,j+j0,k,1.-swfrac(k),1.-swfrac(k+1),dtemp
!diag           call flush(lp)
!diag         endif
            else !.not.pensol
              dtemp=0.0
              dsaln=0.0
            endif
          else !k.gt.kmax
            dtemp=0.0
            dsaln=0.0
          endif
        else !.not.thermo ...
          dtemp=0.0
          dsaln=0.0
        endif !thermo.or.sstflg.gt.0.or.srelax:else
!
        t1d(k)=    t1d(k)+dtemp
        s1d(k)=max(s1d(k)+dsaln,0.0)  !must be non-negative
        th1d(k)=sig(t1d(k),s1d(k))-thbase
      enddo !k
!
! ----------------------------------------------
! --- Don't use PWP when relaxing to climatology
! ----------------------------------------------
!
      if     (rmu(i,j).ne.0.0) then
        kmlb=kmax
        do k=2,kmax
          if (p1d(k).gt.thkmin*onem) then
            kmlb=k-1
            exit !k
          endif
        enddo !k
        dpbl(i,j)=p1d(kmlb+1)
        do k=1,kmax
          temp(i,j,k,n)=t1d(k)
          saln(i,j,k,n)=s1d(k)
          th3d(i,j,k,n)=sig(t1d(k),s1d(k))-thbase
        enddo !k
        return
      endif
!
! ------------------------------------------
! --- relieve mixed layer static instability
! ------------------------------------------
!
      kmlb=1
      tsum=t1d(1)*dp1d(1)
      ssum=s1d(1)*dp1d(1)
      dpsum=dp1d(1)
      do k=2,kmax
        if (locsig) then
          tup=tsum/dpsum
          sup=ssum/dpsum
          alfadt=0.5*(dsiglocdt(tup,sup,dpsum)+ &
                      dsiglocdt(t1d(k),s1d(k),dpsum))*(tup-t1d(k))
          betads=0.5*(dsiglocds(tup,sup,dpsum)+ &
                      dsiglocds(t1d(k),s1d(k),dpsum))*(sup-s1d(k))
          if (alfadt+betads.gt.0.0) then
            kmlb=k
            tsum=tsum+t1d(k)*dp1d(k)
            ssum=ssum+s1d(k)*dp1d(k)
            dpsum=dpsum+dp1d(k)
          else
            exit !k
          endif
        else
          thup=sig(tsum/dpsum,ssum/dpsum)-thbase
          if (th1d(k).lt.thup) then
            kmlb=k
            tsum=tsum+t1d(k)*dp1d(k)
            ssum=ssum+s1d(k)*dp1d(k)
            dpsum=dpsum+dp1d(k)
          else
            exit !k
          endif
        endif
      enddo !k
!
      if (kmlb.gt.1) then
        t1d(1)=tsum/dpsum
        s1d(1)=ssum/dpsum
        th1d(1)=sig(t1d(1),s1d(1))-thbase
        do k=2,kmlb
          t1d(k)=t1d(1)
          s1d(k)=s1d(1)
          th1d(k)=th1d(1)
          do ktr= 1,ntracr
            tr1d(k,ktr)=1.0
          enddo !ktr
!
!diag     if (i.eq.itest .and. j.eq.jtest) then
!diag       write (lp,101) nstep,i+i0,j+j0,k,kmlb, &
!diag         ' relieve static instability - t,s,th:', &
!diag         t1d(k),s1d(k),tr1d(k,1)
!diag       call flush(lp)
!diag     endif
 101      format (i9,2i5,2i3,a/9x,3f9.4)
        enddo !k
      endif !kmlb>1
!
! --- diagnose depth of mixed layer base and homogenize
      call mlbdep(t1d,s1d,th1d,tr1d,u1d,v1d,p1d,dp1d,kmlb,kmax)
!
! ---------------------------------
! --- bulk richardson number mixing
! ---------------------------------
!
! --- mixing within the layer containing the mixed layer base
      kmlb1=kmlb+1
      tsum=t1d(1)*p1d(kmlb1)
      ssum=s1d(1)*p1d(kmlb1)
      usum=u1d(1)*p1d(kmlb1)
      vsum=v1d(1)*p1d(kmlb1)
      k10=kmlb
      do k=kmlb1,kmax
        k1=k-1
        k2=k+1
        if (locsig) then
          alfadt=dsiglocdt(ahalf*(t1d(k1)+t1d(k)), &
                           ahalf*(s1d(k1)+s1d(k)),p1d(k))* &
                                 (t1d(k1)-t1d(k))
          betads=dsiglocds(ahalf*(t1d(k1)+t1d(k)), &
                           ahalf*(s1d(k1)+s1d(k)),p1d(k))* &
                                 (s1d(k1)-s1d(k))
          rib=-g*svref*p1d(k)*min(0.0,alfadt+betads)/ &
             (onem*max(1.e-8,(u1d(k)-u1d(k1))**2+(v1d(k)-v1d(k1))**2))
        else
          rib=g*svref*p1d(k)*max(0.0,th1d(k)-th1d(k1))/ &
             (onem*max(1.e-8,(u1d(k)-u1d(k1))**2+(v1d(k)-v1d(k1))**2))
        endif
!
! --- if rib indicates instability, mix downward to the next interface
        if (rib.lt.ribc.and.p1d(kk+1)-p1d(k+1).ge.tencm) then
!
          tsum=tsum+t1d(k)*dp1d(k)
          ssum=ssum+s1d(k)*dp1d(k)
          do ktr= 1,ntracr
            tr1d(k,ktr)=1.0
          enddo
          usum=usum+u1d(k)*dp1d(k)
          vsum=vsum+v1d(k)*dp1d(k)
!
          t1d(1)=tsum/p1d(k2)
          s1d(1)=ssum/p1d(k2)
          th1d(1)=sig(t1d(1),s1d(1))-thbase
          u1d(1)=usum/p1d(k2)
          v1d(1)=vsum/p1d(k2)
!
          do k3=2,k
            t1d (k3)=t1d (1)
            s1d (k3)=s1d (1)
            th1d(k3)=th1d(1)
            do ktr= 1,ntracr
              tr1d(k3,ktr)=1.0
            enddo
            u1d (k3)=u1d (1)
            v1d (k3)=v1d (1)
!
!diag       if (i.eq.itest .and. j.eq.jtest .and. k3.eq.k10) then
!diag         write (lp,102) nstep,i+i0,j+j0,k,k1,k2,k3,kmlb, &
!diag           ' bulk ri mixing - rib,t,s,th:',min(1000.0,rib), &
!diag           t1d(k3),s1d(k3),th1d(k3)
!diag         call flush(lp)
!diag       endif
 102        format (i9,2i5,5i3,a/9x,4f9.4)
!
          enddo
          kmlb=k
        else
          exit !k
        endif
      enddo !k
!
! --- diagnose depth of mixed layer base and homogenize
      call mlbdep(t1d,s1d,th1d,tr1d,u1d,v1d,p1d,dp1d,kmlb,kmax)
!
! -------------------------------------
! --- gradient richardson number mixing
! -------------------------------------
!
! --- use array 'vcty' to store gradient Ri mixing factor for u,v mixing
!
      do k=1,kk+1
        vcty(i,j,k)=0.0
      enddo
!
! --- perform up to 5 iterations
      do iter=1,5
!
! --- calculate rig array
!
        do k=kmlb+1,kmax
          k1=k-1
          if (locsig) then
            alfadt=dsiglocdt(ahalf*(t1d(k1)+t1d(k)), &
                             ahalf*(s1d(k1)+s1d(k)),p1d(k))* &
                                   (t1d(k1)-t1d(k))
            betads=dsiglocds(ahalf*(t1d(k1)+t1d(k)), &
                             ahalf*(s1d(k1)+s1d(k)),p1d(k))* &
                                   (s1d(k1)-s1d(k))
            rig(k)=-g*min(dp1d(k1),dp1d(k))*svref* &
                  min(-1.0e-3,alfadt+betads)/(onem* &
                  max( 1.0e-6,(u1d(k1)-u1d(k))**2+(v1d(k1)-v1d(k))**2))
          else
            rig(k)=g*min(dp1d(k1),dp1d(k))*svref* &
                  max(1.0e-3,th1d(k)-th1d(k1))/(onem* &
                  max(1.0e-6,(u1d(k1)-u1d(k))**2+ &
                             (v1d(k1)-v1d(k))**2))
          endif
!diag     if (i.eq.itest .and. j.eq.jtest) then
!diag       write(6,103) nstep,i+i0,j+j0,k,iter,th1d(k1)+thbase, &
!diag                    th1d(k)+thbase, &
!diag                    (u1d(k1)-u1d(k))**2+ &
!diag                    (v1d(k1)-v1d(k))**2,rig(k), &
!diag                     dp1d(k1)/onem,dp1d(k)/onem
!diag       call flush(lp)
!diag     endif
 103      format('rig(k)',i9,2i5,2i3,1p,6e13.5)
        enddo !k
!
! --- identify interface where rig has a vertical minimum at each grid point
        kintf=0
        rig2=huge(rig2)
        do k=kmlb+1,kmax
          if(rig(k).lt.rig2) then
            kintf=k
            rig2=rig(k)
          end if
        enddo !k
!
! --- if selected layer pair is unstable, mix to bring rig up to rigc
! --- store factor rig1 in array vcty for u,v mixing
        if(rig2.lt.rigc) then
          k=kintf
          rig1=1.-rig2/rigc
          vcty(i,j,k)=rig1
!
          rigf=rig1*(t1d(k-1)-t1d(k))
          told=t1d(k-1)
          t1d(k-1)=t1d(k-1)-rigf*dp1d(k  )/max(epsil,dp1d(k-1)+dp1d(k))
          t1d(k  )=t1d(k  )+rigf*dp1d(k-1)/max(epsil,dp1d(k-1)+dp1d(k))
!diag     if (i.eq.itest .and. j.eq.jtest.and.mnproc.eq.1) then
!diag     if(k.gt.15.and.k.lt.22) then
!diag       write(6,104) nstep,i+i0,j+j0,k,rigf,rig1,t1d(k-1),t1d(k), &
!diag                    dp1d(k-1)/onem,dp1d(k)/onem, &
!diag                    dp1d(min(kk,k+1))/onem
!diag       call flush(lp)
!diag     endif
 104      format('rig mixing',i9,2i5,i3,1p,7e13.5)
!
          rigf=rig1*(s1d(k-1)-s1d(k))
          sold=s1d(k-1)
          s1d(k-1)=s1d(k-1)-rigf*dp1d(k  )/max(epsil,dp1d(k-1)+dp1d(k))
          s1d(k  )=s1d(k  )+rigf*dp1d(k-1)/max(epsil,dp1d(k-1)+dp1d(k))
!
          th1d(k-1)=sig(t1d(k-1),s1d(k-1))-thbase
          th1d(k  )=sig(t1d(k  ),s1d(k  ))-thbase
!
          do ktr= 1,ntracr
            rigf=rig1*(tr1d(k-1,ktr)-tr1d(k,ktr))
            trold=tr1d(k-1,ktr)
            tr1d(k-1,ktr)=tr1d(k-1,ktr)-rigf*dp1d(k  )/ &
                      max(epsil,dp1d(k-1)+dp1d(k))
            tr1d(k  ,ktr)=tr1d(k  ,ktr)+rigf*dp1d(k-1)/ &
                      max(epsil,dp1d(k-1)+dp1d(k))
          enddo !ktr
!
          rigf=rig1*(u1d(k-1)-u1d(k))
          uold=u1d(k-1)
          u1d(k-1)=u1d(k-1)-rigf*dp1d(k  )/max(epsil,dp1d(k-1)+dp1d(k))
          u1d(k  )=u1d(k  )+rigf*dp1d(k-1)/max(epsil,dp1d(k-1)+dp1d(k))
!
          rigf=rig1*(v1d(k-1)-v1d(k))
          vold=v1d(k-1)
          v1d(k-1)=v1d(k-1)-rigf*dp1d(k  )/max(epsil,dp1d(k-1)+dp1d(k))
          v1d(k  )=v1d(k  )+rigf*dp1d(k-1)/max(epsil,dp1d(k-1)+dp1d(k))
!
        end if !rig2<rigc
!
      enddo !iter
!
! --- diagnose depth of mixed layer base and homogenize
      call mlbdep(t1d,s1d,th1d,tr1d,u1d,v1d,p1d,dp1d,kmlb,kmax)
!
! ------------------------------------
!     reset mixed layer and 3-d arrays
! ------------------------------------
!
      dpbl(i,j)=p1d(kmlb+1)
!
      do k=1,kmax
        temp(i,j,k,n)=t1d(k)
        saln(i,j,k,n)=s1d(k)
        th3d(i,j,k,n)=sig(t1d(k),s1d(k))-thbase
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tr1d(k,ktr)
        enddo
      enddo
!
      return
      end
!
      subroutine mxpwpbiju(m,n, i,j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- hycom version 2.1
!
      real dpm,usum,rig1,rigf
      integer m,n, i,j
      integer k,kintf
!
! ----------------------------------------------------------------------------
! --- pwp vertical diffusion, single j-row (part A), momentum at u grid points
! ----------------------------------------------------------------------------
!
! --- bulk richardson number mixing
! --- homogenize u between the surface and interface kintf, the closest
! --- interface to the interpolated mixed layer thickness
!
      dpm=0.5*(dpbl(i,j)+dpbl(i-1,j))
      usum=0.0
      kintf=2
      pu(i,j,1)=0.0
!
      do k=1,kk
        pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,k,n)
        if (abs(dpm-pu(i,j,k+1)) .lt. abs(dpm-pu(i,j,k)) .and. &
            depthu(i,j)-pu(i,j,k+1) .gt. tencm) then
          kintf=k+1
          usum=usum+u(i,j,k,n)*dpu(i,j,k,n)
        endif
      enddo
!
      if (kintf.gt.2) then
        u(i,j,1,n)=usum/pu(i,j,kintf)
        do k=2,kintf-1
          u(i,j,k,n)=u(i,j,1,n)
        enddo
      endif
!
! --- gradient richardson number mixing
!
      do k=2,kk
        rig1=0.5*(vcty(i,j,k)+vcty(i-1,j,k))
        if (rig1.gt.0.0) then
          rigf=rig1*(u(i,j,k-1,n)-u(i,j,k,n))
          u(i,j,k-1,n)=u(i,j,k-1,n)-rigf*dpu(i,j,k  ,n)/ &
                       max(epsil,dpu(i,j,k-1,n)+dpu(i,j,k,n))
          u(i,j,k  ,n)=u(i,j,k  ,n)+rigf*dpu(i,j,k-1,n)/ &
                       max(epsil,dpu(i,j,k-1,n)+dpu(i,j,k,n))
        endif
      enddo
!
      return
      end
!
      subroutine mxpwpbijv(m,n, i,j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- hycom version 2.1
!
      real dpm,vsum,rig1,rigf
      integer m,n, i,j
      integer k,kintf
!
! ----------------------------------------------------------------------------
! --- pwp vertical diffusion, single j-row (part A), momentum at v grid points
! ----------------------------------------------------------------------------
!
! --- bulk richardson number mixing
! --- homogenize v between the surface and interface kintf, the closest
! --- interface to the interpolated mixed layer thickness
!
      dpm=0.5*(dpbl(i,j)+dpbl(i,j-1))
      vsum=0.0
      kintf=1
      pv(i,j,1)=0.0
!
      do k=1,kk
        pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,k,n)
        if (abs(dpm-pv(i,j,k+1)) .lt. abs(dpm-pv(i,j,k)) .and. &
            depthv(i,j)-pv(i,j,k+1) .gt. tencm) then
          kintf=k+1
          vsum=vsum+v(i,j,k,n)*dpv(i,j,k,n)
        endif
      enddo
!
      if (kintf.gt.2) then
        v(i,j,1,n)=vsum/pv(i,j,kintf)
        do k=2,kintf-1
          v(i,j,k,n)=v(i,j,1,n)
        enddo
      endif
!
! --- gradient richardson number mixing
!
      do k=2,kk
        rig1=0.5*(vcty(i,j,k)+vcty(i,j-1,k))
        if (rig1.gt.0.0) then
          rigf=rig1*(v(i,j,k-1,n)-v(i,j,k,n))
          v(i,j,k-1,n)=v(i,j,k-1,n)-rigf*dpv(i,j,k  ,n)/ &
                       max(epsil,dpv(i,j,k-1,n)+dpv(i,j,k,n))
          v(i,j,k  ,n)=v(i,j,k  ,n)+rigf*dpv(i,j,k-1,n)/ &
                       max(epsil,dpv(i,j,k-1,n)+dpv(i,j,k,n))
        endif
      enddo
!

      return
      end
!
      subroutine mlbdep(t1d,s1d,th1d,tr1d,u1d,v1d,p1d,dp1d,kmlb,kmax)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- hycom version 2.1
!
      integer k,kmlb,kmax,ktr
      real t1d(kdm),s1d(kdm),th1d(kdm),tr1d(kdm,mxtrcr), &
           u1d(kdm),v1d(kdm),p1d(kdm+1),dp1d(kdm)
      real tsum,ssum,usum,vsum,dpsum
!
# include "stmt_fns.h"
!
! --- -----------------------------------------------------------------------
! --- diagnose depth of the PWP mixed layer base and homogenize to that depth
! --- -----------------------------------------------------------------------
!
! --- set to depth of first interface deeper than thkmin across which
! --- the density jump exceeds 1.0e-4
      kmlb=kmax
      do k=2,kmax
        if ((th1d(k)-th1d(k-1)).ge.1.0e-4 .and. &
             p1d(k).gt.thkmin*onem) then
          kmlb=k-1
          exit !k
        endif
      enddo !k
!
      if (kmlb.gt.1) then
        tsum=t1d(1)*dp1d(1)
        ssum=s1d(1)*dp1d(1)
        usum=u1d(1)*dp1d(1)
        vsum=v1d(1)*dp1d(1)
        dpsum=dp1d(1)
        do k=2,kmlb
          tsum=tsum+t1d(k)*dp1d(k)
          ssum=ssum+s1d(k)*dp1d(k)
          usum=usum+u1d(k)*dp1d(k)
          vsum=vsum+v1d(k)*dp1d(k)
          dpsum=dpsum+dp1d(k)
        enddo
!
        t1d(1)=tsum/dpsum
        s1d(1)=ssum/dpsum
        th1d(1)=sig(t1d(1),s1d(1))-thbase
        do ktr= 1,ntracr
          tr1d(1,ktr)=1.0
        enddo
        u1d(1)=usum/dpsum
        v1d(1)=vsum/dpsum
        do k=2,kmlb
          t1d(k)=t1d(1)
          s1d(k)=s1d(1)
          th1d(k)=th1d(1)
          do ktr= 1,ntracr
            tr1d(k,ktr)=1.0
          enddo
          u1d(k)=u1d(1)
          v1d(k)=v1d(1)
        enddo
      endif
!
      return
      end
!
!
!> Revision history:
!>
!> Mar. 2004:  minimum layer thickness used to calculate gradient Ri
!> Oct. 2010:  replaced two calls to dsiglocdX with one call at mid-point
!> Oct. 2013 - added jerlv0=-1 and call to swfrac_ij
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Aug. 2018 - added wtrflx, salflx now only actual salt flux
!> Nov. 2018 - allow for oneta in swfrac and surface fluxes
