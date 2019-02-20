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
      subroutine mxkrtm(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n
!
! --- hycom version 1.0 (adapted from micom version 2.8)
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) ::&
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::&
#endif
       sdot
!
      integer i,j,k
      real delp,q,thk
!cc   integer kmax
!cc   real totem,tosal,tndcyt,tndcys,work(3)
#if defined(RELO)
!
      if     (.not.allocated(sdot)) then
        allocate( &
                  sdot(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy) )
                  sdot = r_init
      endif
#endif
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
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
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,0p,f8.2,f8.1))
!diag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest, &
!diag   '  entering mxlayr:  temp    saln    dens    thkns    dpth', &
!diag   (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n), &
!diag   th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)*qonem, &
!diag   p(itest,jtest,k+1)*qonem,k=1,kk)
!
      if (thermo .or. sstflg.gt.0 .or. srelax) then
!
! --- -----------------------------------
! --- mixed layer entrainment/detrainment
! --- -----------------------------------
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n,sdot) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxkrtmaj(m,n, sdot, j)
      enddo
!$OMP END PARALLEL DO
!
      else !.not.thermo ...
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do  j=1,jj
        do i=1,ii
          if (SEA_P) then
            surflx(i,j)=0.
            salflx(i,j)=0.
            wtrflx(i,j)=0.
            sdot(i,j)=dp(i,j,1,n)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      end if !thermo.or.sstflg.gt.0.or.srelax:else
!
!diag if (itest.gt.0.and.jtest.gt.0.and.turgen(itest,jtest).lt.0.) &
!diag   write (lp,'(i9,2i5,a,f8.2)') nstep,itest,jtest, &
!diag   '  monin-obukhov length (m):',sdot(itest,jtest)*qonem
!
! --- store 'old' t/s column integral in totem/tosal (diagnostic use only)
!cc   totem=0.
!cc   tosal=0.
!cc   do k=1,kk
!cc   if (max(dp(itest,jtest,1,n)+sdot(itest,jtest),thkmin*onem).gt.
!cc  .  p(itest,jtest,k) .or. max(th3d(itest,jtest,1,m),th3d(itest,
!cc  .  jtest,1,n)) +sigjmp.ge.th3d(i,j,k,n)) then
!cc     kmax=k
!cc     totem=totem+temp(itest,jtest,k,n)*dp(itest,jtest,k,n)
!cc     tosal=tosal+saln(itest,jtest,k,n)*dp(itest,jtest,k,n)
!cc   end if
!cc   end do
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n,sdot) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call mxkrtmbj(m,n, sdot, j)
      enddo
!$OMP END PARALLEL DO
!
! --- compare 'old' with 'new' t/s column integral (diagnostic use only)
!
!cc   tndcyt=-totem
!cc   tndcys=-tosal
!cc   do k=kmax,1,-1
!cc     tndcyt=tndcyt+temp(itest,jtest,k,n)*dp(itest,jtest,k,n)
!cc     tndcys=tndcys+saln(itest,jtest,k,n)*dp(itest,jtest,k,n)
!cc   end do
!cc   write (lp,'(i9,2i5,i3,3x,a,1p,3e10.2/25x,a,3e10.2)') nstep,itest,
!cc  .  jtest,kmax,'total saln,srf.flux,tndcy:',tosal/g,salflx(itest,
!cc  .  jtest)*delt1,tndcys/g,'total temp,srf.flux,tndcy:',totem/g,
!cc  .  surflx(itest,jtest)*delt1,tndcyt*spcifh/g
!
! --- store 'old' interface pressures in -pu,pv-
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
      do k=2,kk+1
      do i=1,ii
      if (SEA_U) then
      pu(i,j,k)=min(depthu(i,j),.5*(p(i,j,k)+p(i-1,j,k)))
      endif !iu
      enddo !i
!
      do i=1,ii
      if (SEA_V) then
      pv(i,j,k)=min(depthv(i,j),.5*(p(i,j,k)+p(i,j-1,k)))
      endif !iv
      enddo !i
      enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
! --- store 'new' layer thicknesses in -dpu,dpv-
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
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
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                  dpv(1-nbdy,1-nbdy,1,n), &
                  p,depthu,depthv, 0,0)
!
! --- redistribute momentum in the vertical.
! --- homogenize (u,v) over depth range defined in -util1,util2-
!
! --- thk>0 activates momentum diffusion across mixed-layer interface
      thk=vertmx*onem*delt1
!
!$OMP PARALLEL DO PRIVATE(j,i,k,delp,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
!
      do i=1,ii
      if (SEA_U) then

      util1(i,j)=max(dpu(i,j,1,n),pu(i,j,2)+thk)
      uflux(i,j)=0.
      util3(i,j)=0.
!
      do k=1,kk
      delp=max(0.,min(util1(i,j),pu(i,j,k+1)) &
                 -min(util1(i,j),pu(i,j,k  )))
      uflux(i,j)=uflux(i,j)+u(i,j,k,n)*delp
      util3(i,j)=util3(i,j)            +delp
      enddo !k
!
      u(i,j,1,n)=uflux(i,j)/util3(i,j)
      endif !iu
      enddo !i
!
      do i=1,ii
      if (SEA_V) then
      util2(i,j)=max(dpv(i,j,1,n),pv(i,j,2)+thk)
      vflux(i,j)=0.
      util4(i,j)=0.
!
      do k=1,kk
      delp=max(0.,min(util2(i,j),pv(i,j,k+1)) &
                 -min(util2(i,j),pv(i,j,k  )))
      vflux(i,j)=vflux(i,j)+v(i,j,k,n)*delp
      util4(i,j)=util4(i,j)            +delp
      enddo !k
!
      v(i,j,1,n)=vflux(i,j)/util4(i,j)
      endif !iv
      enddo !i
!
      do k=2,kk
!
      do i=1,ii
      if (SEA_U) then
      pu(i,j,k)=pu(i,j,k-1)+dpu(i,j,k-1,n)
      q=max(0.,min(1.,(util1(i,j)-pu(i,j,k))/(dpu(i,j,k,n)+epsil)))
      u(i,j,k,n)=u(i,j,1,n)*q+u(i,j,k,n)*(1.-q)
      endif !iu
      enddo !i
!
      do i=1,ii
      if (SEA_V) then
      pv(i,j,k)=pv(i,j,k-1)+dpv(i,j,k-1,n)
      q=max(0.,min(1.,(util2(i,j)-pv(i,j,k))/(dpv(i,j,k,n)+epsil)))
      v(i,j,k,n)=v(i,j,1,n)*q+v(i,j,k,n)*(1.-q)
      endif !iv
      enddo !i
      enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
!diag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest, &
!diag   '  exiting  mxlayr:  temp    saln    dens    thkns    dpth', &
!diag   (k,temp(itest,jtest,k,n),saln(itest,jtest,k,n), &
!diag   th3d(itest,jtest,k,n)+thbase,dp(itest,jtest,k,n)*qonem, &
!diag   p(itest,jtest,k+1)*qonem,k=1,kk)
      return
      end

      subroutine mxkrtmaj(m,n, sdot, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       sdot
!
! --- hycom version 1.0 (adapted from micom version 2.8)
!
      integer i,k,ka
!
      real thknss,ustar3,dpth,ekminv,obuinv,buoyfl,dsgdt,tmn,smn, &
           ex,alf1,alf2,cp1,cp3,ape,cc4,spe,pnew,alfadt,betads,thet
!
      real ea1, ea2, em1, em2, em3, em4, em5
      data ea1, ea2, em1, em2, em3, em4, em5 &
         /0.60,0.30,0.45,2.60,1.90,2.30,0.60/          ! Gaspar coefficients
!
# include "stmt_fns.h"
!
      locsig=.true.
!
! --- -----------------------------------
! --- mixed layer entrainment/detrainment
! --- -----------------------------------
!
      do i=1,ii
      if (SEA_P) then
!
! --- determine turb.kin.energy generation due to wind stirring
! --- ustar computed in subr. -thermf-
! --- buoyancy flux (m**2/sec**3), all fluxes into the ocean
! --- note: surface density increases (column is destabilized) if buoyfl < 0
      thknss=dp(i,j,1,n)
      ustar3=ustar(i,j)**3
      tmn=.5*(temp(i,j,1,m)+temp(i,j,1,n))
      smn=.5*(saln(i,j,1,m)+saln(i,j,1,n))
      dsgdt=dsigdt(tmn,smn)
      buoyfl=-g*svref*(dsigds(tmn,smn)* &
           (-wtrflx(i,j)*saln(i,j,1,n)+salflx(i,j))*svref+ &
                       dsgdt          *surflx(i,j) *svref/spcifh)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- option 1 :   k r a u s  -  t u r n e r    mixed-layer t.k.e.  closure
!
!cc   em=0.8*exp(-p(i,j,2)/(50.*onem))  !   hadley centre choice (orig.: 1.25)
!cc   en=0.15                           !   hadley centre choice (orig.: 0.4)
!cc   thermg=-.5*g*((en+1.)*buoyfl+(en-1.)*abs(buoyfl))*rhoref
!cc   turgen(i,j)=delt1*(2.*em*g*ustar3*rhoref+thknss*thermg)*rhoref**2
!
! --- find monin-obukhov length in case of receding mixed layer (turgen < 0).
! --- the monin-obukhov length is found by stipulating turgen = 0.
! --- store temporarily in 'sdot'.
!
!cc   if (turgen(i,j).lt.0.) then
!cc     sdot(i,j)=-2.*em*g*ustar3/min(-epsil,svref*thermg)
!cc   else
!cc     sdot(i,j)=thknss
!cc   end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- option 2 :    g a s p a r   mixed-layer t.k.e. closure
!
      dpth=thknss*qonem
      ekminv=1./hekman(i,j)
      obuinv=buoyfl/max(epsil,ustar3)
      ex=exp(min(50.,dpth*obuinv))
      alf1=ea1+ea2*max(1.,2.5*dpth*ekminv)*ex
      alf2=ea1+ea2*ex
      cp1=((1.-em5)*(alf1/alf2)+.5*em4)*athird
      cp3=max(0.,(em4*(em2+em3)-(alf1/alf2)*(em2+em3-em3*em5))*athird)
      ape=cp3*ustar3-cp1*dpth*buoyfl
!
      if(ape.lt.0.) then                                       ! detrainment
      turgen(i,j)=(g*delt1*rhoref**3)*ape
      sdot(i,j)=max(thkmin*onem,min(thknss,g*cp3/ &
      (svref*cp1*max(epsil,obuinv))))
!
      else                                                     ! entrainment
      cc4=2.*em4/(em1*em1) * alf1*alf1
      spe=(em2+em3)*ustar3-0.5*dpth*buoyfl
      turgen(i,j)=(g*delt1*rhoref**3)*(sqrt((.5*ape-cp1*spe)**2 &
                       +2.*cc4*ape*spe)-(.5*ape+cp1*spe))/(cc4-cp1)
      sdot(i,j)=thknss
      end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- util1,util2 are used to evaluate pot.energy changes during entrainment
      util1(i,j)=th3d(i,j,1,n)*thknss
      util2(i,j)=th3d(i,j,1,n)*thknss**2
!
! --- find pnew in case of mixed layer deepening (turgen > 0). store in 'sdot'.
! --- entrain as many layers as needed to deplete -turgen-.
!
      do k=2,kk
      ka=k-1
      if (k.eq.2) then
        thstar(i,j,ka,1)=th3d(i,j,ka,n)
      endif
      if (locsig) then
        alfadt=0.5* &
              (dsiglocdt(temp(i,j,ka,n),saln(i,j,ka,n),p(i,j,k))+ &
               dsiglocdt(temp(i,j,k ,n),saln(i,j,k ,n),p(i,j,k)))* &
              (temp(i,j,ka,n)-temp(i,j,k,n))
        betads=0.5* &
              (dsiglocds(temp(i,j,ka,n),saln(i,j,ka,n),p(i,j,k))+ &
               dsiglocds(temp(i,j,k ,n),saln(i,j,k ,n),p(i,j,k)))* &
              (saln(i,j,ka,n)-saln(i,j,k,n))
        thstar(i,j,k,1)=thstar(i,j,ka,1)-alfadt-betads
        thet=thstar(i,j,k,1)
      else
        thet=th3d(i,j,k,n)
      endif
      pnew=(2.*turgen(i,j)+thet*p(i,j,k)**2-util2(i,j))/ &
                 max(epsil,thet*p(i,j,k)   -util1(i,j))
! --- stop iterating for 'pnew' as soon as pnew < k-th interface pressure
      if (pnew.lt.p(i,j,k)) pnew=sdot(i,j)
! --- substitute 'pnew' for monin-obukhov length if mixed layer is deepening
      if (turgen(i,j).ge.0.) sdot(i,j)=pnew
!
      util1(i,j)=util1(i,j)+thet*dp(i,j,k,n)
      util2(i,j)=util2(i,j)+thet*(p(i,j,k+1)**2-p(i,j,k)**2)
      enddo !k
      endif !ip
      enddo !i
      return
      end

      subroutine mxkrtmbj(m,n, sdot, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       sdot
!
! --- hycom version 1.0 (adapted from micom version 2.8)
!
      integer i,k,ktr,num
!
      real tdp(idm),sdp(idm),vsflx(idm)
      real pnew,thknss,t1,s1,tmxl,smxl, &
           dpno,sn,tn,dtemp,dsaln,tnew,snew,z,s_up,a,e,b,f,d,c1msig, &
           cc0,cc3,cc1,cc2,x
!
      real ccubq,ccubr,ccubqr,ccubs1,ccubs2,ccubrl,ccubim,root,root1, &
           root2,root3
!
# include "stmt_fns.h"
!
! --- cubic eqn. solver used in mixed-layer detrainment
      ccubq(s)=athird*(cc1/cc3-athird*(cc2/cc3)**2)
      ccubr(s)=athird*(.5*(cc1*cc2)/(cc3*cc3)-1.5*cc0/cc3) &
             -(athird*cc2/cc3)**3
      ccubqr(s)=sqrt(abs(ccubq(s)**3+ccubr(s)**2))
      ccubs1(s)=sign(abs(ccubr(s)+ccubqr(s))**athird,ccubr(s)+ccubqr(s))
      ccubs2(s)=sign(abs(ccubr(s)-ccubqr(s))**athird,ccubr(s)-ccubqr(s))
      root(s)=ccubs1(s)+ccubs2(s)-athird*cc2/cc3
      ccubrl(s)=sqrt(max(0.,-ccubq(s))) &
                *cos(athird*atan2(ccubqr(s),ccubr(s)))
      ccubim(s)=sqrt(max(0.,-ccubq(s))) &
                *sin(athird*atan2(ccubqr(s),ccubr(s)))
      root1(s)=2.*ccubrl(s)-athird*cc2/cc3
      root2(s)=-ccubrl(s)+sqrt(3.)*ccubim(s)-athird*cc2/cc3
      root3(s)=-ccubrl(s)-sqrt(3.)*ccubim(s)-athird*cc2/cc3
!
      do i=1,ii
      if (SEA_P) then
        if (epmass) then  !only actual salt flux
          vsflx(i)= salflx(i,j)
        else !water flux treated as a virtual salt flux
          vsflx(i)=(salflx(i,j)-wtrflx(i,j)*saln(i,j,1,n))
        endif
! --- store (pnew - pold) in 'sdot'.
! --- don't allow mixed layer to get too deep or too shallow.
      sdot(i,j)=min(p(i,j,kk+1),max(thkmin*onem,sdot(i,j)))- &
                dp(i,j,1,n)
      klist(i,j)=2
      tdp(i)=0.
      sdp(i)=0.
!
      do k=2,kk
      pnew=dp(i,j,1,n)+sdot(i,j)
! --- 'tdp,sdp' will be needed for temp./salin. mixing during entrainment
      tdp(i)=tdp(i)+temp(i,j,k,n)*(min(pnew,p(i,j,k+1)) &
                                 -min(pnew,p(i,j,k  )))
      sdp(i)=sdp(i)+saln(i,j,k,n)*(min(pnew,p(i,j,k+1)) &
                                 -min(pnew,p(i,j,k  )))
!
! --- if sdot > 0, remove water from layers about to be entrained.
      dpo(i,j,k,n)=dp(i,j,k,n)					! diapyc.flux
      dp( i,j,k,n)=max(p(i,j,k+1),pnew)-max(p(i,j,k),pnew)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,k,n)-dpo(i,j,k,n))	! diapyc.flux
      if (pnew.ge.p(i,j,k+1)) then
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=0.
        enddo !ktr
      endif
!
! --- if sdot < 0, mixed layer water will be detrained into isopycnic layer
! --- defined in -klist-. to prevent odd/even time step decoupling of mixed-
! --- layer depth, determine -klist- from layer one -th3d- at 2 consecutive
! --- time levels
!
      if (max(th3d(i,j,1,m),th3d(i,j,1,n))+sigjmp.ge.th3d(i,j,k,n))  &
       klist(i,j)=k+1
!
! --- set t/s in massless layers. step 1: copy salinity from layer(s) above
!
      saln(i,j,k,n)=(saln(i,j,k,n)*dp(i,j,k,n)+saln(i,j,k-1,n)*epsil)/ &
                    (              dp(i,j,k,n)+                epsil)
      enddo !k
!
! --- set t/s in massless layers. step 2: copy salinity from layer(s) below
!
      do k=kk-1,2,-1
      saln(i,j,k,n)=(saln(i,j,k,n)*dp(i,j,k,n)+saln(i,j,k+1,n)*epsil)/ &
                    (              dp(i,j,k,n)+                epsil)
      enddo !k
!
! --- set t/s in massless layers. step 3: increase salinity where water
! --- is too fresh to fit into layer k
! 
      do k=2,kk
      if (saln(i,j,k,n).lt.salmin(k)) then
        saln(i,j,k,n)=salmin(k)
        temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
      end if
      enddo !k
!
! --- redistribute temp. and salin. during both de- and entrainment
!
      thknss=dp(i,j,1,n)
      pnew=thknss+sdot(i,j)
      t1=temp(i,j,1,n)
      s1=saln(i,j,1,n)
!
      tmxl=t1+surflx(i,j)*delt1*g/(spcifh*thknss)
      smxl=s1+ vsflx(i)  *delt1*g/        thknss
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,'(i9,2i5,a,3f7.3,f8.2)') &
!diag   nstep,i,j,'  t,s,sig,dp after diab.forcing',tmxl,smxl, &
!diag   sig(tmxl,smxl),thknss*qonem
!
      if (sdot(i,j).ge.0.) then
!
! --- (mixed layer  d e e p e n s)
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,'(i9,2i5,a,f9.3,a)') &
!diag   nstep,i,j,'  entrain',sdot(i,j)*qonem,' m of water'
!
      tmxl=(tmxl*thknss+tdp(i))/pnew
      smxl=(smxl*thknss+sdp(i))/pnew
      dp(i,j,1,n)=pnew
      diaflx(i,j,1)=diaflx(i,j,1)+sdot(i,j)			! diapyc.flux
!
      else if (sdot(i,j).lt.-onecm.and.surflx(i,j).ge.0.) then  ! sdot < 0
!
! --- (mixed layer  r e c e d e s)
!
      k=klist(i,j)
      if (k.gt.kk) go to 27
!
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag  write (lp,'(i9,2i5,a,i2,a,3p,2f7.3)') nstep,i,j, &
!diag  '  sig\*(1),sig\*(',k,') =',th3d(i,j,1,n)+thbase, &
!diag  th3d(i,j,k,n)+thbase
!
      dpno=max(dp(i,j,k,n),0.)
      sn=saln(i,j,k,n)
      tn=temp(i,j,k,n)
!
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag   write (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,k, &
!diag   '  t,s,dp before detrainment',tn,sn,dpno*qonem
!
! --- distribute last time step's heating and freshwater flux over depth range
! --- 'pnew' (monin-obukhov length). split fossil mixed layer (depth= -sdot=
! --- thknss-pnew) into lower part ('lo') of depth z cooled and detrained into
! --- layer k, and an upper part ('up') heated to match temperature rise in
! --- mixed layer. transfer as much salinity as possible from sublayer 'up' to
! --- sublayer 'lo' without creating new maxima/minima in water column.
!
      dtemp=delt1*g*surflx(i,j)/(spcifh*pnew)
      dsaln=delt1*g* vsflx(i)  /        pnew
!
      tnew=    t1+dtemp
      snew=max(s1+dsaln,0.0)  !must be non-negative
!
      if (s1.le.sn .and. t1.gt.tn) then
!
! --- scenario 1: transfer t/s so as to achieve t_lo = t_k, s_lo = s_k
!
        z=-sdot(i,j)*min(1.,dtemp/max(epsil,tnew-tn))*qonem
        s_up=s1+(s1-sn)*dtemp/max(epsil*dtemp,t1-tn)
! --- is scenario 1 feasible?
        if (s_up.ge.min(snew,s1)) go to 24
      end if				! s_1 < s_n
!
! --- scenario 2: (t_lo,s_lo) differ from (tn,sn). main problem now is in
! --- maintaining density in layer k during detrainment. This requires solving
! --- 3rd deg. polynomial cc3*z**3 + cc2*z**2 + cc1*z + cc0 = 0 for z.
!
      s_up=min(s1,snew)
! --- new (t,s) in layer  k  will be t=(a*z+b)/(z+d), s=(e*z+f)/(z+d).
      a=tnew
      e=s_up
      b=(tn*dpno+    dtemp*sdot(i,j))*qonem
      f=(sn*dpno+(s_up-s1)*sdot(i,j))*qonem
      d=dpno*qonem
!
      c1msig=c1-(th3d(i,j,k,n)+thbase)
      cc0=d*d*(d*c1msig+b*c2+f*c3)+b*(d*f*c5+b*(d*c4+b*c6+f*c7))
      cc3=    (  c1msig+a*c2+e*c3)+a*(  e*c5+a*(  c4+a*c6+e*c7))
      cc1=d*(3.  *d*c1msig+(2.*b  +a*d)*c2+(2.  *f+d*e)*c3)+b*((2.*a*d &
          +b  )*c4+3.*a*b*c6+(2.*a*f+b*e)*c7)+(a*d*f+b*(d*e+  f))*c5
      cc2=  (3.  *d*c1msig+(2.*a*d+b  )*c2+(2.*d*e+  f)*c3)+a*((2.*b &
          +a*d)*c4+3.*a*b*c6+(2.*b*e+a*f)*c7)+(b  *e+a*(  f+d*e))*c5
! --- bound cc3 away from zero
      cc3=sign(max(1.e-6,abs(cc3)),cc3)
!
      x=0.0  ! dummy argument that is never used
      if (ccubq(x)**3+ccubr(x)**2.gt.0.) then
! --- one real root
      num=1
      z=root(x)
      else
! --- three real roots
      num=3
      z=root1(x)
      end if
!
!diag if (i.eq.itest.and.j.eq.jtest) then
!diag   work(1)=z
!diag   if (num.eq.3) then
!diag     work(2)=root2(x)
!diag     work(3)=root3(x)
!diag   end if
!diag   write (lp,100) nstep,i,j,' t,s,dp( 1)=',tnew,snew, &
!diag    thknss*qonem,'sdot,z=',sdot(i,j)*qonem,z,'t,s,dp(',k,')=',tn, &
!diag    sn,dpno*qonem,'real root(s):',(work(nu),nu=1,num)
!diag end if
 100  format (i9,2i5,a,2f7.3,f8.2,3x,a,2f8.2/20x,a,i2,a,2f7.3,f8.2, &
       3x,a,1p3e11.4)
!
! --- does root fall into appropriate range?
      if (z.le.0.005) go to 27
!
! --- ready to detrain lowest 'z' meters from mixed layer
!
      temp(i,j,k,n)=(a*z+b)/(z+d)
      saln(i,j,k,n)=(e*z+f)/(z+d)
!
 24   continue
      sdot(i,j)=max(sdot(i,j),-z*onem)
      dp(i,j,1,n)=thknss+sdot(i,j)
      dp(i,j,k,n) =dpno  -sdot(i,j)
      smxl=(snew*pnew+s_up*(dp(i,j,1,n)-pnew))/dp(i,j,1,n)
      tmxl=tnew
      diaflx(i,j,1)=diaflx(i,j,1)+sdot(i,j)			! diapyc.flux
      diaflx(i,j,k)=diaflx(i,j,k)-sdot(i,j)			! diapyc.flux
!
! --- inject 'ventilation' tracer into layer k
      do ktr= 1,ntracr
        tracer(i,j,k,n,ktr)=(tracer(i,j,k,n,ktr)*dpno-sdot(i,j)) &
                                     /(dpno-sdot(i,j))
      enddo !ktr
!
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag   write (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,k, &
!diag   '  t,s,dp after  detrainment',temp(i,j,k,n),saln(i,j,k,n), &
!diag   dp(i,j,k,n)*qonem
!
      end if				!  sdot > or < 0
!
 27   continue
      temp(i,j,1,n)=tmxl
      saln(i,j,1,n)=smxl
      th3d(i,j,1,n)=sig(tmxl,smxl)-thbase
      do ktr= 1,ntracr
        tracer(i,j,1,n,ktr)=1.0
      enddo !ktr
!
      dpmixl(i,j,n)=dp(  i,j,1,n)
      dpbl(  i,j)  =dp(  i,j,1,n)
      tmix(  i,j)  =temp(i,j,1,n)
      smix(  i,j)  =saln(i,j,1,n)
      thmix( i,j)  =th3d(i,j,1,n)
!
!diag if (i.eq.itest.and.j.eq.jtest) write &
!diag   (lp,'(i9,2i5,i3,a,2f9.4,f8.2)') nstep,i,j,1, &
!diag   ' final mixed-layer t,s,dp ',tmxl,smxl,dp(i,j,1,n)*qonem
!
      endif !ip
      enddo !i
      return
      end
!
!
!> Revision history:
!>
!> June 1995 - removed restriction  'klist(i,j) .le. kk'
!> June 1995 - added code for setting t/s in massless layers below mix.layer
!> Oct. 1995 - removed bug created while changing klist (June 1995 revision):
!>             'if (k.gt.kk) go to 26' now reads 'if (k.gt.kk) go to 27'
!> May  1997 - changed -sdot- into local array
!> Mar. 1998 - added -th3d-
!> Nov. 1998 - fixed bug in computing tnew,snew in situations where z < 0.005
!> Dec. 1998 - replaced dsaln by (s_up-s1) in definition of 'f'
!> Feb. 1999 - limited 'tofsig' call in loop 45 to cases where saln < salmin
!> Aug. 2000 - adapted from micom 2.8 to run within hycom 1.0
!> May  2002 - buoyfl (into the ocean), calculated here
!> Aug. 2011 - replaced dpold,dpoldm with dpo
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Aug. 2018 - added wtrflx, salflx now only actual salt flux
!> Nov. 2018 - allow for wtrflx in buoyancy flux 
