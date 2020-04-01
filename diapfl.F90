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
      subroutine diapf1(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
! --- KPP-style implicit interior diapycnal mixing
      implicit none
!
      integer m,n
!
! --------------------
! --- diapycnal mixing
! --------------------
!
! --- interior diapycnal mixing due to three processes:
! ---   shear instability
! ---   double diffusion
! ---   background internal waves
!
! --- this is essentially the k-profile-parameterization (kpp) mixing model
! --- (mxkpp.f) with all surface boundary layer processes removed
!
! --- uses the same tri-diagonal matrix solution of vertical diffusion 
! --- equation as mxkpp.f
!
      integer j
!
      if     (mod(nstep,  mixfrq).ne.0 .and. &
              mod(nstep+1,mixfrq).ne.0      ) then
        return  ! diapycnal mixing only every mixfrq,mixfrq+1 time steps
      endif
!
      call xctilr(u(1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_uv)
      call xctilr(v(1-nbdy,1-nbdy,1,n),1,kk, 1,1, halo_vv)
      call xctilr(p(1-nbdy,1-nbdy,2  ),1,kk, 1,1, halo_ps)
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call diapf1aj(m,n, j)
      enddo
!$OMP END PARALLEL DO
!
! --- momentum mixing
!
      call xctilr(vcty(1-nbdy,1-nbdy,2),1,kk, 1,1, halo_ps)
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        call diapf1bj(m,n, j)
      enddo
!$OMP END PARALLEL DO
!
      return
      end
      subroutine diapf1aj(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
      integer i
!
      do i=1,ii
        if (SEA_P) then
          call diapf1aij(m,n, i,j)
        endif !ip
      enddo !i
!
      return
      end
      subroutine diapf1bj(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
      integer i
!
      do i=1,ii
        if (SEA_U) then
          call diapf1uij(m,n, i,j)
        endif !iu
        if (SEA_V) then
          call diapf1vij(m,n, i,j)
        endif !iv
      enddo !i
!
      return
      end
      subroutine diapf1aij(m,n, i,j)
      use mod_xc         ! HYCOM communication interface     
      use mod_cb_arrays  ! HYCOM saved arrays
#if defined(STOKES)
      use mod_stokes     ! HYCOM Stokes drift
#endif
      
!
! --- hycom version 1.0
! --- KPP-style implicit interior diapycnal mixing
      implicit none
!
      integer m,n, i,j
!
! -----------------------------------------------
! --- diapycnal mixing, single i,j point (part A)
! -----------------------------------------------
!
! --- interior diapycnal mixing due to three processes:
! ---   shear instability
! ---   double diffusion
! ---   background internal waves
!
! --- this is essentially the k-profile-parameterization (kpp) mixing model
! --- (mxkpp.f) with all surface boundary layer processes removed
!
! --- uses the same tri-diagonal matrix solution of vertical diffusion 
! --- equation as mxkpp.f
!
! local variables for kpp mixing
      real shsq(kdm+1)         ! velocity shear squared
      real alfadt(kdm+1)       ! t contribution to density jump
      real betads(kdm+1)       ! s contribution to density jump
      real dbloc(kdm+1)        ! buoyancy jump across interface
      real hwide(kdm)          ! layer thicknesses in m
      real rrho                ! double diffusion parameter
      real diffdd              ! double diffusion diffusivity scale
      real prandtl             ! prandtl number
      real rigr                ! local richardson number
      real fri                 ! function of Rig for KPP shear instability
      real dflsiw              ! wave diffusivity
      real dflmiw              ! wave viscosity
      real dflbot(kdm+1)       ! bottom intensified background viscosity
!
! --- local 1-d arrays for matrix inversion
      real t1do(kdm+1),t1dn(kdm+1),s1do(kdm+1),s1dn(kdm+1), &
           tr1do(kdm+1,mxtrcr),tr1dn(kdm+1,mxtrcr), &
           difft(kdm+1),diffs(kdm+1),difftr(kdm+1), &
           zm(kdm+1),hm(kdm),dzb(kdm)
!
! --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),          & ! upper coeff for (k-1) on k line of trid.matrix
           tcc(kdm),          & ! central ...     (k  ) ..
           tcl(kdm),          & ! lower .....     (k-1) ..
           rhs(kdm)          ! right-hand-side terms
!
      real ratio,froglp,q,wq,wt,wz0,wz1,wz2
      integer k,k1,ka,kmask,ktr,nlayer,mixflg
      real riv_input
!
      real, parameter :: difriv =   50.0e-4  !river diffusion
!
# include "stmt_fns.h"
      froglp=.5*max(2,mixfrq)
!
! --- internal wave diffusion/viscosity
      dflsiw = diws(i,j)
      dflmiw = diwm(i,j)
!
! --- locate lowest substantial mass-containing layer. avoid near-zero
! --- thickness layers near the bottom
      klist(i,j)=0
      kmask=0
!
      do k=1,kk
        p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n) 
        if (dp(i,j,k,n).lt.onemm) kmask=1
        if (p(i,j,k).lt.p(i,j,kk+1)-onem.and.kmask.eq.0) klist(i,j)=k
      enddo
!
! --- dflbot (Prandtl number of one, i.e. diffusion = viscosity)
      if     (botdiw) then
! ---   diffusion coefficent profile is: K = Kb / (1 + h/h0)**2
! ---    where diwbot = Kb; diwqh0 = 1/h0 (input as h0, see forfun.f)
! ---     Decloedt T. and D.S. Luther, 2009: On a Simple Empirical
! ---     Parameterization of Topography-Catalyzed Diapycnal Mixing
! ---     in the Abyssal Ocean. JPO, 40, pp 487-508.
! ---   use Simpson's rule to estimate the average K across each layer
        wz0 = 1.0 / (1.0 + (p(i,j,kk+1)-     p(i,j,1)  )*diwqh0(i,j))**2
        wz1 = 1.0 / (1.0 + (p(i,j,kk+1)-0.5*(p(i,j,1)+ &
                                             p(i,j,2) ))*diwqh0(i,j))**2
        wz2 = 1.0 / (1.0 + (p(i,j,kk+1)-     p(i,j,2)  )*diwqh0(i,j))**2
        wq  = wz0 + wz2 + 4.0*wz1  !6 times the average value over layer 1
        do k= 2,kk
          wt  = wq
          wz0 = wz2
          wz1 = 1.0 / (1.0 + (p(i,j,kk+1) - &
                              0.5*(p(i,j,k)  + &
                                   p(i,j,k+1) ))*diwqh0(i,j))**2
          wz2 = 1.0 / (1.0 + (p(i,j,kk+1) - &
                                   p(i,j,k+1)  )*diwqh0(i,j))**2
          wq  = wz0 + wz2 + 4.0*wz1  !6 times the average value over layer k
          dflbot(k) = (0.5/6.0)*(wt+wq)*diwbot(i,j)
        enddo
        dflbot(kk+1) = dflbot(kk)
      else
        do k= 2,kk+1
          dflbot(k) = 0.0
        enddo
      endif !botdiw
!
! --- calculate vertical grid and layer widths
      do k=1,kk
        if (k.eq.1) then
          hwide(k)=dp(i,j,k,n)*qonem
          zgrid(i,j,k)=-.5*hwide(k)
        elseif (k.lt.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)*qonem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
        elseif (k.eq.klist(i,j)) then
          hwide(k)=dp(i,j,k,n)*qonem
          zgrid(i,j,k)=zgrid(i,j,k-1)-.5*(hwide(k-1)+hwide(k))
          zgrid(i,j,k+1)=zgrid(i,j,k)-.5*hwide(k)
        else
          hwide(k)=0.
        endif
      enddo
!
! --- calculate interface variables required to estimate interior diffusivities
      do k=1,kk
        k1=    k+1
        ka=min(k+1,kk)
        if (k.le.klist(i,j)) then
#if defined(STOKES)
!DAN==========================================================================
!DAN  U and V Stokes Drift Layer Average Velocities addes to Shear calculation
!DAN
      shsq(k1)=(u(i,j,k, n)+u(i+1,j,k, n)+usd(i,j,k)+usd(i+1,j,k)- &
         u(i,j,ka,n)-usd(i+1,j,ka)-u(i,j,ka,n)-usd(i+1,j,ka))**2+ &
        (v(i,j,k, n)+v(i,j+1,k, n)+vsd(i,j,k)+vsd(i,j+1,k)- &
         v(i,j,ka,n)-v(i,j+1,ka,n)-vsd(i,j,ka)-vsd(i,j+1,ka))**2
#else
          shsq(  k1)=(u(i,j,k, n)+u(i+1,j,k, n)- &
                      u(i,j,ka,n)-u(i+1,j,ka,n))**2+ &
                     (v(i,j,k, n)+v(i,j+1,k, n)- &
                      v(i,j,ka,n)-v(i,j+1,ka,n))**2
#endif
          if (locsig) then
            alfadt(k1)=dsiglocdt(ahalf*(temp(i,j,k ,n)+ &
                                        temp(i,j,ka,n) ), &
                                 ahalf*(saln(i,j,k, n)+ &
                                        saln(i,j,ka,n) ),p(i,j,k1))* &
                                       (temp(i,j,k ,n)- &
                                        temp(i,j,ka,n) )
            betads(k1)=dsiglocds(ahalf*(temp(i,j,k ,n)+ &
                                        temp(i,j,ka,n) ), &
                                 ahalf*(saln(i,j,k, n)+ &
                                        saln(i,j,ka,n) ),p(i,j,k1))* &
                                       (saln(i,j,k ,n)- &
                                        saln(i,j,ka,n) )
          else
            alfadt(k1)=dsigdt(ahalf*(temp(i,j,k ,n)+ &
                                     temp(i,j,ka,n) ), &
                              ahalf*(saln(i,j,k, n)+ &
                                     saln(i,j,ka,n) ) )* &
                                    (temp(i,j,k ,n)- &
                                     temp(i,j,ka,n) )
            betads(k1)=dsigds(ahalf*(temp(i,j,k ,n)+ &
                                     temp(i,j,ka,n) ), &
                              ahalf*(saln(i,j,k, n)+ &
                                     saln(i,j,ka,n) ) )* &
                                    (saln(i,j,k ,n)- &
                                     saln(i,j,ka,n) )
          endif
          dbloc(k1)=-g*svref*(alfadt(k1)+betads(k1))
        endif
      enddo
!
! --- determine interior diffusivity profiles throughout the water column
! --- limit mixing to the stratified interior of the ocean
!
      do k=1,kk+1
        vcty(i,j,k)=0.
        dift(i,j,k)=0.
        difs(i,j,k)=0.
      enddo
!
! --- shear instability plus background internal wave contributions
      do k=2,kk+1
        if (k-1.le.klist(i,j) .and. p(i,j,k).gt.dpmixl(i,j,n)) then
          if     (shinst) then
            q =zgrid(i,j,k-1)-zgrid(i,j,k) !0.5*(hwide(k-1)+hwide(k))
            rigr=max(0.0,dbloc(k)*q/(shsq(k)+epsil))
            ratio=min(rigr*qrinfy,1.0)
            fri=(1.0-ratio*ratio)
            fri=fri*fri*fri
            vcty(i,j,k)=difm0*fri+dflmiw+dflbot(k)
            difs(i,j,k)=difs0*fri+dflsiw+dflbot(k)
            dift(i,j,k)=difs(i,j,k)
          else
            vcty(i,j,k)=dflmiw+dflbot(k)
            difs(i,j,k)=dflsiw+dflbot(k)
            dift(i,j,k)=dflsiw+dflbot(k)
          endif
        endif
      enddo
!
! --- double-diffusion (salt fingering and diffusive convection)
      if (dbdiff) then
        do k=2,kk+1
          if (k-1.le.klist(i,j) .and. p(i,j,k).gt.dpmixl(i,j,n)) then
!
! --- salt fingering case
            if (-alfadt(k).gt.betads(k) .and. betads(k).gt.0.) then
              rrho= min(-alfadt(k)/betads(k),rrho0)
              diffdd=1.-((rrho-1.)/(rrho0-1.))**2
              diffdd=dsfmax*diffdd*diffdd*diffdd
              dift(i,j,k)=dift(i,j,k)+0.7*diffdd
              difs(i,j,k)=difs(i,j,k)+diffdd
!
! --- diffusive convection case
            elseif (alfadt(k).gt.0.0 .and. betads(k).lt.0.0 .and. &
                   -alfadt(k).gt.betads(k)) then
              rrho=-alfadt(k)/betads(k)
              diffdd=1.5e-6*9.*.101*exp(4.6*exp(-.54*(1./rrho-1.)))
              prandtl=.15*rrho
              if (rrho.gt..5) prandtl=(1.85-.85/rrho)*rrho
              dift(i,j,k)=dift(i,j,k)+diffdd
              difs(i,j,k)=difs(i,j,k)+prandtl*diffdd
            endif
          endif
        enddo
      endif
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,101) &
!diag  (nstep,i+i0,j+i0,k, &
!diag  hwide(k),1.e4*vcty(i,j,k),1.e4*dift(i,j,k), &
!diag  1.e4*difs(i,j,k),k=1,kk+1)
!
! --- perform the vertical mixing at p points
!
      mixflg=0
      do k=1,klist(i,j)
        if (dift(i,j,k+1).gt.0. .or. &
            difs(i,j,k+1).gt.0.) mixflg=mixflg+1
        difft( k+1)=froglp*dift(i,j,k+1)
        diffs( k+1)=froglp*difs(i,j,k+1)
        difftr(k+1)=froglp*difs(i,j,k+1)
        t1do(k)=temp(i,j,k,n)
        s1do(k)=saln(i,j,k,n)
        do ktr= 1,ntracr
          tr1do(k,ktr)=tracer(i,j,k,n,ktr)
        enddo
        hm(k)=hwide(k)
        zm(k)=zgrid(i,j,k)
      enddo
!
      if (mixflg.le.1) return
      nlayer=klist(i,j)
      k=nlayer+1
      ka=min(k,kk)
      difft( k)=0.
      diffs( k)=0.
      difftr(k)=0.
      t1do(k)=temp(i,j,ka,n)
      s1do(k)=saln(i,j,ka,n)
      do ktr= 1,ntracr
        tr1do(k,ktr)=tracer(i,j,ka,n,ktr)
      enddo
      zm(k)=zgrid(i,j,k)
!
! --- do rivers here because difs is also used for tracers.
      if(cpl_orivers.and.cpl_irivers) then
         riv_input = imp_orivers(i,j,1)+imp_irivers(i,j,1)
      else
         riv_input = rivers(i,j,1)
      endif
      riv_input = rivers(i,j,1)
      if     (thkriv.gt.0.0 .and. riv_input.ne.0.0) then
        do k=1,nlayer
          if     (-zm(k)+0.5*hm(k).lt.thkriv) then !interface<thkriv
            diffs(k+1) = max(diffs(k+1),froglp*difriv)
          endif
        enddo !k
      endif !river
!
! --- compute factors for coefficients of tridiagonal matrix elements.
!       tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
!       tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
!
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
!
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
!
! --- solve the diffusion equation
!
! --- t solution
      call tridcof(difft,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=t1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,t1do,t1dn,difft, i,j)
      if     ( tofset.eq.0.0 .or. &
              (mod(nstep  ,tsofrq).ne.0 .and. &
               mod(nstep+1,tsofrq).ne.0      ) ) then
        do k=1,nlayer
          temp(i,j,k,n)=t1dn(k)
        enddo
      else  !include tofset drift correction
        do k=1,nlayer
          temp(i,j,k,n)=t1dn(k) + baclin*max(2,tsofrq)*tofset
        enddo
      endif !without:with tofset
!
! --- t-like tracer solution
      do ktr= 1,ntracr
        if     (trcflg(ktr).eq.2) then
          do k=1,nlayer
            rhs(k)=tr1do(k,ktr)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,tr1do,tr1dn,difft, i,j)
          do k=1,nlayer
            tracer(i,j,k,n,ktr)=tr1dn(k,ktr)
          enddo
        endif
      enddo !ktr
!
! --- s solution and th3d reset
      call tridcof(diffs,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=s1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,s1do,s1dn,diffs, i,j)
      if     ( sofset.eq.0.0 .or. &
              (mod(nstep  ,tsofrq).ne.0 .and. &
               mod(nstep+1,tsofrq).ne.0      ) ) then
        do k=1,nlayer
          saln(i,j,k,n)=s1dn(k)
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        enddo
      else  !include sofset drift correction
        do k=1,nlayer
          saln(i,j,k,n)=s1dn(k) + baclin*max(2,tsofrq)*sofset
          th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
        enddo
      endif !without:with sofset
!
! --- standard tracer solution
      if     (ntracr.gt.0) then
        call tridcof(difftr,tri,nlayer,tcu,tcc,tcl)
      endif
      do ktr= 1,ntracr
        if     (trcflg(ktr).ne.2) then
          do k=1,nlayer
            rhs(k)=tr1do(k,ktr)
          enddo
          call tridmat(tcu,tcc,tcl,nlayer, &
                       hm,rhs,tr1do(1,ktr),tr1dn(1,ktr),difftr, i,j)
          do k=1,nlayer
            tracer(i,j,k,n,ktr)=tr1dn(k,ktr)
          enddo
        endif
      enddo !ktr
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,102) &
!diag  (nstep,i+i0,j+j0,k, &
!diag difft(k),t1do(k),t1dn(k),t1dn(k)-t1do(k), &
!diag diffs(k),s1do(k),s1dn(k),s1dn(k)-s1do(k),k=1,nlayer)
!
      return
!
 101  format(25x,'   thick      viscty    t diff    s diff  ' &
           /(i9,2i5,i3,2x,4f10.2))
 102  format(25x, &
      ' diff t  t old   t new   t chng  diff s  s old   s new   s chng' &
           /(i9,2i5,i3,1x,8f8.3))
      end
      subroutine diapf1uij(m,n, i,j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
! --- KPP-style implicit interior diapycnal mixing
      implicit none
!
      integer m,n, i,j 
!
! -----------------------------------------------------------------
! --- diapycnal mixing, single i,j point, momentum at u grid points
! -----------------------------------------------------------------
!
! --- local 1-d arrays for matrix inversion
      real u1do(kdm+1),u1dn(kdm+1), &
           diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
!
! --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),          & ! upper coeff for (k-1) on k line of trid.matrix
           tcc(kdm),          & ! central ...     (k  ) ..
           tcl(kdm),          & ! lower .....     (k-1) ..
           rhs(kdm)          ! right-hand-side terms
!
      real presu,froglp
      integer k,ka,kmask(idm),nlayer,mixflg
!
      froglp=.5*max(2,mixfrq)
!
      presu=0.
      kmask(1)=0
      mixflg=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpu(i,j,ka,n).lt.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presu.lt.depthu(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*froglp*(vcty(i,j,k+1)+vcty(i-1,j,k+1))
          if (diffm(k+1).gt.0.) mixflg=mixflg+1
          u1do(k)=u(i,j,k,n)
          hm(k)=dpu(i,j,k,n)*qonem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presu=presu+dpu(i,j,k,n)
          nlayer=k
        elseif (k.eq.nlayer+1) then
          diffm(k)=0.
          u1do(k)=u1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
      if (mixflg.le.1) return
!
! --- compute factors for coefficients of tridiagonal matrix elements.
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
!
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
!
! --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)= u1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,u1do,u1dn,diffm, i,j)
      do k=1,nlayer
        u(i,j,k,n)=u1dn(k)
      enddo
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,106) &
!diag  (nstep,i+i0,j+j0,k,hm(k),u1do(k),u1dn(k),k=1,nlayer)
!
      return
 106  format(23x,'   thick   u old   u new'/(i9,2i5,i3,1x,f10.3,2f8.3))
      end
      subroutine diapf1vij(m,n, i,j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
! --- KPP-style implicit interior diapycnal mixing
      implicit none
!
      integer m,n, i,j 
!
! -----------------------------------------------------------------
! --- diapycnal mixing, single i,j point, momentum at v grid points
! -----------------------------------------------------------------
!
! --- local 1-d arrays for matrix inversion
      real v1do(kdm+1),v1dn(kdm+1), &
           diffm(kdm+1),zm(kdm+1),hm(kdm),dzb(kdm)
!
! --- tridiagonal matrix solution arrays
      real tri(kdm,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kdm),          & ! upper coeff for (k-1) on k line of trid.matrix
           tcc(kdm),          & ! central ...     (k  ) ..
           tcl(kdm),          & ! lower .....     (k-1) ..
           rhs(kdm)          ! right-hand-side terms
!
      real presv,froglp
      integer k,ka,kmask(idm),nlayer,mixflg
!
      froglp=.5*max(2,mixfrq)
!
      presv=0.
      kmask(1)=0
      mixflg=0
      do k=1,kk+1
        ka=min(k,kk)
        if (dpv(i,j,ka,n).lt.tencm.or.k.eq.kk+1) kmask(1)=1
        if (presv.lt.depthv(i,j)-tencm.and.kmask(1).eq.0) then
          diffm(k+1)=.5*froglp*(vcty(i,j,k+1)+vcty(i,j-1,k+1))
          if (diffm(k+1).gt.0.) mixflg=mixflg+1
          v1do(k)=v(i,j,k,n)
          hm(k)=dpv(i,j,k,n)*qonem
          if (k.eq.1) then
            zm(k)=-.5*hm(k)
          else
            zm(k)=zm(k-1)-.5*(hm(k-1)+hm(k))
          endif
          presv=presv+dpv(i,j,k,n)
          nlayer=k
        elseif (k.eq.nlayer+1) then
          diffm(k)=0.
          v1do(k)=v1do(k-1)
          zm(k)=zm(k-1)-.5*hm(k-1)
        endif
      enddo
      if (mixflg.le.1) return
!
! --- compute factors for coefficients of tridiagonal matrix elements.
!
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
!
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
!
! --- solve the diffusion equation
      call tridcof(diffm,tri,nlayer,tcu,tcc,tcl)
      do k=1,nlayer
        rhs(k)=v1do(k)
      enddo
      call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,v1do,v1dn,diffm, i,j)
      do k=1,nlayer
        v(i,j,k,n)=v1dn(k)
      enddo
!
!diag if (i.eq.itest.and.j.eq.jtest) write (lp,107) &
!diag  (nstep,i+i0,j+j0,k,hm(k),v1do(k),v1dn(k),k=1,nlayer)
!
      return
 107  format(23x,'   thick   v old   v new'/(i9,2i5,i3,1x,f10.3,2f8.3))
      end
!
      subroutine diapf2(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
! --- MICOM-style explict interior diapycnal mixing for hybrid coordinates
      implicit none
!
      integer m,n
!
      integer j
!
      if (diapyc.eq.0. .or. (mod(nstep  ,mixfrq).ne.0 .and. &
                             mod(nstep+1,mixfrq).ne.0)) return
!diag write (lp,'(i9,3x,a)') nstep,'entering   d i a p f l'
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do 31 j=1,jj
        call diapf2j(m,n, j)
   31 continue
!$OMP END PARALLEL DO
!
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                  dpv(1-nbdy,1-nbdy,1,n), &
                  p,depthu,depthv, 0,0)
!
!diag write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
      subroutine diapf2j(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n, j
!
      integer i,k,k1,k2,ka,kmin(idm),kmax(idm),ktr
      real flxu(idm,kdm),flxl(idm,kdm),pdot(idm,kdm),flngth(idm,kdm), &
           ennsq,alfa,beta,q,qmin,qmax,amount,froglp,delp, &
           alfadt1,alfadt2,betads1,betads2,plev, &
           trflxu(idm,0:kdm+1,mxtrcr), &
           trflxl(idm,0:kdm+1,mxtrcr),cliptr(idm,mxtrcr), &
            tflxu(idm,0:kdm+1), tflxl(idm,0:kdm+1),clipt( idm), &
            sflxu(idm,0:kdm+1), sflxl(idm,0:kdm+1),clips( idm), &
           told(idm,2),sold(idm,2),trold(idm,2,mxtrcr)
!     real totem(idm),tosal(idm),tndcyt,tndcys	! col.integrals (diag.use only)
!
      real       small
      parameter (small=1.e-6)
!
# include "stmt_fns.h"
!
! --- -------------------------------
! --- diapycnal mixing, single j-row.
! --- -------------------------------
!
! --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
!
      do i=1,ii
      if (SEA_P) then
!
! --- t/s conservation diagnostics (optional):
!       totem(i)=0.
!       tosal(i)=0.
!       do k=1,kk
!         totem(i)=totem(i)+temp(i,j,k,n)*dp(i,j,k,n)
!         tosal(i)=tosal(i)+saln(i,j,k,n)*dp(i,j,k,n)
!       enddo
!
      do k=1,kk
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo !k
!
      sold(i,1)=saln(i,j,kk,n)
      told(i,1)=temp(i,j,kk,n)
      tflxl(i,   0)=0.
      tflxu(i,kk+1)=0.
      sflxl(i,   0)=0.
      sflxu(i,kk+1)=0.
      do ktr= 1,ntracr
        trold( i,   1,ktr)=tracer(i,j,kk,n,ktr)
        trflxl(i,   0,ktr)=0.
        trflxu(i,kk+1,ktr)=0.
      enddo !ktr
!
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag   write (lp,'(i9,2i5,3x,a/(i36,4f10.3))') nstep,i+i0,j+j0, &
!diag   'before diapf2: thickness  salinity temperature density', &
!diag   (k,dp(i,j,k,n)*qonem,saln(i,j,k,n), &
!diag   temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=1,kk)
!
      kmin(i)=kk+1
      kmax(i)=1
!
      do k=2,kk
!
! --- locate lowest mass-containing layer and upper edge of stratified region
      if (p(i,j,k).lt.p(i,j,kk+1)-onemm)  then
        kmax(i)=k
        if (kmin(i).eq.kk+1 .and. &
            th3d(i,j,k,n).gt.th3d(i,j,k-1,n)+sigjmp) then
          kmin(i)=k
        endif
      endif
      enddo !k
!
!diag if (j.eq.jtest.and.itest.ge.ifp(j,l).and.itest.le.ilp(j,l)) &
!diag   write (lp,'(i9,2i5,a,2i5)') &
!diag     nstep,itest+i0,j+j0,' kmin,kmax =',kmin(itest),kmax(itest)
!
! --- find buoyancy frequency for each layer
!
      do k=2,kk-1
      k1=k-1
      k2=k+1
!
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
! --- ennsq = buoy.freq.^2 / g^2
        if (locsig) then
          alfadt1=dsiglocdt(ahalf*(temp(i,j,k1,n)+ &
                                   temp(i,j,k ,n) ), &
                            ahalf*(saln(i,j,k1,n)+ &
                                   saln(i,j,k, n) ),p(i,j,k))* &
                                  (temp(i,j,k1,n)- &
                                   temp(i,j,k, n) )
          betads1=dsiglocds(ahalf*(temp(i,j,k1,n)+ &
                                   temp(i,j,k ,n) ), &
                            ahalf*(saln(i,j,k1,n)+ &
                                   saln(i,j,k, n) ),p(i,j,k))* &
                                  (saln(i,j,k1,n)- &
                                   saln(i,j,k, n) )
          alfadt2=dsiglocdt(ahalf*(temp(i,j,k ,n)+ &
                                   temp(i,j,k2,n) ), &
                            ahalf*(saln(i,j,k ,n)+ &
                                   saln(i,j,k2,n) ),p(i,j,k2))* &
                                  (temp(i,j,k, n)- &
                                   temp(i,j,k2,n) )
          betads2=dsiglocdt(ahalf*(temp(i,j,k ,n)+ &
                                   temp(i,j,k2,n) ), &
                            ahalf*(saln(i,j,k ,n)+ &
                                   saln(i,j,k2,n) ),p(i,j,k2))* &
                                  (saln(i,j,k, n)- &
                                   saln(i,j,k2,n) )
          ennsq=-min(0.,min(alfadt1+betads1,alfadt2+betads2)) &
               /max(p(i,j,k2)-p(i,j,k),onem)
        else
          ennsq=max(0.,min(th3d(i,j,k2,n)-th3d(i,j,k ,n), &
                           th3d(i,j,k ,n)-th3d(i,j,k1,n))) &
               /max(p(i,j,k2)-p(i,j,k),onem)
        endif
! --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
! --- (dimensions of flngth: length in pressure units)
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(i,k)=diapyc*sqrt(ennsq) * baclin*froglp * onem
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc
!cc        flngth(i,k)=diapyc*ennsq*g * baclin*froglp * onem
! -----------------------------------------------------------------------
!
      endif
      enddo !k
!
! --- find t/s fluxes at the upper and lower interface of each layer
! --- (compute only the part common to t and s fluxes)
!
      do k=1,kk
      flxu(i,k)=0.
      flxl(i,k)=0.
!
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
!
        if (locsig) then
          plev=p(i,j,k)+0.5*dp(i,j,k,n)
          alfa=-svref*dsiglocdt(temp(i,j,k,n),saln(i,j,k,n),plev)
          beta= svref*dsiglocds(temp(i,j,k,n),saln(i,j,k,n),plev)
        else
          alfa=-svref*dsigdt(temp(i,j,k,n),saln(i,j,k,n))
          beta= svref*dsigds(temp(i,j,k,n),saln(i,j,k,n))
        endif
!
        flxu(i,k)=flngth(i,k)/ &
          max(beta*(saln(i,j,k,n)-saln(i,j,k-1,n)) &
             -alfa*(temp(i,j,k,n)-temp(i,j,k-1,n)),small)
        flxl(i,k)=flngth(i,k)/ &
          max(beta*(saln(i,j,k+1,n)-saln(i,j,k,n)) &
             -alfa*(temp(i,j,k+1,n)-temp(i,j,k,n)),small)
!
        q=min(1.,.5*min(p(i,j,k)-p(i,j,k-1),p(i,j,k+2)-p(i,j,k+1))/ &
          max(flxu(i,k),flxl(i,k),epsil))
!
!diag   if (q.ne.1.) write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f5.2)')  &
!diag     nstep,i+i0,j+j0,k,' flxu/l,dpu/l,q=',flxu(i,k),flxl(i,k), &
!diag     (p(i,j,k)-p(i,j,k-1))*qonem,(p(i,j,k+2)-p(i,j,k+1))*qonem,q
!
        flxu(i,k)=flxu(i,k)*q
        flxl(i,k)=flxl(i,k)*q
!
      endif				!  kmin < k < kmax
!
!diag if (i.eq.itest.and.j.eq.jtest.and.k.ge.kmin(i).and.k.le.kmax(i)) &
!diag    write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)') &
!diag    nstep,i+i0,j+j0,k, &
!diag    'thknss   temp   saln    flngth      flxu      flxl', &
!diag    dp(i,j,k,n)*qonem,temp(i,j,k,n),saln(i,j,k,n),flngth(i,k), &
!diag    flxu(i,k)*qonem,flxl(i,k)*qonem
!
      enddo !k
!
! --- determine mass flux -pdot- implied by t/s fluxes.
!
      do k=1,kk
      if (k.gt.kmin(i) .and. k.le.kmax(i)) then
          pdot(i,k)=flxu(i,k)-flxl(i,k-1)
      else
          pdot(i,k)=0.
      endif
      enddo !k
!
! --- convert flxu,flxl into actual t/s (and tracer) fluxes
!
      do k=1,kk
      tflxu(i,k)=0.
      tflxl(i,k)=0.
      sflxu(i,k)=0.
      sflxl(i,k)=0.
      if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
        tflxu(i,k)=flxu(i,k)*temp(i,j,k-1,n)
        sflxu(i,k)=flxu(i,k)*saln(i,j,k-1,n)
!
        tflxl(i,k)=flxl(i,k)*temp(i,j,k+1,n)
        sflxl(i,k)=flxl(i,k)*saln(i,j,k+1,n)
      endif
      do ktr= 1,ntracr
        trflxu(i,k,ktr)=0.
        trflxl(i,k,ktr)=0.
        if (k.gt.kmin(i) .and. k.lt.kmax(i)) then
          trflxu(i,k,ktr)=flxu(i,k)*tracer(i,j,k-1,n,ktr)
          trflxl(i,k,ktr)=flxl(i,k)*tracer(i,j,k+1,n,ktr)
        endif
      enddo !ktr
      enddo !k
!
      do ktr= 1,ntracr
        cliptr(i,ktr)=0.
      enddo !ktr
      clipt( i)=0.
      clips( i)=0.
!
! --- update interface pressure and layer temperature/salinity
      do k=kk,1,-1
      ka=max(1,k-1)
!
      sold(i,2)=sold(i,1)
      sold(i,1)=saln(i,j,k,n)
      told(i,2)=told(i,1)
      told(i,1)=temp(i,j,k,n)
      do ktr= 1,ntracr
        trold(i,2,ktr)=trold( i,1,    ktr)
        trold(i,1,ktr)=tracer(i,j,k,n,ktr)
      enddo
!
      dpo(i,j,k,n)=dp(i,j,k,n)
      p(i,j,k)=p(i,j,k)-pdot(i,k)
      dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
!
      if (k.ge.kmin(i) .and. k.le.kmax(i)) then
        delp=dp(i,j,k,n)
        if (delp.gt.0.) then
          amount=temp(i,j,k,n)*dpo(i,j,k,n) &
            -(tflxu(i,k+1)-tflxu(i,k)+tflxl(i,k-1)-tflxl(i,k))
          q=amount
          qmax=max(temp(i,j,ka,n),told(i,1),told(i,2))
          qmin=min(temp(i,j,ka,n),told(i,1),told(i,2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clipt(i)=clipt(i)+(q-amount)
          temp(i,j,k,n)=amount/delp
!
          amount=saln(i,j,k,n)*dpo(i,j,k,n) &
            -(sflxu(i,k+1)-sflxu(i,k)+sflxl(i,k-1)-sflxl(i,k))
          q=amount
          qmax=max(saln(i,j,ka,n),sold(i,1),sold(i,2))
          qmin=min(saln(i,j,ka,n),sold(i,1),sold(i,2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clips(i)=clips(i)+(q-amount)
          saln(i,j,k,n)=amount/delp
!
          do ktr= 1,ntracr
            amount=tracer(i,j,k,n,ktr)*dpo(i,j,k,n) &
                 -(trflxu(i,k+1,ktr)-trflxu(i,k,ktr)+ &
                   trflxl(i,k-1,ktr)-trflxl(i,k,ktr))
            q=amount
            qmax=max(tracer(i,j,ka,n,ktr),trold(i,1,ktr), &
                                          trold(i,2,ktr))
            qmin=min(tracer(i,j,ka,n,ktr),trold(i,1,ktr), &
                                          trold(i,2,ktr))
            amount=max(qmin*delp,min(amount,qmax*delp))
            cliptr(i,ktr)=cliptr(i,ktr)+(q-amount)
            tracer(i,j,k,n,ktr)=amount/delp
          enddo !ktr
        endif
      endif
      enddo !k
!
      clipt(i)=clipt(i)/pbot(i,j) + baclin*froglp*tofset
      clips(i)=clips(i)/pbot(i,j) + baclin*froglp*sofset
      do ktr= 1,ntracr
        cliptr(i,ktr)=cliptr(i,ktr)/pbot(i,j)
      enddo !ktr
!
      do k=1,kk
!
! --- restore 'clipped' and 'offset' t/s amount to column
      temp(i,j,k,n)=temp(i,j,k,n)+clipt(i)
      saln(i,j,k,n)=saln(i,j,k,n)+clips(i)
      th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
      do ktr= 1,ntracr
        tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)+cliptr(i,ktr)
      enddo !ktr
!
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,k,n)-dpo(i,j,k,n))	! diapyc.flx.
! --- make sure p is computed from dp, not the other way around (roundoff!)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo !k
!
! --- t/s conservation diagnostics (optional):
!       tndcyt=-totem(i)
!       tndcys=-tosal(i)
!       do k=1,kk
!         tndcyt=tndcyt+temp(i,j,k,n)*dp(i,j,k,n)
!         tndcys=tndcys+saln(i,j,k,n)*dp(i,j,k,n)
!       enddo
!       if (abs(tndcyt/totem(i)).gt.1.e-11)
!    .  write (lp,100) i,j,'  diapf2 temp.col.intgl.:',totem(i),tndcyt,
!    .  clipt(i)
!       if (abs(tndcys/tosal(i)).gt.1.e-11)
!    .  write (lp,100)
!    .  i+i0,j+i0,'  diapf2 saln.col.intgl.:',tosal(i),tndcys,clips(i)
!100    format(2i5,a,1p,e16.8,2e13.5)
!
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag   write (lp,'(i9,2i5,3x,a/(i36,0p,4f10.3))') &
!diag   nstep,i+i0,j+j0, &
!diag   'after  diapf2: thickness  salinity temperature density', &
!diag   (k,dp(i,j,k,n)*qonem,saln(i,j,k,n), &
!diag   temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=1,kk)
!
      endif !ip
      enddo !i
!
      return
      end
!
      subroutine diapf3(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0 (adapted from micom version 2.8)
! --- MICOM-style explict interior diapycnal mixing for isopycnal coords
      implicit none
!
      integer m,n
!
      integer j
!
      if (diapyc.eq.0. .or. (mod(nstep  ,mixfrq).ne.0 .and. &
                             mod(nstep+1,mixfrq).ne.0)) return
!diag write (lp,'(i9,3x,a)') nstep,'entering   d i a p f l'
!
!$OMP PARALLEL DO PRIVATE(j) &
!$OMP              SHARED(m,n) &
!$OMP          SCHEDULE(STATIC,jblk)
      do 31 j=1,jj
        call diapf3j(m,n, j)
   31 continue
!
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                  dpv(1-nbdy,1-nbdy,1,n), &
                  p,depthu,depthv, 0,0)
!
!diag write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
      subroutine diapf3j(m,n, j)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0 (adapted from micom version 2.8)
! --- MICOM-style explict interior diapycnal mixing for isopycnic coordinates
      implicit none
!
      integer m,n,j
!
      integer i,k,k1,k2,ktr,l
      integer kmin(idm),kmax(idm)
      real flxu(idm,kdm),flxl(idm,kdm),pdot(idm,kdm),flngth(idm,kdm), &
           ennsq,alfa,beta,smax,smin,sold(idm,2),q,salt,froglp, &
           alfadt1,alfadt2,betads1,betads2,plev, &
           flxtru(idm,kdm,mxtrcr), &
           flxtrl(idm,kdm,mxtrcr),trold(idm,2,mxtrcr),trmax,trmin
!
      real       small
      parameter (small=1.e-6)
!
# include "stmt_fns.h" 
!
! --- ------------------------------
! --- diapycnal mixing, single j-row.
! --- -------------------------------
!
! --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
!
!cc      salt=0.
!
      do i=1,ii
      if (SEA_P) then
!
      do k=1,kk
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo !k
!
      sold(i,1)=saln(i,j,kk,n)
      do ktr= 1,ntracr
        trold(i,1,ktr)=tracer(i,j,kk,n,ktr)
      enddo !ktr
      flxl(i,1)=0.
!
!diag if (i.eq.itest .and. j.eq.jtest) &
!diag   write (lp,'(i9,2i5,3x,a/(i36,0p,3f10.3,3p,f10.3))') &
!diag   nstep,i+i0,j+j0, &
!diag   'before diapfl: thickness  salinity temperature density', &
!diag   1,dp(i,j,1,n)*qonem,saln(i,j,1,n),temp(i,j,1,n), &
!diag   th3d(i,j,1,n)+thbase,(k,dp(i,j,k,n)*qonem,saln(i,j,k,n), &
!diag   temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=2,kk)
!
      kmin(i)=kk+1
      kmax(i)=1
!
      do k=2,kk
!
! --- locate lowest mass-containing layer
      if (p(i,j,k).lt.p(i,j,kk+1))  then
        kmax(i)=k
!
! --- make sure salinity is compatible with density in layer k
!cc     q=sofsig(th3d(i,j,k,n)+thbase,temp(i,j,k,n))
!cc     salt=salt+(q-saln(i,j,k,n))*dp(i,j,k,n)*scp2(i)
!cc     saln(i,j,k,n)=q
!
! --- locate uppermost isopycnic layer heavier than mixed layer
        if (kmin(i).eq.kk+1 .and. &
          max(th3d(i,j,1,m),th3d(i,j,1,n))+sigjmp.le.th3d(i,j,k,n)) &
          kmin(i)=k
      end if
      enddo !k
!
!diag if (j.eq.jtest.and.itest.ge.ifp(j,l).and.itest.le.ilp(j,l)) &
!diag   write (lp,'(i9,2i5,a,2i5)') &
!diag     nstep,itest+i0,j+j0,' kmin,kmax =',kmin(itest),kmax(itest)
!
! --- temporarily swap layers  1  and  kmin-1
!
      do k=3,kk
      if (k.eq.kmin(i)) then
        q=temp(i,j,k-1,n)
        temp(i,j,k-1,n)=temp(i,j,1,n)
        temp(i,j,1,n)=q
!
        q=saln(i,j,k-1,n)
        saln(i,j,k-1,n)=saln(i,j,1,n)
        saln(i,j,1,n)=q
!
        q=th3d(i,j,k-1,n)
        th3d(i,j,k-1,n)=th3d(i,j,1,n)
        th3d(i,j,1,n)=q
!
        do ktr= 1,ntracr
          q=tracer(i,j,k-1,n,ktr)
          tracer(i,j,k-1,n,ktr)=tracer(i,j,1,n,ktr)
          tracer(i,j,  1,n,ktr)=q
        enddo
!
        flxl(i,k-1)=0.
      end if
      enddo !k
!
! --- find buoyancy frequency for each layer
!
      do k=2,kk
      k1=k-1
      k2=k+1
!
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
! --- ennsq = buoy.freq.^2 / g^2
        if (locsig) then
          alfadt1=dsiglocdt(ahalf*(temp(i,j,k1,n)+ &
                                   temp(i,j,k ,n) ), &
                            ahalf*(saln(i,j,k1,n)+ &
                                   saln(i,j,k, n) ),p(i,j,k))* &
                                  (temp(i,j,k1,n)- &
                                   temp(i,j,k, n) )
          betads1=dsiglocds(ahalf*(temp(i,j,k1,n)+ &
                                   temp(i,j,k ,n) ), &
                            ahalf*(saln(i,j,k1,n)+ &
                                   saln(i,j,k, n) ),p(i,j,k))* &
                                  (saln(i,j,k1,n)- &
                                   saln(i,j,k, n) )
          alfadt2=dsiglocdt(ahalf*(temp(i,j,k ,n)+ &
                                   temp(i,j,k2,n) ), &
                            ahalf*(saln(i,j,k ,n)+ &
                                   saln(i,j,k2,n) ),p(i,j,k2))* &
                                  (temp(i,j,k, n)- &
                                   temp(i,j,k2,n) )
          betads2=dsiglocdt(ahalf*(temp(i,j,k ,n)+ &
                                   temp(i,j,k2,n) ), &
                            ahalf*(saln(i,j,k ,n)+ &
                                   saln(i,j,k2,n) ),p(i,j,k2))* &
                                  (saln(i,j,k, n)- &
                                   saln(i,j,k2,n) )
          ennsq=-min(0.,min(alfadt1+betads1,alfadt2+betads2)) &
               /max(p(i,j,k2)-p(i,j,k),onem)
        else
          ennsq=max(0.,min(th3d(i,j,k2,n)-th3d(i,j,k  ,n), &
                           th3d(i,j,k ,n)-th3d(i,j,k1,n))) &
                /max(p(i,j,k2)-p(i,j,k),onem)
        endif
! --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc / buoyancy frequency
        flngth(i,k)=(diapyc*sqrt(ennsq) * baclin*froglp * onem)
! -----------------------------------------------------------------------
! --- use the following if exch.coeff. = diapyc
!cc        flngth(i,k)=diapyc*ennsq*g * baclin*froglp * onem
! -----------------------------------------------------------------------
!
      end if
      enddo !k
!
! --- find t/s fluxes at the upper and lower interface of each layer
! --- (compute only the part common to t and s fluxes)
!
      do k=2,kk
      flxu(i,k)=0.
      flxl(i,k)=0.
!
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
!
        if (locsig) then
          plev=p(i,j,k)+0.5*dp(i,j,k,n)
          alfa=-svref*dsiglocdt(temp(i,j,k,n),saln(i,j,k,n),plev)
          beta= svref*dsiglocds(temp(i,j,k,n),saln(i,j,k,n),plev)
        else
          alfa=-svref*dsigdt(temp(i,j,k,n),saln(i,j,k,n))
          beta= svref*dsigds(temp(i,j,k,n),saln(i,j,k,n))
        endif
!
        flxu(i,k)=flngth(i,k)/ &
          max(beta*(saln(i,j,k,n)-saln(i,j,k-1,n)) &
             -alfa*(temp(i,j,k,n)-temp(i,j,k-1,n)),small)
        flxl(i,k)=flngth(i,k)/ &
          max(beta*(saln(i,j,k+1,n)-saln(i,j,k,n)) &
             -alfa*(temp(i,j,k+1,n)-temp(i,j,k,n)),small)
!
        q=min(1.,min(p(i,j,k)-p(i,j,k-1),p(i,j,k+2)-p(i,j,k+1))/ &
          max(flxu(i,k),flxl(i,k),epsil))
!
!diag   if (i.eq.itest .and. j.eq.jtest .and. q.ne.1.) &
!diag     write (lp,'(i9,2i5,i3,a,1p,2e10.2,0p,2f7.2,f4.2)')  &
!diag     nstep,i+i0,j+j0,k,' flxu/l,dpu/l,q=',flxu(i,k),flxl(i,k), &
!diag     (p(i,j,k)-p(i,j,k-1))*qonem,(p(i,j,k+2)-p(i,j,k+1))*qonem,q
!
        flxu(i,k)=flxu(i,k)*q
        flxl(i,k)=flxl(i,k)*q
!
      end if				!  kmin < k < kmax
!
!diag if (i.eq.itest.and.j.eq.jtest.and.k.ge.kmin(i).and.k.le.kmax(i)) &
!diag    write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)') &
!diag    nstep,i+i0,j+j0,k, &
!diag    ' thknss   temp   saln   flngth      flxu      flxl', &
!diag    dp(i,j,k,n)*qonem,temp(i,j,k,n),saln(i,j,k,n),flngth(i,k), &
!diag    flxu(i,k)*qonem,flxl(i,k)*qonem
!
      enddo !k
!
! --- determine mass flux -pdot- implied by t/s fluxes.
!
      do k=2,kk
      if (k.ge.kmin(i) .and. k.le.kmax(i)) &
          pdot(i,k)=flxu(i,k)-flxl(i,k-1)
      enddo !k
!
! --- convert flxu,flxl into actual salt fluxes - calculate tracer fluxes
!
      do k=2,kk
      if (k.ge.kmin(i) .and. k.lt.kmax(i)) then
        do ktr= 1,ntracr
          flxtru(i,k,ktr)=-flxu(i,k)*(tracer(i,j,k  ,n,ktr)- &
                                      tracer(i,j,k-1,n,ktr))
          flxtrl(i,k,ktr)=-flxl(i,k)*(tracer(i,j,k+1,n,ktr)- &
                                      tracer(i,j,k  ,n,ktr))
        enddo
        flxu(i,k)=-flxu(i,k)*(saln(i,j,k  ,n)-saln(i,j,k-1,n))
        flxl(i,k)=-flxl(i,k)*(saln(i,j,k+1,n)-saln(i,j,k  ,n))
      endif
      enddo !k
!
! --- update interface pressure and layer salinity
      do k=kk,2,-1
      sold(i,2)=sold(i,1)
      sold(i,1)=saln(i,j,k,n)
      do ktr= 1,ntracr
        trold(i,2,ktr)=trold( i,1,    ktr)
        trold(i,1,ktr)=tracer(i,j,k,n,ktr)
      enddo
!
      if (k.ge.kmin(i) .and. k.le.kmax(i)) then
        p(i,j,k)=p(i,j,k)-pdot(i,k)
        saln(i,j,k,n)=saln(i,j,k,n)-(flxl(i,k)-flxu(i,k)) &
            /max(p(i,j,k+1)-p(i,j,k),small)
        do ktr= 1,ntracr
          tracer(i,j,k,n,ktr)=tracer(i,j,k,n,ktr)- &
            (flxtrl(i,k,ktr)-flxtru(i,k,ktr)) &
              /max(p(i,j,k+1)-p(i,j,k),small)
        enddo
!
! --- avoid excessive salinity (and tracer) values in totally eroded layers
        smin=min(saln(i,j,k-1,n),sold(i,1),sold(i,2))
        smax=max(saln(i,j,k-1,n),sold(i,1),sold(i,2))
        saln(i,j,k,n)=max(smin,min(saln(i,j,k,n),smax))
!
        do ktr= 1,ntracr
          trmin=min(tracer(i,j,k-1,n,ktr),trold(i,1,ktr), &
                                          trold(i,2,ktr))
          trmax=max(tracer(i,j,k-1,n,ktr),trold(i,1,ktr), &
                                          trold(i,2,ktr))
          tracer(i,j,k,n,ktr)=max(trmin,min(tracer(i,j,k,n,ktr),trmax))
        enddo
!
        temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase,saln(i,j,k,n))
      else if (kmin(i).ne.kk+1 .and. k.lt.kmin(i)) then
        p(i,j,k)=min(p(i,j,k),p(i,j,k+1))
      end if
      enddo !k
!
! --- undo effect of loop 40 (i.e., restore layers  1  and  kmin-1)
!
      do k=3,kk
      if (k.eq.kmin(i)) then
        q=temp(i,j,k-1,n)
        temp(i,j,k-1,n)=temp(i,j,1,n)
        temp(i,j,1,n)=q
!
        q=saln(i,j,k-1,n)
        saln(i,j,k-1,n)=saln(i,j,1,n)
        saln(i,j,1,n)=q
!
        q=th3d(i,j,k-1,n)
        th3d(i,j,k-1,n)=th3d(i,j,1,n)
        th3d(i,j,1,n)=q
!
        do ktr= 1,ntracr
          q=tracer(i,j,k-1,n,ktr)
          tracer(i,j,k-1,n,ktr)=tracer(i,j,1,n,ktr)
          tracer(i,j,  1,n,ktr)=q
        enddo
      endif
      enddo !k
!
      do k=kk,2,-1
!
! --- if layer 1 has been totally eroded, transfer layer -kmin- to layer 1
      if (k.eq.kmin(i) .and. p(i,j,k).lt..1*onemm) then
!
!diag   if (i.eq.itest .and. j.eq.jtest) then
!diag   write (lp,'(i9,2i5,3x,a,i3,a)') nstep,i,j,'diapfl -- layer',k, &
!diag   ' erodes mixed layer'
!diag   write (lp,'(i9,2i5,i3,3x,a/22x,f9.3,2f7.3,1p,3e10.3)') &
!diag    nstep,i+i0,j+j0,k, &
!diag    ' thknss   temp   saln   flngth      flxu      flxl', &
!diag    (p(i,j,k+1)-p(i,j,k))*qonem,temp(i,j,k,n), &
!diag    saln(i,j,k,n),flngth(i,k),flxu(i,k),flxl(i,k)
!diag   end if
!
        kmin(i)=kmin(i)+1
        temp(i,j,1,n)=temp(i,j,k,n)
        saln(i,j,1,n)=saln(i,j,k,n)
        th3d(i,j,1,n)=th3d(i,j,k,n)
      end if
      enddo !k
!
      do k=1,kk
      p(i,j,k+1)=max(p(i,j,k),p(i,j,k+1))
      dpo(i,j,k,n)=dp(i,j,k,n)					! diapyc.flx.
      dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,k,n)-dpo(i,j,k,n))	! diapyc.flx.
! --- make sure p is computed from dp, not the other way around (roundoff!)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
      enddo !k
!
      dpmixl(i,j,n)=dp(i,j,1,n)
!diag if (i.eq.itest.and.j.eq.jtest) &
!diag   write (lp,'(i9,2i5,3x,a/(i36,0p,3f10.3,3p,f10.3))') &
!diag   nstep,i+i0,j+j0, &
!diag   'after  diapfl: thickness  salinity temperature density', &
!diag   1,dp(i,j,1,n)*qonem,saln(i,j,1,n),temp(i,j,1,n), &
!diag   th3d(i,j,1,n)+thbase,(k,dp(i,j,k,n)*qonem,saln(i,j,k,n), &
!diag   temp(i,j,k,n),th3d(i,j,k,n)+thbase,k=2,kk)
!
      endif !ip
      enddo !i
!
!diag write (lp,'(i9,7x,1p,e9.2,a)') nstep,salt*1.e-6/g, &
!diag   ' kg salt added in diapfl'
!
      call dpudpv(dpu(1-nbdy,1-nbdy,1,n), &
                  dpv(1-nbdy,1-nbdy,1,n), &
                  p,depthu,depthv, 0,0)
!
!cc      write (lp,'(i9,3x,a)') nstep,'exiting    d i a p f l'
      return
      end
!
!> Revision history:
!>
!> Mar. 2000 - conversion to SI units
!> May  2000 - converted T/S advection equations to flux form
!> Jul. 2000 - added diapf2j for OpenMP parallelization.
!> Aug. 2000 - adapted diapf3 from micom 2.8 to run within hycom 1.0
!> Jan. 2004 - added latdiw to diapf1
!> Mar. 2004 - added thkriv river support to diapf1
!> Mar. 2005 - added [ts]ofset to diapf1 and diapf2
!> Feb. 2009 - modified latdiw to use array diwlat
!> Mar. 2010 - added diwbot
!> Apr. 2010 - botdiw based on Decloedt&Luther in KPP&MY, removed latdiw
!> Oct  2010 - replaced two calls to dsiglocdX with one call at mid-point
!> Aug  2011 - replaced dpold,dpoldm with dpo
!> July 2013 - added diws and diwm
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> July 2017 - added needed halo updates (xctilr)
!> Dec. 2018 - add /* USE_NUOPC_CESMBETA */ macro and riv_input for coupled simulation
