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
      module mod_barotp
      use mod_xc    ! HYCOM communication interface
!
! --- module for barotp and related routines
!
      public  :: barotp, barotp_init

      contains

      subroutine barotp(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
      use mod_tides      ! HYCOM tides
#if defined(STOKES)
      use mod_stokes     !    HYCOM Stokes Drift
#endif
      implicit none
!
      integer m,n
!
! --- ------------------------------------------------------------------------
! --- advance barotropic equations.
! ---   on entry: -n- is time t-dt, -m- is time t
! ---   on exit:                    -m- is time t, -n- is time t+dt
! ---   time level 3 is only used internally (n and m are always 1 or 2).
!
! --- LeapFrog version based on:
! ---   Y. Morel, Baraille, R., Pichon A. (2008) "Time splitting and
! ---   linear stability of the slow part of the barotropic component", 
! ---   Ocean Modeling, 23, pp 73-81.
! --- ------------------------------------------------------------------------
!
      logical, parameter ::  lpipe_barotp=.true.       !usually .false.
      logical, parameter :: ldebug_barotp=.false.      !usually .false.
!
      character text*12
!
      real    q,pbudel,pbvdel,utndcy,vtndcy,wblpf
      real    d11,d12,d21,d22,ubp,vbp,z1
      real    xmin(2)
      real    sminny(jdm,2)
      real*8  sump
      integer i,j,l,lll,ml,nl,mn,lstep1,margin,mbdy,k,icof
!	 & ,iffstep
      logical ldrag
!	  data iffstep/0/
!	  save iffstep
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
              pbavo,ubavo,vbavo,displd,gslpr,gtide, &
              flxloc,flyloc,uflxba,vflxba
!
      if     (.not.allocated(pbavo)) then
        allocate( &
                pbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                ubavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                vbavo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               displd(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                gslpr(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy),  &
                gtide(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 6*(idm+2*nbdy)*(jdm+2*nbdy) )
                pbavo = r_init
                ubavo = r_init
                vbavo = r_init
               displd = r_init
                gslpr = r_init
                gtide = r_init
        if     (btrmas) then
          allocate( &
                 flxloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 flyloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 uflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 vflxba(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
          call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) )
                 flxloc = r_init
                 flyloc = r_init
                 uflxba = r_init
                 vflxba = r_init
        endif
      endif
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
              pbavo,ubavo,vbavo,displd,gslpr,gtide, &
              flxloc,flyloc,uflxba,vflxba
#endif
!
      mbdy = 6
!
      margin = mbdy
! --- atmospheric pressure forcing
!
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if     (mslprf) then
            if     (natm.eq.2) then
              gslpr(i,j) =( mslprs(i,j,l0)*w0+ &
                            mslprs(i,j,l1)*w1 )*svref
            else
              gslpr(i,j) =( mslprs(i,j,l0)*w0+ &
                            mslprs(i,j,l1)*w1+ &
                            mslprs(i,j,l2)*w2+ &
                            mslprs(i,j,l3)*w3 )*svref
            endif !natm
          else
            gslpr(i,j) = 0.0
          endif !mslprf
!
          if (SEA_P) then
!
! ---       tidal body forcing, including Scalar SAL,
! ---       SAL should only be applied to the mass anomally, i.e. to
! ---       non-steric SSH, so use sshflg>0 except in tides-only cases.
            if     (tidflg.gt.0 .and. sshflg.eq.0) then !tides
              gtide(i,j)=-g*etide(i,j) &
                         -salfac(i,j)* srfhgt(i,j)
            elseif (tidflg.gt.0 .and. sshflg.eq.1) then !tides
              gtide(i,j)=-g*etide(i,j) &
                         -salfac(i,j)*(srfhgt(i,j)-steric(i,j))
            elseif (tidflg.gt.0 .and. sshflg.eq.2) then !tides
              gtide(i,j)=-g*etide(i,j) &
                         -salfac(i,j)*(srfhgt(i,j)-montg1(i,j))
            else
              gtide(i,j)=0.0
            endif !tides (sshflg)
          endif !ip
        enddo !i
      enddo !j
!
! --- utotn,vtotn from momtum is time step t-1 to t+1 barotropic tendency
      call xctilr(utotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_uv)
      call xctilr(vtotn(  1-nbdy,1-nbdy    ),1, 1, 6,6, halo_vv)
!
      if     (lpipe .and. lpipe_barotp) then
! ---   compare two model runs.
        call pipe_compare_sym2(utotn, iu,'barotp:utotn', &
                               vtotn, iv,'barotp:vtotn')
        call pipe_compare_sym1(pvtrop,iq,'barotp:pvtrp')
      endif
!
! --- explicit time integration of barotropic flow (forward-backward scheme)
! --- in order to combine forward-backward scheme with leapfrog treatment of
! --- coriolis term, v-eqn must be solved before u-eqn every other time step
!
      if (btrmas) then
        uflxba(:,:)   = 0.0
        vflxba(:,:)   = 0.0
        flxloc(:,:)   = 0.0
        flyloc(:,:)   = 0.0
      endif
!
      if     (btrlfr) then  !always true for btrmas
!
        if     (delt1.ne.baclin) then  !not on very 1st time step
! ---     start at time level t-dt and go to t+dt.
          lstep1 = lstep + lstep  !more stable, but also more expensive
          icof   = 2
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = pbavg(i,j,n)  !save t-1 for RA filter
              ubavo(i,j)   = ubavg(i,j,n)  !save t-1 for RA filter
              vbavo(i,j)   = vbavg(i,j,n)  !save t-1 for RA filter
!
              pbavg(i,j,3) = pbavg(i,j,n)
              ubavg(i,j,3) = ubavg(i,j,n)
              vbavg(i,j,3) = vbavg(i,j,n)
            enddo !i
          enddo !j
        else !1st time step
! ---     start at time level t and go to t+dt.
          lstep1 = lstep
          icof   = 1
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1,jj
            do i=1,ii
              pbavo(i,j)   = 0.0 !makes correct mean height safe
              pbavg(i,j,n) = pbavg(i,j,m)
              ubavg(i,j,n) = ubavg(i,j,m)
              vbavg(i,j,n) = vbavg(i,j,m)
              pbavg(i,j,3) = pbavg(i,j,m)
              ubavg(i,j,3) = ubavg(i,j,m)
              vbavg(i,j,3) = vbavg(i,j,m)
            enddo !i
          enddo !j
        endif !usual:1st time step
      else
! ---   start at time level t    and go to t+dt.
        lstep1 = lstep          !original, less stable, method
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            pbavo(i,j)   = 0.0 !makes correct mean height safe
            pbavg(i,j,n) = pbavg(i,j,m)
            ubavg(i,j,n) = ubavg(i,j,m)
            vbavg(i,j,n) = vbavg(i,j,m)
          enddo !i
        enddo !j
      endif !btrlfr
!
      ldrag = tidflg.gt.0 .and. drgscl.ne.0.0 .and. thkdrg.eq.0.0
!
      if     (ldrag) then
        displd(:,:) = 0.0
      endif
!
! --- time step loop
!
      if     (btrlfr) then
        wblpf = 0.0   !1st minor time step, lll=1, only
      else
        wblpf = wbaro
      endif
!
      do 840 lll=1,lstep1,2
!
      call xctilr(pbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_ps)
      call xctilr(ubavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,1  ),1, 3, 6,6, halo_vv)
!
! --- odd minor time step.
!
      ml=n
      nl=3
!
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1( &
          pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2( &
          ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn', &
          vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1( &
          pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2( &
          ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm', &
          vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
!
! --- continuity equation, and tidal drag on p-grid
!
! --- rhs: pbavg, ubavg+, vbavg+
! --- lhs: pbavg
!
      if     (btrmas) then
!
        margin = mbdy
!
!$OMP   PARALLEL DO PRIVATE(j,l,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
              uflxba(i,j) = uflxba(i,j) + &
                  coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
            endif !iu
            if (SEA_V) then
              flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
              vflxba(i,j) = vflxba(i,j) + &
                  coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
            endif !iv
          enddo !i
          do l=1,isu(j)  !ok
            i=ifu(j,l)-1
            if (i.ge.1-margin) then
                if (iuopn(i,j).ne.0) then
                    flxloc(i,j) = ubavg(i,j,ml)*scuy(i,j)* &
                       max(0.,pbot(i,j))               
                    uflxba(i,j) = uflxba(i,j) + &
                        coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
                endif !iuopn
            endif !i
            i=ilu(j,l)+1
            if (i.le.ii+margin) then
                if (iuopn(i,j).ne.0) then
                    flxloc(i,j) = ubavg(i,j,ml)*scuy(i,j)* &
                       max(0.,pbot(i-1,j))            
                    uflxba(i,j) = uflxba(i,j) + &
                        coeflx(lll+1,icof)*(1.0+wblpf)*flxloc(i,j)
                endif !iuopn
            endif !i
          enddo !l
        enddo !j
!
        do i=1-margin,ii+margin
          do l=1,jsv(i) !ok
            j=jfv(i,l)-1
            if (j.ge.1-margin) then
                if (ivopn(i,j).ne.0) then
                    flyloc(i,j) = vbavg(i,j,ml)*scvx(i,j)* &
                       max(0.,pbot(i,j))               
                    vflxba(i,j) = vflxba(i,j) + &
                        coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
                endif !ivopn
            endif !j
            j=jlv(i,l)+1
            if (j.le.jj+margin) then
                if (ivopn(i,j).ne.0) then
                    flyloc(i,j) = vbavg(i,j,ml)*scvx(i,j)* &
                       max(0.,pbot(i,j-1))               
                    vflxba(i,j) = vflxba(i,j) + &
                        coeflx(lll+1,icof)*(1.0+wblpf)*flyloc(i,j)
                endif !ivopn
            endif !j
          enddo !l
        enddo  !i
!
        margin = mbdy - 1
!

!$OMP PARALLEL DO PRIVATE(j,l,i,ubp,vbp,d11,d12,d21,d22,q) &
!$OMP            SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do i=1-margin,ii+margin
             if (SEA_P) then
                pbavg(i,j,nl)= &
                  ((1.0-wblpf)*pbavg(i,j,ml)+ &
                        wblpf *pbavg(i,j,nl) )- &
                   (1.0+wblpf)*dlt*(flxloc(i+1,j)-flxloc(i,j) + &
                                    flyloc(i,j+1)-flyloc(i,j)  )* &
                                   scp2i(i,j)
!
                if     (ldrag) then
!
! ---             tidal drag tensor on p-grid:
! ---               ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
! ---               vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
! ---             solve implicitly by inverting the matrix:
! ---                1+(dlt/H)*t.11    (dlt/H)*t.12
! ---                  (dlt/H)*t.21  1+(dlt/H)*t.22
! ---             use depths (H) rather than onem*pbavg (h) for stability.
!
                  ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
                  vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
                  d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
                  d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
                  d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
                  d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
                  q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
! ---             set util5,util6 to the ubavg,vbavg drag increment
                  util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
                  util6(i,j) = q*(ubp*d21+vbp*(1.0-d12)) - vbp
! ---             add an explicit antidrag correction
!                 util5(i,j) = util5(i,j) - (d11*untide(i,j)+
!    &                                       d12*vntide(i,j) )
!                 util6(i,j) = util6(i,j) - (d21*untide(i,j)+
!    &                                       d22*vntide(i,j) )
!
!                 if (ldebug_barotp .and.
!    &                i.eq.itest.and.j.eq.jtest) then
!                   write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
!    &                nstep,i+i0,j+j0,lll,
!    &                'ubp,new,vbp,new =',
!    &              ubp,ubp+util5(i,j),
!    &              vbp,vbp+util6(i,j)
!                 endif !debug
                else
                  util5(i,j) = 0.0
                  util6(i,j) = 0.0
                endif !ldrag
            endif 
          enddo !i
        enddo !j
!
      else !.not.btrmas
!
      margin = mbdy - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel, &
!$OMP                     ubp,vbp,d11,d12,d21,d22,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
!
!   Barotropic Stokes flow included here
!
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))* &
                          (depthu(i+1,j)*scuy(i+1,j)) &
                    -(ubavg(i,  j,ml)+usdbavg(i,  j))* &
                          (depthu(i,  j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))* &
                          (depthv(i,j+1)*scvx(i,j+1)) &
                    -(vbavg(i,j,  ml)+vsdbavg(i,j  ))* &
                          (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j)) &
                     -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1)) &
                     -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))
#endif
            pbavg(i,j,nl)= &
              ((1.-wblpf)*pbavg(i,j,ml)+ &
                   wblpf *pbavg(i,j,nl) )- &
               (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
!
            if     (ldrag) then
!
! ---         tidal drag tensor on p-grid:
! ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
! ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
! ---         solve implicitly by inverting the matrix:
! ---            1+(dlt/H)*t.11    (dlt/H)*t.12
! ---              (dlt/H)*t.21  1+(dlt/H)*t.22
! ---         use depths (H) rather than onem*pbavg (h) for stability.
!
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
! ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
! ---         add an explicit antidrag correction
!             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
!    &                                   d12*vntide(i,j) )
!             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
!    &                                   d22*vntide(i,j) )
! ---         dissipation per m^2
              displd(i,j) = displd(i,j) + &
                            (ubp*util5(i,j) + vbp*util6(i,j))* &
                            depths(i,j)*rhoref/dlt
!
!             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
!    &            nstep,i+i0,j+j0,lll,
!    &            'ubp,new,vbp,new =',
!    &          ubp,ubp+util5(i,j),
!    &          vbp,vbp+util6(i,j)
!             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
      endif !btrmas:else
!
      mn=ml
!
! --- u momentum equation, 1st
!
! --- rhs: pbavg+, vbavg+, pvtrop+
! --- lhs: ubavg
!
      margin = margin - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,utndcy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-svref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)- &
                          (gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)- &
                          (gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)+ &
             ((vbavg(i  ,j,  mn)*depthv(i  ,j) &
              +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+ &
              (vbavg(i-1,j,  mn)*depthv(i-1,j) &
              +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))* &
             (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
!
            ubavg(i,j,nl)= &
              ((1.-wblpf)*ubavg(i,j,ml)+ &
                   wblpf *ubavg(i,j,nl) )+ &
               (1.+wblpf)*dlt*(utndcy+utotn(i,j))+ &
                      0.5*(util5(i,j)+util5(i-1,j))
!
!           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
!    &          nstep,i+i0,j+j0,lll,
!    &          'u_old,u_new,p_grad,t_g,m_g,corio,u_star,drag =',
!    &          ubavg(i,j,ml),ubavg(i,j,nl),
!    &           -svref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
!    &                 -(gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)*dlt,
!    &                 -(gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)*dlt,
!    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
!    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
!    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
!    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
!    &          *(pvtrop(i,j)+pvtrop(i,j+1))
!    &          *.125 * dlt,utotn(i,j) * dlt,
!    &          0.5*(util5(i,j)+util5(i-1,j))
!           endif !debug
          endif !iu
        enddo !i
      enddo !j
!
      mn = nl
!
! --- v momentum equation, 2nd
! --- rhs: pbavg+, ubavg+, pvtrop+
! --- lhs: vbavg
!
      margin = margin - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-svref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)- &
                          (gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)- &
                          (gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)- &
             ((ubavg(i,  j  ,mn)*depthu(i,  j  ) &
              +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+ &
              (ubavg(i,  j-1,mn)*depthu(i,  j-1) &
              +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))* &
             (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
!
            vbavg(i,j,nl)= &
              ((1.-wblpf)*vbavg(i,j,ml)+ &
                   wblpf *vbavg(i,j,nl) )+ &
               (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+ &
                      0.5*(util6(i,j)+util6(i,j-1))
!
!           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
!    &          nstep,i+i0,j+j0,lll,
!    &          'v_old,v_new,p_grad,t_g,m_g,corio,v_star,drag =',
!    &          vbavg(i,j,ml),vbavg(i,j,nl),
!    &          -svref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
!    &                -(gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)*dlt,
!    &                -(gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)*dlt,
!    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
!    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
!    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
!    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
!    &          *(pvtrop(i,j)+pvtrop(i+1,j))
!    &          *.125 * dlt, vtotn(i,j) * dlt,
!    &          0.5*(util6(i,j)+util6(i,j-1))
!           endif !debug
          endif !iv
        enddo !i
      enddo !j
!
!     if     (ldebug_barotp) then
!       call xcsync(flush_lp)
!     endif
!
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll)
      endif
!
! --- even minor time step.
!
      ml=3
      nl=n
      wblpf = wbaro  !used for all subsequent time steps: lll=2,lstep1
!
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare_sym1( &
          pbavg(1-nbdy,1-nbdy,nl),ip,'barot+:pbavn')
        call pipe_compare_sym2( &
          ubavg(1-nbdy,1-nbdy,nl),iu,'barot+:ubavn', &
          vbavg(1-nbdy,1-nbdy,nl),iv,'barot+:vbavn')
        call pipe_compare_sym1( &
          pbavg(1-nbdy,1-nbdy,ml),ip,'barot+:pbavm')
        call pipe_compare_sym2( &
          ubavg(1-nbdy,1-nbdy,ml),iu,'barot+:ubavm', &
          vbavg(1-nbdy,1-nbdy,ml),iv,'barot+:vbavm')
      endif
!
! --- continuity equation
!
! --- rhs: pbavg, ubavg+, vbavg+
! --- lhs: pbavg
!
      if     (btrmas) then
!
         margin = mbdy
!
!$OMP    PARALLEL DO PRIVATE(j,l,i) &
!$OMP             SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do i=1-margin,ii+margin
             if (SEA_U) then
                 flxloc(i,j) = ubavg(i,j,ml)*(depthu(i,j)*scuy(i,j))
                 uflxba(i,j) = uflxba(i,j) + &
                     coeflx(lll+2,icof)*(1.0+wblpf)*flxloc(i,j)
             endif 
           enddo !i
         enddo !j
!
!$OMP    PARALLEL DO PRIVATE(j,l,i) &
!$OMP             SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do i=1-margin,ii+margin
             if (SEA_V) then
                flyloc(i,j) = vbavg(i,j,ml)*(depthv(i,j)*scvx(i,j))
                vflxba(i,j) = vflxba(i,j) + &
                              coeflx(lll+2,icof)*(1.0+wblpf)*flyloc(i,j)
            endif 
          enddo !i
        enddo !j
!
         margin = mbdy - 1
!
!$OMP    PARALLEL DO PRIVATE(j,l,i,ubp,vbp,d11,d12,d21,d22,q) &
!$OMP             SCHEDULE(STATIC,jblk)
         do j=1-margin,jj+margin
           do i=1-margin,ii+margin
             if (SEA_P) then
                pbavg(i,j,nl)= &
                  ((1.0-wblpf)*pbavg(i,j,ml)+ &
                        wblpf *pbavg(i,j,nl) )- &
                   (1.0+wblpf)*dlt*(flxloc(i+1,j)-flxloc(i,j) + &
                                    flyloc(i,j+1)-flyloc(i,j)  )* &
                                   scp2i(i,j)
!
                if     (ldrag) then
!
! ---             tidal drag tensor on p-grid:
! ---               ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
! ---               vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
! ---             solve implicitly by inverting the matrix:
! ---                1+(dlt/H)*t.11    (dlt/H)*t.12
! ---                  (dlt/H)*t.21  1+(dlt/H)*t.22
! ---             use depths (H) rather than onem*pbavg (h) for stability.
!
                  ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
                  vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
                  d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
                  d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
                  d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
                  d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
                  q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
! ---             set util5,util6 to the ubavg,vbavg drag increment
                  util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
                  util6(i,j) = q*(ubp*d21+vbp*(1.0-d12)) - vbp
! ---             add an explicit antidrag correction
!                 util5(i,j) = util5(i,j) - (d11*untide(i,j)+
!    &                                       d12*vntide(i,j) )
!                 util6(i,j) = util6(i,j) - (d21*untide(i,j)+
!    &                                       d22*vntide(i,j) )
!
!                 if (ldebug_barotp .and.
!    &                i.eq.itest.and.j.eq.jtest) then
!                   write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
!    &                nstep,i+i0,j+j0,lll,
!    &                'ubp,new,vbp,new =',
!    &              ubp,ubp+util5(i,j),
!    &              vbp,vbp+util6(i,j)
!                 endif !debug
                else
                  util5(i,j) = 0.0
                  util6(i,j) = 0.0
                endif !ldrag
            endif 
          enddo !i
        enddo !j
!
        if     (lpipe .and. lpipe_barotp) then
! ---     compare two model runs.
          text = 'uflxba.ml      '
          call pipe_compare_sym1(uflxba,iu,text)
          text = 'vflxba.ml      '
          call pipe_compare_sym1(vflxba,iv,text)
          text = 'flxloc.ml   '
          call pipe_compare_sym1(flxloc,iu,text)
          text = 'flyloc.ml   '
          call pipe_compare_sym1(flyloc,iv,text)
        endif  !lpipe
!
      else !.not.btrmas
!
      margin = mbdy - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,pbudel,pbvdel, &
!$OMP                     ubp,vbp,d11,d12,d21,d22,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
#if defined(STOKES)
!
!           Barotropic Stokes flow included here
!
            pbudel = (ubavg(i+1,j,ml)+usdbavg(i+1,j))* &
                          (depthu(i+1,j)*scuy(i+1,j)) &
                    -(ubavg(i,  j,ml)+usdbavg(i,  j))* &
                           (depthu(i ,j)*scuy(i,  j))
            pbvdel = (vbavg(i,j+1,ml)+vsdbavg(i,j+1))* &
                          (depthv(i,j+1)*scvx(i,j+1)) &
                    -(vbavg(i,j,  ml)+vsdbavg(i,j  ))* &
                          (depthv(i,j  )*scvx(i,j  ))
#else
            pbudel =  ubavg(i+1,j,ml)*(depthu(i+1,j)*scuy(i+1,j)) &
                     -ubavg(i  ,j,ml)*(depthu(i  ,j)*scuy(i  ,j))
            pbvdel =  vbavg(i,j+1,ml)*(depthv(i,j+1)*scvx(i,j+1)) &
                     -vbavg(i,j  ,ml)*(depthv(i,j  )*scvx(i,j  ))

#endif
            pbavg(i,j,nl)= &
              ((1.-wblpf)*pbavg(i,j,ml)+ &
                   wblpf *pbavg(i,j,nl) )- &
               (1.+wblpf)*dlt*(pbudel + pbvdel)*scp2i(i,j)
!
            if     (ldrag) then
! ---         tidal drag tensor on p-grid:
! ---           ub = ub - (dlt/H)*(t.11*ub + t.12*vb)
! ---           vb = vb - (dlt/H)*(t.21*ub + t.22*vb)
! ---         solve implicitly by inverting the matrix:
! ---            1+(dlt/H)*t.11    (dlt/H)*t.12
! ---              (dlt/H)*t.21  1+(dlt/H)*t.22
! ---         use depths (H) rather than onem*pbavg (h) for stability.
!
              ubp = 0.5*(ubavg(i+1,j,nl)+ubavg(i,j,nl))
              vbp = 0.5*(vbavg(i,j+1,nl)+vbavg(i,j,nl))
              d11 = -dlt/depths(i,j) * drgten(1,1,i,j)
              d12 = -dlt/depths(i,j) * drgten(1,2,i,j)
              d21 = -dlt/depths(i,j) * drgten(2,1,i,j)
              d22 = -dlt/depths(i,j) * drgten(2,2,i,j)
              q   = 1.0/((1.0-d11)*(1.0-d22)-d12*d21)
! ---         set util5,util6 to the ubavg,vbavg drag increment
              util5(i,j) = q*(ubp*(1.0-d22)+vbp*d12) - ubp
              util6(i,j) = q*(ubp*d21+vbp*(1.0-d11)) - vbp
! ---         add an explicit antidrag correction
!             util5(i,j) = util5(i,j) - (d11*untide(i,j)+
!    &                                   d12*vntide(i,j) )
!             util6(i,j) = util6(i,j) - (d21*untide(i,j)+
!    &                                   d22*vntide(i,j) )
! ---         dissipation per m^2
              displd(i,j) = displd(i,j) + &
                            (ubp*util5(i,j) + vbp*util6(i,j))* &
                            depths(i,j)*rhoref/dlt
!
!             if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!               write (lp,'(i9,2i5,i3,3x,a,4g15.6)')
!    &            nstep,i+i0,j+j0,lll+1,
!    &            'ubp,new,vbp,new =',
!    &          ubp,ubp+util5(i,j),
!    &          vbp,vbp+util6(i,j)
!             endif !debug
            else
              util5(i,j) = 0.0
              util6(i,j) = 0.0
            endif !ldrag
          endif !ip
        enddo !i
      enddo !j
      endif !btrmas:else
!
      mn=ml
!
! --- v momentum equation, 1st
!
! --- rhs: pbavg+, ubavg+, pvtrop+
! --- lhs: vbavg
!
      margin = margin - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,vtndcy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_V) then
            vtndcy=-svref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)- &
                          (gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)- &
                          (gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)- &
             ((ubavg(i,  j  ,mn)*depthu(i,  j  ) &
              +ubavg(i+1,j  ,mn)*depthu(i+1,j  ))+ &
              (ubavg(i,  j-1,mn)*depthu(i,  j-1) &
              +ubavg(i+1,j-1,mn)*depthu(i+1,j-1)))* &
             (0.125*(pvtrop(i,j)+pvtrop(i+1,j)))
!
            vbavg(i,j,nl)= &
              ((1.-wblpf)*vbavg(i,j,ml)+ &
                   wblpf *vbavg(i,j,nl))+ &
               (1.+wblpf)*dlt*(vtndcy+vtotn(i,j))+ &
                      0.5*(util6(i,j)+util6(i,j-1))
!
!           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!             write (lp,'(i9,2i5,i3,3x,a,8g15.6)')
!    &          nstep,i+i0,j+j0,lll+1,
!    &          'v_old,v_new,p_grad,t_g,m_g,corio,v_star,drag =',
!    &          vbavg(i,j,ml),vbavg(i,j,nl),
!    &          -svref*(pbavg(i,j,nl)-pbavg(i,j-1,nl))*scvyi(i,j)*dlt,
!    &                -(gtide(i,j)   -gtide(i,j-1)   )*scvyi(i,j)*dlt,
!    &                -(gslpr(i,j)   -gslpr(i,j-1)   )*scvyi(i,j)*dlt,
!    &          -(ubavg(i,  j  ,mn)*depthu(i,j  )
!    &           +ubavg(i+1,j  ,mn)*depthu(i+1,j  )
!    &           +ubavg(i,  j-1,mn)*depthu(i,j-1)
!    &           +ubavg(i+1,j-1,mn)*depthu(i+1,j-1))
!    &          *(pvtrop(i,j)+pvtrop(i+1,j))
!    &          *.125 * dlt, vtotn(i,j) * dlt,
!    &          0.5*(util6(i,j)+util6(i,j-1))
!           endif !debug
          endif !iv
        enddo !i
      enddo !j
!
      mn=nl
!
! --- u momentum equation, 2nd
!
! --- rhs: pbavg+, vbavg+, pvtrop+
! --- lhs: ubavg
!
      margin = margin - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,utndcy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            utndcy=-svref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)- &
                          (gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)- &
                          (gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)+ &
             ((vbavg(i  ,j,  mn)*depthv(i  ,j) &
              +vbavg(i  ,j+1,mn)*depthv(i  ,j+1))+ &
              (vbavg(i-1,j,  mn)*depthv(i-1,j) &
              +vbavg(i-1,j+1,mn)*depthv(i-1,j+1)))* &
             (0.125*(pvtrop(i,j)+pvtrop(i,j+1)))
!
            ubavg(i,j,nl)= &
              ((1.-wblpf)*ubavg(i,j,ml)+ &
                   wblpf *ubavg(i,j,nl) )+ &
               (1.+wblpf)*dlt*(utndcy+utotn(i,j))+ &
                      0.5*(util5(i,j)+util5(i-1,j))
!
!           if (ldebug_barotp .and. i.eq.itest.and.j.eq.jtest) then
!             write (lp,'(i9,2i5,i3,3x,a,7g15.6)')
!    &          nstep,i+i0,j+j0,lll+1,
!    &          'u_old,u_new,p_grad,t_g,m_g,corio,u_star,drag =',
!    &          ubavg(i,j,ml),ubavg(i,j,nl),
!    &          -svref*(pbavg(i,j,nl)-pbavg(i-1,j,nl))*scuxi(i,j)*dlt,
!    &                -(gtide(i,j)   -gtide(i-1,j)   )*scuxi(i,j)*dlt,
!    &                -(gslpr(i,j)   -gslpr(i-1,j)   )*scuxi(i,j)*dlt,
!    &          (vbavg(i  ,j,  mn)*depthv(i  ,j)
!    &          +vbavg(i  ,j+1,mn)*depthv(i  ,j+1)
!    &          +vbavg(i-1,j,  mn)*depthv(i-1,j)
!    &          +vbavg(i-1,j+1,mn)*depthv(i-1,j+1))
!    &          *(pvtrop(i,j)+pvtrop(i,j+1))
!    &          *.125 * dlt,utotn(i,j) * dlt,
!    &          0.5*(util5(i,j)+util5(i-1,j))
!           endif !debug
          endif !iu
        enddo !i
      enddo !j
!
!     if     (ldebug_barotp) then
!       call xcsync(flush_lp)
!     endif
!
#if ! defined(RELO)
      if     (lbflag.eq.1) then
        call latbdp( nl)
      elseif (lbflag.eq.3) then
        call latbdf( nl,lll+1)
      endif
#endif
      if     (lbflag.eq.2) then
        call latbdt( nl,lll+1)
      elseif (lbflag.eq.4) then
        call latbdtf(nl,lll+1)
      endif
!
 840  continue  ! lll=1,lstep1,2
!
      if     (ldrag) then  !disp_count updated in momtum
        displd_mn(:,:) = displd_mn(:,:) + displd(:,:)/real(lstep1)
      endif
!
      if     (lbflag.eq.1) then
!
! ---   correct mean height.
! ---   this should not be required - so there may be a bug in the bc.
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1-margin,ii+margin
            if (SEA_P) then
              util1(i,j) = pbavg(i,j,n)*scp2(i,j)
            endif !ip
          enddo !i
        enddo !j
        call xcsum(sump, util1,ipa)
        q = sump/area
!
! ---   rhs: pbavg
! ---   lhs: pbavg
!
        margin = 0
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pbavo(i,j)   = pbavo(i,j)   - q
              pbavg(i,j,1) = pbavg(i,j,1) - q
              pbavg(i,j,2) = pbavg(i,j,2) - q
              pbavg(i,j,3) = pbavg(i,j,3) - q
            endif !ip
          enddo !i
        enddo !j
      endif
      if     (lpipe .and. lpipe_barotp) then
        call pipe_compare(pbavg(1-nbdy,1-nbdy,1), ip,'barotp:pbav1')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,2), ip,'barotp:pbav2')
        call pipe_compare(pbavg(1-nbdy,1-nbdy,3), ip,'barotp:pbav3')
      endif
!
      if     (btrlfr .and. delt1.ne.baclin) then  !not on very 1st time step
! ---   Robert-Asselin time filter 
!$OMP   PARALLEL DO PRIVATE(j,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            q = 0.5*ra2fac*(    pbavo(i,j)   +    & !t-1
                                pbavg(i,j,n) -    & !t+1
                            2.0*pbavg(i,j,m)  )  !t
            pbavg(i,j,m) = pbavg(i,j,m) + q
            q = 0.5*ra2fac*(    ubavo(i,j)   +    & !t-1
                                ubavg(i,j,n) -    & !t+1
                            2.0*ubavg(i,j,m)  )  !t
            ubavg(i,j,m) = ubavg(i,j,m) + q
            q = 0.5*ra2fac*(    vbavo(i,j)   +    & !t-1
                                vbavg(i,j,n) -    & !t+1
                            2.0*vbavg(i,j,m)  )  !t
            vbavg(i,j,m) = vbavg(i,j,m) + q
          enddo !i
        enddo !j
      endif !btrlfr & not on very 1st time step
!
! --- always update 1 + eta
!
!$OMP PARALLEL DO PRIVATE(j,l,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            oneta(i,j,n)  = max( oneta0, 1.0 + pbavg(i,j,n)/pbot(i,j) )  !t+1
            oneta(i,j,m)  = max( oneta0, 1.0 + pbavg(i,j,m)/pbot(i,j) )  !t&RA
          endif 
        enddo !i
      enddo !j
!
      if (btrmas) then
!
        uflxba(:,:) = dlt*uflxba(:,:)/delt1
        vflxba(:,:) = dlt*vflxba(:,:)/delt1
!
        call xctilr(uflxba(  1-nbdy,1-nbdy),  1, 1, 2,2, halo_uv)
        call xctilr(vflxba(  1-nbdy,1-nbdy),  1, 1, 2,2, halo_vv)
        call xctilr(onetacnt(1-nbdy,1-nbdy),  1, 1, 2,2, halo_ps)
!
!                Be sure that the sum over the vertical of the baroclinic
!                mass fluxes is exactly the same as the barotropic fluxes
!                --------------------------------------------------------
!
        if     (lpipe .and. lpipe_barotp) then
! ---     compare two model runs.
          text = 'onetacnt    '
          call pipe_compare_sym1(onetacnt,ip,text)
          text = 'uflxba      '
          call pipe_compare_sym1(uflxba,iu,text)
          text = 'vflxba      '
          call pipe_compare_sym1(vflxba,iv,text)
        endif  !lpipe
!
        margin = 1
        util1(:,:) = 0.0
        util2(:,:) = 0.0
        uflux(:,:) = 0.0
        vflux(:,:) = 0.0
        do k =1,kk
!$OMP     PARALLEL DO PRIVATE(j,l,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_P) then
                dp(i,j,k,n)=onetacnt(i,j)*dp(i,j,k,n)  !onetacnt from cnuity
                p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              endif  !ip
              if (SEA_U) then
                util1(i,j) = util1(i,j) + uflx(i,j,k)
              endif  !iu
              if (SEA_V) then
                util2(i,j) = util2(i,j) + vflx(i,j,k)
              endif  !iv
            enddo !i
!
            do l=1,isu(j) !ok
              i=ifu(j,l)-1
              if (i.ge.1-margin) then
                  if (iuopn(i,j).ne.0) then
                      util1(i,j) = util1(i,j)+uflx(i,j,k)
                  endif
              endif
              i=ilu(j,l)+1
              if (i.le.ii+margin) then
                  if (iuopn(i,j).ne.0) then
                      util1(i,j) = util1(i,j)+uflx(i,j,k)
                  endif
              endif
            enddo !l
          enddo !j
!$OMP     END PARALLEL DO
!
          do i=1-margin,ii+margin
            do l=1,jsv(i) !ok
              j=jfv(i,l)-1
              if (j.ge.1-margin) then
                if (ivopn(i,j).ne.0) then
                    util2(i,j) = util2(i,j)+vflx(i,j,k)
                endif
              endif
              j=jlv(i,l)+1
              if (j.le.jj+margin) then
                if (ivopn(i,j).ne.0) then
                    util2(i,j) = util2(i,j)+vflx(i,j,k)
                endif
              endif
            enddo !l
          enddo !i
!
          if     (lpipe .and. lpipe_barotp) then
! ---       compare two model runs.
            write (text,'(a9,i3)') 'dpcnt kn=',k
            call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
            write (text,'(a9,i3)') 'util1  k=',k
            call pipe_compare_sym1(util1,iu,text)
            write (text,'(a9,i3)') 'util2  k=',k
            call pipe_compare_sym1(util2,iv,text)
          endif  !lpipe
!
        enddo  !k
!
        do k= 1,kk
!
!$OMP   PARALLEL DO PRIVATE(j,l,i,z1) &
!$OMP            SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
!
            do i=1-margin,ii+margin
              if (SEA_U) then
                if (uflxba(i,j)-util1(i,j).ge.0.) then
                    z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)
                else
                    z1=dp(i  ,j,k,n)/p(i  ,j,kk+1)
                endif
                util3(i,j)=z1  !lpipe
                uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
              endif  !iu
              if (SEA_V) then
                if (vflxba(i,j)-util2(i,j).ge.0.) then
                    z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)
                else
                    z1=dp(i,j  ,k,n)/p(i,j  ,kk+1)
                endif
                util4(i,j)=z1  !lpipe
                vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
              endif  !iv
            enddo !i
!
            do l=1,isu(j) !ok
              i=ifu(j,l)-1
              if (i.ge.1-margin  ) then
                if (iuopn(i,j).ne.0) then
                  z1=dp(i,j,k,n)/p(i,j,kk+1)
                  util3(i,j)=z1  !lpipe
                  uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                  uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
                endif
              endif
              i=ilu(j,l)+1
              if (i.le.ii+margin ) then
                if (iuopn(i,j).ne.0) then
                  z1=dp(i-1,j,k,n)/p(i-1,j,kk+1)
                  util3(i,j)=z1  !lpipe
                  uflux(i,j)=(uflxba(i,j)-util1(i,j))*z1
                  uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
                endif
              endif
            enddo  !l
!
          enddo !j
!$OMP     END PARALLEL DO
!
          do i=1-margin,ii+margin
            do l=1,jsv(i) !ok
              j=jfv(i,l)-1
              if (j.ge.1-margin  ) then
                if (ivopn(i,j).ne.0) then
                  z1=dp(i,j,k,n)/p(i,j,kk+1)
                  util4(i,j)=z1  !lpipe
                  vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                  vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
                endif
              endif
              j=jlv(i,l)+1
              if (j.le.jj+margin  ) then
                if (ivopn(i,j).ne.0) then
                  z1=dp(i,j-1,k,n)/p(i,j-1,kk+1)
                  util4(i,j)=z1  !lpipe
                  vflux(i,j)=(vflxba(i,j)-util2(i,j))*z1
                  vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
                endif
              endif
            enddo !l
          enddo !i
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_P) then
                  dp(i,j,k,n)=dp(i,j,k,n)- &
                         ((uflux(i+1,j)-uflux(i,j))+ &
                          (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
                  p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              endif 
            enddo !i
          enddo !j
!
          if     (lpipe .and. lpipe_barotp) then
! ---       compare two model runs.
            write (text,'(a9,i3)') 'dp.fl kn=',k
            call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
            write (text,'(a9,i3)') 'z1u    k=',k
            call pipe_compare_sym1(util3,iu,text)
            write (text,'(a9,i3)') 'uflux  k=',k
            call pipe_compare_sym1(uflux,iu,text)
            write (text,'(a9,i3)') 'z1v    k=',k
            call pipe_compare_sym1(util4,iv,text)
            write (text,'(a9,i3)') 'vflux  k=',k
            call pipe_compare_sym1(vflux,iv,text)
          endif  !lpipe
!
        enddo !k
!
!$OMP   PARALLEL DO PRIVATE(j,l,k,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
            do k=1,kk
              do i=1-margin,ii+margin
                if (SEA_P) then
                    if (p(i,j,kk+1) .gt. 0.0) then
                        dp(i,j,k,n)=dp(i,j,k,n)*pbot(i,j)/p(i,j,kk+1)
                    else
                        dp(i,j,k,n)=0.
                        dp(i,j,1,n)=pbot(i,j)
                    endif
                    p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                    if (isopyc .and. k.eq.1) then
                        dpmixl(i,j,n)=dp(i,j,k,n)
                    endif
                endif 
              enddo !i
            enddo !k
          enddo !j
!$OMP   END PARALLEL DO
!
        call xctilr(oneta(1-nbdy,1-nbdy,1),  1, 2, 6,6, halo_ps)
        call xctilr(   dp(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
!
      endif  !btrmas
!
! --- check for clipped oneta
!
      if     (mod(nstep,3).eq.0 .or. diagno) then
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk) !NOCSD
        do j=1,jj
          sminny(j,1)= 999.  !simplifies OpenMP parallelization
          sminny(j,2)= 999.  !simplifies OpenMP parallelization
!DIR$     PREFERVECTOR
          do i=1,ii
            if (SEA_P) then
              sminny(j,1)=min(sminny(j,1),oneta(i,j,1))
              sminny(j,2)=min(sminny(j,2),oneta(i,j,2))
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO !NOCSD
        xmin(1) = minval(sminny(1:jj,1))
        xmin(2) = minval(sminny(1:jj,2))
        call xcminr(xmin(1:2))
!
        do mn= 1,2
          if     (xmin(mn).eq.oneta0) then
           do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if (oneta(i,j,mn).le.oneta0) then
                    write (lp,'(i9,a,2i5,i3,a,f9.6)')  &
                      nstep,' i,j,mn =',i+i0,j+j0,mn, &
                      ' clipped oneta after barotp call ', &
                      oneta(i,j,mn)
                  endif !oneta0
                endif !ip
              enddo !i
            enddo !j
            call xcsync(flush_lp)
          endif !oneta0
        enddo !mn
!
        if (diagno) then
          if     (mnproc.eq.1) then
            write (lp,'(i9,a,f9.6)') &
              nstep, &
              ' min of oneta after barotp:',min(xmin(1),xmin(2))
            call flush(lp)
          endif
        endif !diagno

      endif !every 3 time steps or diagno
!
      return
      end subroutine barotp

      subroutine barotp_init
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
! --- ------------------------------------------------------------------------
! --- calculate coeflx if needed
! --- ------------------------------------------------------------------------
!
      integer i,j,l,ind
!
      real,    allocatable :: xm(:,:),xs(:),worksg(:)
      integer, allocatable :: iworksg(:)
!
      if (btrmas) then
        allocate(      xm(2*lstep+1,2*lstep+1), &
                       xs(2*lstep+1), &
                   worksg(2*lstep+1), &
                  iworksg(2*lstep+1) )
!
        do l= 1,2
          xm(:,:) =  0.0
          xm(1,1) =  1.0
          xm(2,2) =  1.0
          xm(2,1) = -1.0
          do i= 3,l*lstep+1
            xm(i,i)   = 1.0
            xm(i,i-1) = wbaro-1.0
            xm(i,i-2) = -wbaro
          enddo !i
          do j= 1,l*lstep+1
            xs(:) = 0.0
            xs(j) = 1.0
            call s8gefs(xm,2*lstep+1,l*lstep+1,xs,j,ind, &
                        worksg,iworksg)  !local copy of a standard routine
            coeflx(j,l) = xs(l*lstep+1)
            if     (mnproc.eq.1) then
              write(6,'(a,i2,i4,2f10.6)') &
                 'l,j,coef,*1+wb =', &
                  l,j,coeflx(j,l),coeflx(j,l)*(1.0+wbaro)
            endif
          enddo !j
        enddo !l
!
        deallocate( xm, xs, worksg, iworksg )
      else !.not.btrmas
        coeflx(:,:) = 0.0  !not used
      endif !btrmas:else
!
      return
      end subroutine barotp_init

      end module mod_barotp

!
!> Revision history:
!>
!> Mar. 1995 - changed vertical velocity averaging interval from 10 cm to 1 m
!>             (loops 33,35)
!> Mar. 1995 - changed order of loop nesting in loop 842
!> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
!> Aug. 1997 - transferred loops preceding loop 840 to momeq2.f
!> Jan. 2000 - added latbdp for lateral boundary ports
!> Aug. 2001 - two barotropic time steps per loop, for halo efficiency
!> Nov. 2006 - added lbflag==3 (latbdf) and svref_bt (mod_tides)
!> Nov. 2006 - removed svref_bt (and mod_tides)
!> Apr. 2007 - added btrlfr: leapfrog time step; see also momtum
!> Apr. 2010 - bugfixes for 1st time step and 1st miner time step
!> Apr  2011 - added    Robert-Asselin filtering for btrlfr
!> Aug  2011 - reworked Robert-Asselin filtering for btrlfr
!> Mar. 2012 - added latbdtf for nesting with Flather b.c.'s.
!> Jan. 2013 - added tidal drag tensor
!> June 2013 - added   lbflag==6 for latbdtc
!> Apr. 2014 - replace ip with ipa for mass sums
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> May  2014 - removed lbflag==6 for latbdtc
!> Apr. 2015 - added atmospheric pressure forcing
!> Apr. 2015 - added tidal body forcing
!> Aug. 2018 - added btrmas and barotp_init, converted to a module
!> Feb. 2019 - replaced onetai by 1.0
!> Sep. 2019 - added oneta0, and oneta diagnostic test
