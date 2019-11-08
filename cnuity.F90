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
      subroutine cnuity(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      use mod_pipe       ! HYCOM debugging interface
      use mod_floats     ! HYCOM synthetic floats, drifters and moorings
#if defined(STOKES)
      use mod_stokes     ! HYCOM Stokes drift
#endif
!
! --- hycom version 1.0
      implicit none
!
      integer m,n
!
! --- ------------------------------------------------------
! --- continuity equation (flux-corrected transport version)
! --- ------------------------------------------------------
! --- on entry:
! ---  dp( :,:,:,n) = time step t-1
! ---  dp( :,:,:,m) = time step t
!
! --- dpv( :,:,:,n) = time step t-1
! --- dpv( :,:,:,m) = time step t
! --- dpu( :,:,:,n) = time step t-1
! --- dpu( :,:,:,m) = time step t
!
! --- on exit:
! ---  dpo(:,:,:,n) = time step t-1
! ---  dpo(:,:,:,m) = time step t
! ---  dp( :,:,:,m) = time step t   with RA time smoothing
! ---  dp( :,:,:,n) = time step t+1
!
! --- onetacnt(:,:) = 1+eta afer cnuity
! --- ------------------------------------------------------
!
      logical, parameter :: lpipe_cnuity=.false. !usually .false.
!
      real,    parameter ::      dpfatal=-10.0   !fatal negative dp in meters
      real,    parameter :: epsil_cnuity=1.e-14
!
#if defined(RELO)
      integer, save, allocatable, dimension (:,:) :: &
       masku,maskv
      real,    save, allocatable, dimension (:,:) :: &
       pold,oneta_u,oneta_v
      real,    save, allocatable, dimension (:) :: &
       dpmn
#else
      integer, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       masku,maskv
      real,    save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       pold,oneta_u,oneta_v
      real,    save, dimension(1-nbdy:jdm+nbdy) :: &
       dpmn
#endif
!
      integer i,iflip,iprint,isave,j,jsave,k,l,ia,ib,ja,jb,margin,mbdy
      real    q,dpmin,clip,flxhi,flxlo,dtinv,dpup,dpdn,thkdfu,thkdfv
      real    dpold,dpmid,dpnew
      real    dpkmin(2*kdm)
!
      character*12 text,textu,textv
!
#if defined(RELO)
      if     (.not.allocated(masku)) then
        allocate( &
                  masku(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                  maskv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                   pold(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                oneta_u(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                oneta_v(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 4*(idm+2*nbdy)*(jdm+2*nbdy) ) !real=2*int
                  masku = -99
                  maskv = -99
                   pold = r_init
                oneta_u = r_init
                oneta_v = r_init
        allocate( &
                 dpmn(1-nbdy:jdm+nbdy) )
        call mem_stat_add( jdm+2*nbdy )
                 dpmn = r_init
      endif
!
#endif
      mbdy = 6
!
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1,   1, 6,6, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,1),1,2*kk, 6,6, halo_ps)
      call xctilr(dpu(    1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_us)
      call xctilr(dpv(    1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_vs)
      call xctilr(u(      1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_uv)
      call xctilr(v(      1-nbdy,1-nbdy,1,m),1,  kk, 6,6, halo_vv)
      call xctilr(ubavg(  1-nbdy,1-nbdy,  m),1,   1, 6,6, halo_uv)
      call xctilr(vbavg(  1-nbdy,1-nbdy,  m),1,   1, 6,6, halo_vv)
!
! --- rhs: dpmixl.n
! --- lhs: util3, utotn, vtotn
!
      margin = mbdy
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if     (btrmas) then
! ---       use dp, rather than dp'
            onetamas(i,j,:) = oneta(i,j,:)
            if (SEA_U) then
! ---         depthu is either pbot(i,j) or pbot(i-1,j)
              if     (pbot(i,j).eq.pbot(i-1,j)) then
                oneta_u(i,j) = 0.5*(onetamas(i,j,m)+onetamas(i-1,j,m))
              elseif (pbot(i,j).eq.depthu(i,j)) then
                oneta_u(i,j) =      onetamas(i,j,m)
              else
                oneta_u(i,j) =                      onetamas(i-1,j,m)
              endif
            endif !iu
            if (SEA_V) then
! ---         depthv is either pbot(i,j) or pbot(i,j-1)
              if     (pbot(i,j).eq.pbot(i,j-1)) then
                oneta_v(i,j) = 0.5*(onetamas(i,j,m)+onetamas(i,j-1,m))
              elseif (pbot(i,j).eq.depthv(i,j)) then
                oneta_v(i,j) =      onetamas(i,j,m)
              else
                oneta_v(i,j) =                      onetamas(i,j-1,m)
              endif
            endif !iv
          else
! ---       use dp'
            onetamas(i,j,:) = 1.0
             oneta_u(i,j)   = 1.0
             oneta_v(i,j)   = 1.0
          endif !btrmas:else
           utotn(i,j) = 0.0
           vtotn(i,j) = 0.0
           util3(i,j) = 0.0
          dpmold(i,j) = dpmixl(i,j,n)  ! save for Robert-Asselin filter
          do k= 1,kk
            dpo(i,j,k,n) = dp(i,j,k,n)  !t-1
          enddo !k
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      do 76 k=1,kk
!
! --- uflux/vflux = low-order (diffusive) mass fluxes at old time level.
! --- uflux2/vflux2 = 'antidiffusive' fluxes, defined as high-order minus low-
! --- order fluxes. high-order fluxes are second-order in space, time-centered.
!
! ---   rhs: depthu+, util3, dp.n, ubavg.m
! ---   lhs: uflux
!
        margin = mbdy - 1
!       write(6,*)'cnuity.F line 96'
!
        do j=1-margin,jj+margin
          do l=1,isu(j) !ok
            i=ifu(j,l)-1
            if (i.ge.1-margin) then
              if (iuopn(i,j).ne.0) then
                q=min(dp(i  ,j,k,n),max(0.,depthu(i+1,j)-util3(i  ,j))) &
                    * onetamas(i  ,j,n)
                utotm(i,j)=(u(i+1,j,k,m)+ubavg(i,j,m))*scuy(i,j)
#if defined(STOKES)
                utotm(i,j)=utotm(i,j)+usd(i+1,j,k)*scuy(i,j)
#endif                
                uflux(i,j)=utotm(i,j)*q
              endif
            endif
            i=ilu(j,l)+1
            if (i.le.ii+margin) then
              if (iuopn(i,j).ne.0) then
                q=min(dp(i-1,j,k,n),max(0.,depthu(i-1,j)-util3(i-1,j))) &
                       * onetamas(i-1,j,n)
                utotm(i,j)=(u(i-1,j,k,m)+ubavg(i,j,m))*scuy(i,j)
#if defined(STOKES)
                utotm(i,j)=utotm(i,j)+usd(i-1,j,k)*scuy(i,j)
#endif                
                uflux(i,j)=utotm(i,j)*q
              endif
            endif
          enddo
        enddo
!        write(6,*)'cnuity.F line 128'
!
! ---   rhs: depthv+, util3, dp.n, vbavg.m
! ---   lhs: vflux
!
        margin = mbdy - 1
!
        do i=1-margin,ii+margin
          do l=1,jsv(i) !ok
            j=jfv(i,l)-1
            if (j.ge.1-margin) then
              if (ivopn(i,j).ne.0) then
                q=min(dp(i,j  ,k,n),max(0.,depthv(i,j+1)-util3(i,j  ))) &
                      * onetamas(i,j  ,n)
                vtotm(i,j)=(v(i,j+1,k,m)+vbavg(i,j,m))*scvx(i,j)
#if defined(STOKES)
                vtotm(i,j)=vtotm(i,j)+vsd(i,j+1,k)*scvx(i,j)
#endif                
                vflux(i,j)=vtotm(i,j)*q
              endif
            endif
            j=jlv(i,l)+1
            if (j.le.jj+margin) then
              if (ivopn(i,j).ne.0) then
                q=min(dp(i,j-1,k,n),max(0.,depthv(i,j-1)-util3(i,j-1))) &
                         * onetamas(i,j-1,n)
                vtotm(i,j)=(v(i,j-1,k,m)+vbavg(i,j,m))*scvx(i,j)
#if defined(STOKES)
                vtotm(i,j)=vtotm(i,j)+vsd(i,j-1,k)*scvx(i,j)
#endif                
                vflux(i,j)=vtotm(i,j)*q
              endif
            endif
          enddo
        enddo
!      write(6,*)'cnuity.F line 165'
!
! ---   rhs: u.m, ubavg.m, depthu, dp.n+, util3+, dpu.m, uflux
! ---   rhs: v.m, vbavg.m, depthv, dp.n+, util3+, dpv.m, vflux
! ---   lhs: utotm,uflux,uflux2,uflx
! ---   lhs: vtotm,vflux,vflux2,vflx
!
        margin = mbdy - 1
!
!$OMP   PARALLEL DO PRIVATE(j,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              utotm(i,j)=(u(i,j,k,m)+ubavg(i,j,m))*scuy(i,j)
#if defined(STOKES)
              utotm(i,j)=utotm(i,j)+usd(i,j,k)*scuy(i,j)
#endif              
              if (utotm(i,j).ge.0.) then
                q=min(dp(i-1,j,k,n),max(0.,depthu(i,j)-util3(i-1,j))) &
                   * onetamas(i-1,j,n)
              else
                q=min(dp(i  ,j,k,n),max(0.,depthu(i,j)-util3(i  ,j))) &
                   * onetamas(i  ,j,n)
              endif
               uflux(i,j)=utotm(i,j)*q
              uflux2(i,j)=utotm(i,j)*dpu(i,j,k,m)*oneta_u(i,j)- &
                                                    uflux(i,j)
              uflx(i,j,k)=uflux(i,j)
            endif !iu
          enddo !i
!
!      write(6,*)'cnuity.F line 197, j = ',j
          do i=1-margin,ii+margin
            if (SEA_V) then
              vtotm(i,j)=(v(i,j,k,m)+vbavg(i,j,m))*scvx(i,j)
#if defined(STOKES)
              vtotm(i,j)=vtotm(i,j)+vsd(i,j,k)*scvx(i,j)
#endif              
              if (vtotm(i,j).ge.0.) then
                q=min(dp(i,j-1,k,n),max(0.,depthv(i,j)-util3(i,j-1))) &
                   * onetamas(i,j-1,n)
              else
                q=min(dp(i,j  ,k,n),max(0.,depthv(i,j)-util3(i,j  ))) &
                   * onetamas(i,j  ,n)
              endif
               vflux(i,j)=vtotm(i,j)*q
              vflux2(i,j)=vtotm(i,j)*dpv(i,j,k,m)*oneta_v(i,j)- &
                                                    vflux(i,j)
              vflx(i,j,k)=vflux(i,j)
            endif !iv
          enddo !i
!        write(6,*)'cnuity.F line 216, j = ',j
        enddo !j
!$OMP   END PARALLEL DO
!       write(6,*)'cnuity.F line 216'
!
! ---   advance -dp- field using low-order (diffusive) flux values
! ---   rhs: dp.n, dp.m, util3, uflux+, vflux+
! ---   lhs: dpo,util3,dp.n
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i,dpmin) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          dpmin=999.
          do i=1-margin,ii+margin
            if (SEA_P) then
              util3(i,j)=util3(i,j)+dp(i,j,k,n)
              dp(i,j,k,n)=dp(i,j,k,n)* &
                           onetamas(i,j,n)- &
                          ((uflux(i+1,j)-uflux(i,j))+ &
                           (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
              dpo(i,j,k,m)=dp(i,j,k,n)  ! save for loop 19 test
              dpmin=min(dpmin,dp(i,j,k,n))
            endif !ip
          enddo !i
          dpmn(j)=dpmin   ! minimizes false sharing
        enddo  !j, loop 19
!$OMP   END PARALLEL DO
!
        dpmin=999.
        do j=1,jj
          dpmin=min(dpmin,dpmn(j))
        enddo
        dpkmin(k)=dpmin
!
        if     (lpipe .and. lpipe_cnuity) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'dp.low k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        endif
!
        do j=1-margin,jj+margin
          do l=1,isu(j) !ok
            i=ifu(j,l)-1
            if (i.ge.1-margin) then
              if (iuopn(i,j).ne.0) then
                uflux(i,j)=0.0
              endif
            endif
            i=ilu(j,l)+1
            if (i.le.ii+margin) then
              if (iuopn(i,j).ne.0) then
                uflux(i,j)=0.0
              endif
            endif
          enddo
        enddo
!
        do i=1-margin,ii+margin
          do l=1,jsv(i) !ok
            j=jfv(i,l)-1
            if (j.ge.1-margin) then
              if (ivopn(i,j).ne.0) then
                vflux(i,j)=0.0
              endif
            endif
            j=jlv(i,l)+1
            if (j.le.jj+margin) then
              if (ivopn(i,j).ne.0) then
                vflux(i,j)=0.0
              endif
            endif
          enddo
        enddo
!
        if     (lpipe .and. lpipe_cnuity) then
! ---     compare two model runs.
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              masku(i,j)=iu(i,j)
              if (i.gt. 1) masku(i,j)=masku(i,j)+iu(i-1,j)
              if (i.lt.ii) masku(i,j)=masku(i,j)+iu(i+1,j)
              maskv(i,j)=iv(i,j)
              if (j.gt. 1) maskv(i,j)=maskv(i,j)+iv(i,j-1)
              if (j.lt.jj) maskv(i,j)=maskv(i,j)+iv(i,j+1)
            enddo
          enddo
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux,masku,textu, &
                                 vflux,maskv,textv)
          write (textu,'(a9,i3)') 'uflux2 k=',k
          write (textv,'(a9,i3)') 'vflux2 k=',k
          call pipe_compare_sym2(uflux2,masku,textu, &
                                 vflux2,maskv,textv)
        endif
!
!diag if (mod(k,15).eq.1) then
!diag   do i=itest-1,itest+1
!diag   do j=jtest-1,jtest+1
!diag   write (lp,101) nstep,i+i0,j+j0,k,dpo(i-1,j,k,n),uflux(i,j), &
!diag    'old dp''s, fluxes:',dpo(i,j-1,k,n),dpo(i,j,k,n),dpo(i,j+1,k,n) &
!diag    ,vflux(i,j),dp(i,j,k,n),vflux(i,j+1),dpo(i+1,j,k,n),uflux(i+1,j)
!diag   enddo
!diag   enddo
!diag endif
 101  format (i9,2i5,i3,1p,e15.2,e30.2/a17,6e10.2/e37.2,e30.2)
!
! --- at each grid point, determine the ratio of the largest permissible
! --- pos. (neg.) change in -dp- to the sum of all incoming (outgoing) fluxes
!
! ---   rhs: dp.n+, uflux2+, vflux2+
! ---   lhs: util1,util2
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
! ---         assume margin<nblk
              ia=ipim1(i,j) !i-1 if sea and i if land
              ib=ipip1(i,j) !i+1 if sea and i if land
              ja=ipjm1(i,j) !j-1 if sea and j if land
              jb=ipjp1(i,j) !j+1 if sea and j if land
              util1(i,j)=max(dp(i,j,k,n),dp(ia,j,k,n),dp(ib,j,k,n), &
                                         dp(i,ja,k,n),dp(i,jb,k,n))
              util2(i,j)=max(0., &
                         min(dp(i,j,k,n),dp(ia,j,k,n),dp(ib,j,k,n), &
                                         dp(i,ja,k,n),dp(i,jb,k,n)))
!
              util1(i,j)=(util1(i,j)-dp(i,j,k,n)) &
              /(((max(0.,uflux2(i,j))-min(0.,uflux2(i+1,j))) &
                +(max(0.,vflux2(i,j))-min(0.,vflux2(i,j+1)))+epsil) &
              *delt1*scp2i(i,j))
!
              util2(i,j)=(util2(i,j)-dp(i,j,k,n)) &
              /(((min(0.,uflux2(i,j))-max(0.,uflux2(i+1,j))) &
                +(min(0.,vflux2(i,j))-max(0.,vflux2(i,j+1)))-epsil) &
              *delt1*scp2i(i,j))
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! --- limit antidiffusive fluxes
! --- (keep track in -utotn,vtotn- of discrepancy between high-order
! --- fluxes and the sum of low-order and clipped antidiffusive fluxes.
! --- this will be used later to restore nondivergence of barotropic flow)
!
! ---   rhs: uflux2+,util1+,util2+
! ---   rhs: vflux2+,util1+,util2+
! ---   lhs: utotn,uflux,uflx
! ---   lhs: vtotn,vflux,vflx
!
        margin = mbdy - 3
!
!$OMP   PARALLEL DO PRIVATE(j,i,clip) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              if (uflux2(i,j).ge.0.) then
                clip=min(1.,util1(i,j),util2(i-1,j))
              else
                clip=min(1.,util2(i,j),util1(i-1,j))
              endif
              utotn(i,j)=utotn(i,j)+uflux2(i,j)*(1.-clip)
              uflux(i,j)=uflux2(i,j)*clip
              uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
            endif !iu
          enddo !i
!
          do i=1-margin,ii+margin
            if (SEA_V) then
              if (vflux2(i,j).ge.0.) then
                clip=min(1.,util1(i,j),util2(i,j-1))
              else
                clip=min(1.,util2(i,j),util1(i,j-1))
              endif
              vtotn(i,j)=vtotn(i,j)+vflux2(i,j)*(1.-clip)
              vflux(i,j)=vflux2(i,j)*clip
              vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! --- evaluate effect of antidiffusive fluxes on -dp- field
!
! ---   rhs: dp.n, p, uflux+,vflux+
! ---   lhs: dp.n, p
!
        margin = mbdy - 4
!
!$OMP   PARALLEL DO PRIVATE(j,i,dpmin) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          dpmin=999.
          do i=1-margin,ii+margin
            if (SEA_P) then
              dp(i,j,k,n)=dp(i,j,k,n)- &
                           ((uflux(i+1,j)-uflux(i,j))+ &
                            (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              dpmin=min(dpmin,dp(i,j,k,n))
            endif !ip
          enddo !i
          dpmn(j)=dpmin   ! minimizes false sharing
        enddo  !j, loop 15
!$OMP   END PARALLEL DO
!
        dpmin=999.
        do j=1,jj
          dpmin=min(dpmin,dpmn(j))
        enddo
        dpkmin(k+kk)=dpmin
!
        if     (lpipe .and. lpipe_cnuity) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'dp-dif k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        endif
!
 76   continue  ! k=1,kk
!
! --- check for negative thicknesses.
!
 100  format (i9,' i,j,k=',2i5,i3,' neg. dp (m) in loop ',i2,g15.2,a)
!
      if     (mod(nstep,3).eq.0) then  !skip some time steps for efficiency
        call xcminr(dpkmin(1:2*kk))
        do k= 1,kk
          dpmin=dpkmin(k)
          if ((k.eq.1 .and. dpmin.le.0.0) .or. dpmin.lt.-onecm) then
            iprint = 0
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if (dpo(i,j,k,m).eq.dpmin .and. iprint.le.5) then
                    write (lp,100) nstep,i+i0,j+j0,k,19, &
                                   dpmin*qonem,' max'
                    iprint = iprint + 1
                  endif
                endif !ip
              enddo !i
            enddo !j
            call xcsync(flush_lp)
            if     (dpmin.lt.dpfatal*onem) then
              do j=1,jj
                do i=1,ii
                  if (SEA_P) then
                    if (dpo(i,j,k,m).lt.dpfatal*onem .and. &
                        iprint.le.5) then
                      write (lp,100) nstep,i+i0,j+j0,k,19, &
                                     dpo(i,j,k,m)*qonem,' fatal'
                      iprint = iprint + 1
                    endif
                  endif !ip
                enddo !i
              enddo !j
              call xcsync(flush_lp)
            endif !dpfatal
          endif
        enddo !k
        do k= 1,kk
          dpmin=dpkmin(k+kk)
          if ((k.eq.1 .and. dpmin.le.0.0) .or. dpmin.lt.-onecm) then
            iprint = 0
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if (dp(i,j,k,n).eq.dpmin .and. iprint.le.5) then
                    write (lp,100) nstep,i+i0,j+j0,k,15, &
                                   dpmin*qonem,' max'
                    iprint = iprint + 1
                  endif
                endif !ip
              enddo !i
            enddo !j
            call xcsync(flush_lp)
            if     (dpmin.lt.dpfatal*onem) then
              do j=1,jj
                do i=1,ii
                  if (SEA_P) then
                    if (dp(i,j,k,n).lt.dpfatal*onem .and. &
                        iprint.le.5) then
                      write (lp,100) nstep,i+i0,j+j0,k,15, &
                                     dp(i,j,k,n)*qonem,' fatal'
                      iprint = iprint + 1
                    endif
                  endif !ip
                enddo !i
              enddo !j
              call xcsync(flush_lp)
            endif !dpfatal
          endif
        enddo !k
        if     (minval(dpkmin(1:2*kk)).lt.dpfatal*onem) then
          if     (mnproc.eq.1) then
            write(lp,'(/ a,f9.2 /)') &
              'error: neg. dp (m) < ',dpfatal
!           write(lp,*) 'dpkmin =',qonem*dpkmin(1:2*kk)
          endif
          call xcstop('cnuity')
                 stop 'cnuity'
        endif !dpfatal
      endif !every 3 time steps
!
      if     (.not.btrmas) then
!
! --- restore nondivergence of vertically integrated mass flow by
! --- recovering fluxes lost in the flux limiting process.
! --- treat these fluxes as an 'upstream' barotropic correction to
! --- the sum of diffusive and antidiffusive fluxes obtained so far.
!
      do 77 k=1,kk
!
! ---   rhs: utotn, dp.n+, p+
! ---   rhs: vtotn, dp.n+, p+
! ---   lhs: uflux, uflx
! ---   lhs: vflux, vflx
!
        margin = mbdy - 5
!
!$OMP   PARALLEL DO PRIVATE(j,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              if (utotn(i,j).ge.0.) then
!               if     (p(i-1,j,kk+1).eq.0.0) then
!                 write(lp,*) 'error: i,j,p.i-1 = ',
!    &                        i,j,p(i-1,j,kk+1)
!                 call xcstop('(cnuity)')
!                        stop
!               endif
                q=dp(i-1,j,k,n)/p(i-1,j,kk+1)
              else
!               if     (p(i,  j,kk+1).eq.0.0) then
!                 write(lp,*) 'error: i,j,p.i   = ',
!    &                        i,j,p(i,  j,kk+1)
!                 call xcstop('(cnuity)')
!                        stop
!               endif
                q=dp(i  ,j,k,n)/p(i  ,j,kk+1)
              endif
              uflux(i,j)=utotn(i,j)*q
              uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
            endif !iu
          enddo !i
!
          do i=1-margin,ii+margin
            if (SEA_V) then
              if (vtotn(i,j).ge.0.) then
!               if     (p(i,j-1,kk+1).eq.0.0) then
!                 write(lp,*) 'error: i,j,p.j-1 = ',
!    &                        i,j,p(i,j-1,kk+1)
!                 call xcstop('(cnuity)')
!                        stop
!               endif
                q=dp(i,j-1,k,n)/p(i,j-1,kk+1)
              else
!               if     (p(i,j,  kk+1).eq.0.0) then
!                 write(lp,*) 'error: i,j,p.j   = ',
!    &                        i,j,p(i,j,  kk+1)
!                 call xcstop('(cnuity)')
!                        stop
!               endif
                q=dp(i,j  ,k,n)/p(i,j  ,kk+1)
              endif
              vflux(i,j)=vtotn(i,j)*q
              vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! ---   rhs: dp.n, p, uflux+, vflux+
! ---   lhs: dp.n, p
!
        margin = mbdy - 6
!
!$OMP   PARALLEL DO PRIVATE(j,i,dpmin) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          dpmin=999.
          do i=1-margin,ii+margin
            if (SEA_P) then
              dp( i,j,k,n)=dp(i,j,k,n)- &
                           ((uflux(i+1,j)-uflux(i,j))+ &
                            (vflux(i,j+1)-vflux(i,j)))*delt1*scp2i(i,j)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
              dpmin=min(dpmin,dp(i,j,k,n))
            endif !ip
          enddo !i
          dpmn(j)=dpmin   ! minimizes false sharing
        enddo  !j, loop 14
!
        dpmin=999.
        do j=1,jj
          dpmin=min(dpmin,dpmn(j))
        enddo
        dpkmin(k)=dpmin
!
        if     (lpipe .and. lpipe_cnuity) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'dp.res k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        endif
!
 77   continue  ! k=1,kk
!
      call pipe_comparall(m,n, 'cnui77, step')
!
      endif !.not.btrmas
!
! --- check for negative thicknesses.
!
      if     (mod(nstep,3).eq.0) then  !skip some time steps for efficiency
        call xcminr(dpkmin(1:kk))
        do k= 1,kk
          dpmin=dpkmin(k)
          if (dpmin.lt.-onecm) then
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if (dp(i,j,k,n).eq.dpmin) then
                    write (lp,100) nstep,i+i0,j+j0,k,14, &
                                   dpmin*qonem,' max'
                  endif
                endif !ip
              enddo !i
            enddo !j
            call xcsync(flush_lp)
          endif
        enddo !k
      endif !every 3 time steps
!
      if  (.not.btrmas) then
!
! --- add bottom-pressure restoring term arising from split-explicit treatment
! --- of continuity equation (step 4 in appendix b to 1992 brhs paper)
!
! ---   rhs: dp.n, p, pbot
! ---   lhs: dp.n, p, dpmixl.n
!
        margin = mbdy - 6
!
!$OMP PARALLEL DO PRIVATE(j,k,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            q = pbot(i,j)/p(i,j,kk+1)
            do k=1,kk
              dp(i,j,k,n)=dp(i,j,k,n)*q
              p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
            enddo !k
            if (isopyc) then
              dpmixl(i,j,n)=dp(i,j,1,n)
            endif
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_cnuity) then
! ---   compare two model runs.
        do k=1,kk
          write (text,'(a9,i3)') 'dp.bot k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        enddo
      endif
!
      endif  !.not.btrmas
!
! --- ---------------------------------------------------------------------
! --- biharmonic thickness diffusion (literally, interface depth diffusion)
! --- ---------------------------------------------------------------------
!
      if (thkdf4.eq.0.) go to 800  ! only one of thkdf2 and thkdf4 is non-zero
                                   ! apply thkdf4
!
      mbdy = 6
!
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1, 1, 6,6, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
      call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 6,6, halo_ps)
!
      dtinv=1./delt1
      iflip=mod(nstep,2)
!
! --- rhs: p
! --- lhs: uflux, vflux, pold
!
      margin = mbdy
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            uflux(i,j)=0.
          endif !iu
          if (SEA_V) then
            vflux(i,j)=0.
          endif !iv
          if (SEA_P) then
            if (iflip.eq.1) then
              pold(i,j)=p(i,j,kk+1)
            else
              pold(i,j)=0.
            endif
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
! --- alternate between upward and downward direction in k loop
!
      do 13 k=2*(1-iflip)+kk*iflip,kk*(1-iflip)+2*iflip,1-2*iflip
!
! ---   rhs: p+, dp.n+
! ---   lhs: util1, util2
!
        margin = mbdy - 1
!
!$OMP   PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              if (min(dp(i,j,k,n),dp(i,j,k-1,n)).lt.onecm) then
                util1(i,j)=0.0
                util2(i,j)=0.0
              else
                ia=ipim1x(i,j) !i-1 if sea; else i+1 if sea; otherwise i
                ib=ipip1x(i,j) !i+1 if sea; else i-1 if sea; otherwise i
                ja=ipjm1x(i,j) !j-1 if sea; else j+1 if sea; otherwise j
                jb=ipjp1x(i,j) !j+1 if sea; else j-1 if sea; otherwise j
                util1(i,j)=p(i,j,k)-.5*(p(ia,j,k)+p(ib,j,k))
                util2(i,j)=p(i,j,k)-.5*(p(i,ja,k)+p(i,jb,k))
                if (util1(i,j).gt.0.0) then
                  if (min(dp(ia,j,k,  n),dp(ib,j,k,  n)).lt.onecm) then
                    util1(i,j)=0.0
                  endif
                else
                  if (min(dp(ia,j,k-1,n),dp(ib,j,k-1,n)).lt.onecm) then
                    util1(i,j)=0.0
                  endif
                endif
                if (util2(i,j).gt.0.0) then
                  if (min(dp(i,ja,k,  n),dp(i,jb,k,  n)).lt.onecm) then
                    util2(i,j)=0.0
                  endif
                else
                  if (min(dp(i,ja,k-1,n),dp(i,jb,k-1,n)).lt.onecm) then
                    util2(i,j)=0.0
                  endif
                endif
              endif !min(dp):else
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! ---   limit fluxes to avoid intertwining interfaces
!
! ---   rhs: p+, pold+, uflux, util1+, uflx+
! ---   rhs: p+, pold+, vflux, util2+, vflx+
! ---   lhs: uflux, uflx+
! ---   lhs: vflux, vflx+
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i,flxhi,flxlo,thkdfu,thkdfv) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
!
          do i=1-margin,ii+margin
            if (SEA_U) then
              flxhi= .25*(p(i  ,j,kk+1)-p(i  ,j,k))*scp2(i  ,j)
              flxlo=-.25*(p(i-1,j,kk+1)-p(i-1,j,k))*scp2(i-1,j)
!
              if (iflip.eq.0) then		!  downward k loop
                flxhi=min(flxhi, &
                          uflux(i,j)+ &
                            .25*(p(i-1,j,k)-pold(i-1,j))*scp2(i-1,j))
                flxlo=max(flxlo, &
                          uflux(i,j)- &
                            .25*(p(i  ,j,k)-pold(i  ,j))*scp2(i  ,j))
              else				!  upward k loop
                flxhi=min(flxhi, &
                          uflux(i,j)+ &
                            .25*(pold(i  ,j)-p(i  ,j,k))*scp2(i  ,j))
                flxlo=max(flxlo, &
                          uflux(i,j)- &
                            .25*(pold(i-1,j)-p(i-1,j,k))*scp2(i-1,j))
              endif
!
              thkdfu = delt1*thkdf4u(i,j)
              uflux(i,j)=min( flxhi, &
                              max( flxlo, &
                                   thkdfu* &
                                     (util1(i-1,j)-util1(i,j)) ) )
              uflx(i,j,k-1)=uflx(i,j,k-1)+uflux(i,j)*dtinv
              uflx(i,j,k  )=uflx(i,j,k  )-uflux(i,j)*dtinv
            endif !iu
          enddo !i
!
          do i=1-margin,ii+margin
            if (SEA_V) then
              flxhi= .25*(p(i,j  ,kk+1)-p(i,j  ,k))*scp2(i,j  )
              flxlo=-.25*(p(i,j-1,kk+1)-p(i,j-1,k))*scp2(i,j-1)
!
              if (iflip.eq.0) then              !  downward k loop
                flxhi=min(flxhi, &
                          vflux(i,j)+ &
                            .25*(p(i,j-1,k)-pold(i,j-1))*scp2(i,j-1))
                flxlo=max(flxlo, &
                          vflux(i,j)- &
                            .25*(p(i,j  ,k)-pold(i,j  ))*scp2(i,j  ))
              else                              !  upward k loop
                flxhi=min(flxhi, &
                          vflux(i,j)+ &
                            .25*(pold(i,j  )-p(i,j  ,k))*scp2(i,j  ))
                flxlo=max(flxlo, &
                          vflux(i,j)- &
                            .25*(pold(i,j-1)-p(i,j-1,k))*scp2(i,j-1))
              endif
!
              thkdfv = delt1*thkdf4v(i,j)
              vflux(i,j)=min( flxhi, &
                              max( flxlo, &
                                   thkdfv* &
                                     (util2(i,j-1)-util2(i,j)) ) )
              vflx(i,j,k-1)=vflx(i,j,k-1)+vflux(i,j)*dtinv
              vflx(i,j,k  )=vflx(i,j,k  )-vflux(i,j)*dtinv
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! ---   rhs: p, uflux+, vflux+
! ---   lhs: pold, p
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pold(i,j)=p(i,j,k)
              p(i,j,k)=p(i,j,k)-((uflux(i+1,j)-uflux(i,j))+ &
                                 (vflux(i,j+1)-vflux(i,j)))*scp2i(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
!diag if (itest.gt.0.and.jtest.gt.0) &
!diag write (lp,'(i9,2i5,i3,'' intfc.depth diffusion -- p_old,p_new ='', &
!diag 2f9.3)') nstep,itest,jtest,k,pold(itest,jtest)*qonem,p(itest, &
!diag jtest,k)*qonem
!
 13   continue  ! k=2*(1-iflip)+kk*iflip,kk*(1-iflip)+2*iflip,1-2*iflip
!
! --- rhs: p
! --- lhs: p, dp.n, dpmixl.n
!
      margin = mbdy - 2
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do k=1,kk
          do i=1-margin,ii+margin
            if (SEA_P) then
              if (p(i,j,k+1).lt.p(i,j,k)) then
!diag           write (lp,'(i9,2i5,i3,a,g15.2,i4)') nstep,i+i0,j+j0,k, &
!diag           '  neg. dp after thknss smoothing', &
!diag           qonem*(p(i,j,k+1)-p(i,j,k)),iflip
                p(i,j,k+1)=p(i,j,k)
              endif
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
              if (isopyc .and. k.eq.1) then
                dpmixl(i,j,n)=dp(i,j,k,n)
              endif
            endif !ip
          enddo !i
        enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_cnuity) then
! ---   compare two model runs.
        do k=1,kk
          write (text,'(a9,i3)') 'dp-df4 k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        enddo
      endif
!
 800  continue  ! end of biharmonic thickness diffusion
!
! --- ---------------------------------------------------------------------
! --- Laplacian thickness diffusion (literally, interface depth diffusion)
! --- ---------------------------------------------------------------------
!
      if (thkdf2.eq.0.) go to 850  ! only one of thkdf2 and thkdf4 is non-zero
                                   ! apply thkdf2
!
      mbdy = 6
!
      call xctilr(dpmixl( 1-nbdy,1-nbdy,  n),1, 1, 6,6, halo_ps)
      call xctilr(dp(     1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
      call xctilr(p(      1-nbdy,1-nbdy,2  ),1,kk, 6,6, halo_ps)
!
      dtinv=1./delt1
!
! --- rhs: p
! --- lhs: uflux, vflux, pold
!
      margin = mbdy
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            uflux(i,j)=0.
          endif !iu
          if (SEA_V) then
            vflux(i,j)=0.
          endif !iv
          if (SEA_P) then
            vflux(i,j)=0.
            if (iflip.eq.1) then
              pold(i,j)=p(i,j,kk+1)
            else
              pold(i,j)=0.
            endif
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      do 113 k=2,kk
!
! ---   rhs: p+, dp.n+
! ---   lhs: util1, util2
!
        margin = mbdy - 1
!
! ---   limit fluxes to avoid intertwining interfaces
!
! ---   rhs: p+, pold+, uflux, util1+, uflx+
! ---   rhs: p+, pold+, vflux, util2+, vflx+
! ---   lhs: uflux, uflx+
! ---   lhs: vflux, vflx+
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i,flxhi,flxlo,thkdfu,thkdfv) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
!
          do i=1-margin,ii+margin
            if (SEA_U) then
              flxhi= .25*(p(i  ,j,kk+1)-p(i  ,j,k))*scp2(i  ,j)
              flxlo=-.25*(p(i-1,j,kk+1)-p(i-1,j,k))*scp2(i-1,j)
              thkdfu = delt1*thkdf4u(i,j) !thkdf2u
              uflux(i,j)=min(flxhi, &
                             max(flxlo, &
                                 thkdfu*(p(i-1,j,k)-p(i,j,k)) ))
              uflx(i,j,k-1)=uflx(i,j,k-1)+uflux(i,j)*dtinv
              uflx(i,j,k  )=uflx(i,j,k  )-uflux(i,j)*dtinv
            endif !iu
          enddo !i
!
          do i=1-margin,ii+margin
            if (SEA_V) then
              flxhi= .25*(p(i,j  ,kk+1)-p(i,j  ,k))*scp2(i,j  )
              flxlo=-.25*(p(i,j-1,kk+1)-p(i,j-1,k))*scp2(i,j-1)
              thkdfv = delt1*thkdf4v(i,j) !thkdf2v
              vflux(i,j)=min(flxhi, &
                             max(flxlo, &
                                 thkdfv*(p(i,j-1,k)-p(i,j,k)) ))
              vflx(i,j,k-1)=vflx(i,j,k-1)+vflux(i,j)*dtinv
              vflx(i,j,k  )=vflx(i,j,k  )-vflux(i,j)*dtinv
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
! ---   rhs: p, uflux+, vflux+
! ---   lhs: pold, p
!
        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pold(i,j)=p(i,j,k)
              p(i,j,k) =p(i,j,k)-((uflux(i+1,j)-uflux(i,j))+ &
                                  (vflux(i,j+1)-vflux(i,j)))*scp2i(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
!diag if (itest.gt.0.and.jtest.gt.0) &
!diag write (lp,'(i9,2i5,i3," intfc.depth diffusion -- p_old,p_new =", &
!diag 2f9.3)') nstep,itest+i0,jtest+j0,k,pold(itest,jtest)*qonem, &
!diag p(itest,jtest,k)*qonem
!
 113  continue  ! k=2,kk
!
! --- rhs: p
! --- lhs: p, dp.n, dpmixl.n
!
      margin = mbdy - 2
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do k=1,kk
          do i=1-margin,ii+margin
            if (SEA_P) then
              if (p(i,j,k+1).lt.p(i,j,k)) then
!diag           write (lp,'(i9,2i5,i3,a,g15.2,i4)') nstep,i+i0,j+j0,k, &
!diag           '  neg. dp after thknss smoothing', &
!diag           qonem*(p(i,j,k+1)-p(i,j,k)),iflip
                p(i,j,k+1)=p(i,j,k)
              endif
              dp(i,j,k,n)=p(i,j,k+1)-p(i,j,k)
              if (isopyc .and. k.eq.1) then
                dpmixl(i,j,n)=dp(i,j,k,n)
              endif
            endif !ip
          enddo !i
        enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_cnuity) then
! ---   compare two model runs.
        do k=1,kk
          write (text,'(a9,i3)') 'dp-df2 k=',k
          call pipe_compare_sym1(dp(1-nbdy,1-nbdy,k,n),ip,text)
        enddo
      endif
!
 850  continue  ! end of Laplacian thickness diffusion
!
! --- estimate wveli for synthetic floats as the change in pressure
! --- interface depth in meters, positive upward
      if (synflt .and. wvelfl) then
!$OMP   PARALLEL DO PRIVATE(j,i,k) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              pold(i,j)=0.0
              do k=1,kk
                pold(i,j)=pold(i,j)+dpo(i,j,k,n)
                wveli(i,j,k+1)=-(p(i,j,k+1)-pold(i,j))/onem
              enddo !k
            endif !ip
          enddo !i
        enddo !j
      endif !synflt+wvelfl
!
! --- account for vertical advection of dpmixl. calculate the vertical
! --- excursions of the coordinates immediately above and below the mixed
! --- layer base, then vertically interpolate this motion to dpmixl.
! --- also apply biharmonic thickness diffusion to the mixed layer.
      if(hybrid .and. mxlkta) then
!
! ---   rhs: util1, util2, dpmixl.n, dpo, p
! ---   lhs: util1, util2, dpmixl.n

        margin = mbdy - 2
!
!$OMP   PARALLEL DO PRIVATE(j,k,i,dpup,dpdn,q) &
!$OMP          SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              util1(i,j)=0.
              util2(i,j)=0.
            endif !ip
          enddo !i
          do k=1,kk
            do i=1-margin,ii+margin
              if (SEA_P) then
                util1(i,j)=util2(i,j)
                util2(i,j)=util2(i,j)+dpo(i,j,k,n)
                if (util2(i,j).ge.dpmixl(i,j,n).and. &
                    util1(i,j).lt.dpmixl(i,j,n)     ) then
                  dpup=p(i,j,k  )-util1(i,j)
                  dpdn=p(i,j,k+1)-util2(i,j)
                  q=(util2(i,j)-dpmixl(i,j,n))/max(onemm,dpo(i,j,k,n))
                  dpmixl(i,j,n)=dpmixl(i,j,n)+(dpdn+q*(dpup-dpdn))
                endif
              endif !ip
            enddo !i
          enddo !k
        enddo !j
!$OMP   END PARALLEL DO
!
        if (thkdf4.ne.0.) then
!
! ---     lhs: uflux, vflux
!
          margin = mbdy - 3
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_U) then
                uflux(i,j)=0.
              endif !iu
              if (SEA_V) then
                vflux(i,j)=0.
              endif !iu
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
! ---     rhs: dpmixl.n+
! ---     lhs: util1, util2
!
          margin = mbdy - 4
!
!$OMP     PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_P) then
                ia=ipim1x(i,j) !i-1 if sea; else i+1 if sea; otherwise i
                ib=ipip1x(i,j) !i+1 if sea; else i-1 if sea; otherwise i
                ja=ipjm1x(i,j) !j-1 if sea; else j+1 if sea; otherwise j
                jb=ipjp1x(i,j) !j+1 if sea; else j-1 if sea; otherwise j
                util1(i,j)=dpmixl(i,j,n)- &
                            0.5*(dpmixl(ia,j,n)+dpmixl(ib,j,n))
                util2(i,j)=dpmixl(i,j,n)- &
                            0.5*(dpmixl(i,ja,n)+dpmixl(i,jb,n))
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
! ---     rhs: util1+, util2+
! ---     lhs: uflux, vflux
!
          margin = mbdy - 5
!
!$OMP     PARALLEL DO PRIVATE(j,i,thkdfu,thkdfv) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_U) then
                thkdfu = delt1*thkdf4u(i,j)
                uflux(i,j)=thkdfu*(util1(i-1,j)-util1(i,j))
              endif !iu
              if (SEA_V) then
                thkdfv = delt1*thkdf4v(i,j)
                vflux(i,j)=thkdfv*(util2(i,j-1)-util2(i,j))
              endif !iv
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
! ---     rhs: dpmixl.n, uflux+, vflux+
! ---     lhs: dpmixl.n
!
          margin = mbdy - 6
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_P) then
                dpmixl(i,j,n)=dpmixl(i,j,n)- &
                               ((uflux(i+1,j)-uflux(i,j))+ &
                                (vflux(i,j+1)-vflux(i,j)))*scp2i(i,j)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
        endif  ! end: (thkdf4.ne.0.)
!
        if (thkdf2.ne.0.) then
          margin = mbdy-3
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_U) then
                uflux(i,j)=0.
              endif !iu
              if (SEA_V) then
                vflux(i,j)=0.
              endif !iu
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
! ---     rhs: dpmixl.n+
! ---     lhs: util1, util2
!
          margin = mbdy - 4
!
!$OMP     PARALLEL DO PRIVATE(j,i,thkdfu,thkdfv) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_U) then
                thkdfu = delt1*thkdf4u(i,j) !thkdf2u
                uflux(i,j)=thkdfu*(dpmixl(i-1,j,n)-dpmixl(i,j,n))
              endif !iu
              if (SEA_V) then
                thkdfv = delt1*thkdf4v(i,j) !thkdf2v
                vflux(i,j)=thkdfv*(dpmixl(i,j-1,n)-dpmixl(i,j,n))
              endif !iv
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
! ---     rhs: dpmixl.n, uflux+, vflux+
! ---     lhs: dpmixl.n
!
          margin = mbdy - 6
!
!$OMP     PARALLEL DO PRIVATE(j,i) &
!$OMP              SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              if (SEA_P) then
                dpmixl(i,j,n)=dpmixl(i,j,n)- &
                               ((uflux(i+1,j)-uflux(i,j))+ &
                                (vflux(i,j+1)-vflux(i,j)))*scp2i(i,j)
              endif !ip
            enddo !i
          enddo !j
!$OMP     END PARALLEL DO
!
        endif  ! end: (thkdf2.ne.0.)
!
      endif  ! end: (hybrid .and. mxlkta)
!
! --- cumalative fluxes
!
! --- rhs: uflxav, vflxav, dpav, uflx, vflx, dp
! --- lhs: uflxav, vflxav, dpav
!
      margin = 0
!
!$OMP PARALLEL DO PRIVATE(j,k,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do k=1,kk
          do i=1-margin,ii+margin
            if (SEA_U) then
              uflxav(i,j,k)=uflxav(i,j,k)+uflx(i,j,k)
            endif !iu
            if (SEA_V) then
              vflxav(i,j,k)=vflxav(i,j,k)+vflx(i,j,k)
            endif !iv
            if (SEA_P) then
              dpav(i,j,k)=dpav(i,j,k)+dp(i,j,k,n)
            endif !ip
          enddo !i
        enddo !k
      enddo !j
!$OMP END PARALLEL DO
!
      if     (btrmas) then
!
! ---   add bottom-pressure restoring term arising from split-explicit treatment
! ---   of continuity equation (step 4 in appendix b to 1992 brhs paper)
! ---   conversion from dp back to dp'
!
! ---   rhs: dp.n, p, pbot
! ---   lhs: dp.n, p, dpmixl.n
!
        margin = 1 !!Alex (2.2.99)
!
!$OMP   PARALLEL DO PRIVATE(j,l,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              onetacnt(i,j) = p(i,j,kk+1)/pbot(i,j)   !used in barotp
              if (onetacnt(i,j) .gt. epsil_cnuity) then
                q=1.0/onetacnt(i,j)
                k=1
                  dp(i,j,k,n)=dp(i,j,k,n)*q
                   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                do k=2,kk
                  dp(i,j,k,n)=dp(i,j,k,n)*q
                   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                enddo !k
              else !dry point?
                k=1
                  dp(i,j,1,n)=pbot(i,j)
                   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                do k=2,kk
                  dp(i,j,k,n)=0.0
                   p(i,j,k+1)=p(i,j,k)+dp(i,j,k,n)
                enddo !k
              endif !onetacnt
              if (isopyc) then
                dpmixl(i,j,n)=dp(i,j,1,n)
              endif
            endif ! ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
      endif !btrmas
!
! --- Robert-Asselin time filter of thickness field
!
      mbdy = 6
!
! --- dp.m and dpo.n halos are already up to date
      call xctilr(dp(1-nbdy,1-nbdy,1,n),1,kk, 6,6, halo_ps)
!
      margin = mbdy
!
!$OMP PARALLEL DO PRIVATE(j,i,k, &
!$OMP                     dpold,dpmid,dpnew,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            do k=1,kk
              dpold = dpo(i,j,k,n)  !t-1
              dpmid = dp( i,j,k,m)  !t
              dpnew = dp( i,j,k,n)  !t+1
              q     = 0.5*ra2fac*(dpold+dpnew-2.0*dpmid)
              dpo(i,j,k,m) = dp(i,j,k,m)  !t
              dp( i,j,k,m) = dpmid + q    !t & RA
            enddo !k
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      return
      end subroutine cnuity
!
!
!> Revision history:
!>
!> July 1997 - combined diff. and antidiff.flux calc. (eliminated loops 20,21)
!> Aug. 1997 - set u/vflux=0 before entering thickness smoothing k-loop 13
!> Jul. 1998 - reversed i/j loop nesting in loop 26
!> Nov. 1999 - added code for vertical advection of the mixed layer base for the
!>             kraus-turner mixed layer model (117 and 118 loops)
!> Apr. 2000 - changed i/j loop nesting to j/i
!> May  2000 - added code to eliminate neg. dp resulting from intfc.smoothing
!> Aug. 2000 - loop 119 executed only when hybrid vertical coordinate is used
!> Dec. 2000 - added biharmonic diffusion of KTa mixed layer
!> May  2002 - thickness diffusion coefficent based on max(sc?x,sc?y)
!> Oct  2003 - allow spacially varying thkdf4
!> Apr  2004 - check for neg. dp every 3 timesteps, fatal if < dpfatal
!> Aug  2011 - replaced dpold,dpoldm with dpo
!> Aug  2011 - apply Robert-Asselin filter to dp here (used to be in tsadvc)
!> Jan  2012 - check for zero thickness top layer
!> Aug  2012 - fixed model hang on neg. dp bug
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> May  2014 - use ipim1, ipip1, ipjm1, ipjp1  as sea-only neighbors
!> May  2014 - use ipim1x,ipip1x,ipjm1x,ipjp1x as sea-only neighbors
!> Jan  2018 - spacially varying thkdf2 is now allowed
!> Aug. 2018 - btrmas added, use onetamas to simplify logic
!> Nov. 2018 - added oneta_u and oneta_v to correct and simplify logic
