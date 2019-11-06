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
      module mod_tsadvc
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none
!
! --- module for tsadvc and related routines
!
      private !! default is private
      public  :: tsadvc
!
      integer, parameter, dimension (0:4) :: &
        mbdy_advtyp = (/ 2,     & !PCM
                         5,     & !MPDATA
                         5,     & !FCT2
                         0,     & !N/A
                         5 /)     !FCT4
!
      logical, parameter :: lpipe_advem=.false.  !extra checking (when pipe on)
      logical, parameter :: lconserve  =.false.  !explicitly conserve the field
!
      logical, save ::  ldebug_tsdif   !switch to debug diffusion - usually .false.
      logical, save ::  ldebug_advem   !switch to debug advection - usually .false.
      integer, save ::  itests,jtests  !local copy of itest,jtest
!
#if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
#else
      real, save, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
#endif
         fmx,fmn       & ! local max,min
        ,flx,fly       & ! fluxes
        ,fldlo         & ! lo order solution
        ,fmxlo,fmnlo   & ! local min
        ,fax,fay       & ! fluxes
        ,rp,rm         & ! FCT/MPDATA terms
        ,flxdiv        & ! flux divergence
        ,tx1,ty1       & ! MPDATA terms
        ,fldao,fldan     ! total field quantity (old/center, new)

#if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
#else
      real,    save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: &
#endif
          uloc,vloc,hloc,dtloc,ucumdt,vcumdt,flxcum,flycum

#if defined(RELO)
      logical, save, allocatable, dimension(:,:) :: &
#else
      logical, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy):: &
#endif
          lcalc

      contains

      subroutine advem(advtyp,fld,fldc,u,v,fco,fcn,posdef, &
                       scal,scali,dt2,btrmas)
      implicit none
!
      logical, intent(in)    :: btrmas
      integer, intent(in)    :: advtyp
      real,    intent(in)    :: posdef,dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: fldc,u,v,fco,fcn,scal,scali
!
! --- wrapper for advection schemes
!
! --- a recent text on advection schemes is:
! --- D.R. Durran (1999): Numerical Methods for wave equations in
! --- geophysical fluid dynamics, Springer.
!
! --- advtyp= 0 for 1st order PCM    (Donor Cell)
! --- advtyp= 1 for 2nd order MPDATA (old to new, as in 2.1.03)
! --- advtyp= 2 for 2nd order FCT    (Leapfrog time step)
! --- advtyp= 4 for 4th order FCT    (Leapfrog time step)
!
! --- time steps are "old", "center" and "new".
!
! --- fld    - scalar field, at old time step on input but new on output
! --- fldc   - scalar field, at center time step
! --- u,v    - mass fluxes satisfying continuity equation (old to new)
! --- fco    - thickness of the layer at old time step
! --- fcn    - thickness of the layer at new time step
! --- posdef - offset for MPDATA to make the field positive
! --- scal   - spatial increments (squared)
! --- scali  - inverse of scal
! --- dt2    - temporal increment (from old to new, i.e. two time steps)
!
!  on return, fld's valid halo will be 0 wide.
!
      real    offset
      real*8  sumold,sumnew,sumcor
      integer i,j
!
#if defined(RELO)
      if     (.not.allocated(fmx)) then
        allocate( &
                 fmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 fmn(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 flx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 fly(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               fldlo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               fmxlo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               fmnlo(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 fax(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 fay(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                  rp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                  rm(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              flxdiv(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 tx1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 ty1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               fldao(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               fldan(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 16*(idm+2*nbdy)*(jdm+2*nbdy) )
                 fmx = r_init
                 fmn = r_init
                 flx = r_init
                 fly = r_init
               fldlo = r_init
               fmxlo = r_init
               fmnlo = r_init
                 fax = r_init
                 fay = r_init
                  rp = r_init
                  rm = r_init
              flxdiv = r_init
                 tx1 = r_init
                 ty1 = r_init
               fldao = r_init
               fldan = r_init
      endif !.not.allocated
#endif
!
      if     (advtyp.eq.0) then
        call advem_pcm(   fld,     u,v,fco,fcn,       scal,scali,dt2)
      elseif (advtyp.eq.1) then
        call advem_mpdata(fld,     u,v,fco,fcn,posdef,scal,scali,dt2)
      elseif (advtyp.eq.2 .and. btrmas) then
        call advem_fct2c( fld,fldc,u,v,fco,fcn,       scal,scali,dt2)
      elseif (advtyp.eq.2) then !.not.btrmas
        call advem_fct2(  fld,fldc,u,v,fco,fcn,       scal,scali,dt2)
      elseif (advtyp.eq.4) then
        call advem_fct4(  fld,fldc,u,v,fco,fcn,       scal,scali,dt2)
      else
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i4 /)') &
            'error: advem called with advtyp =',advtyp
        endif
        call xcstop('advem')
               stop 'advem'
      endif
!
      if     (lconserve) then !usually .false.
!
! ---   explicit conservation of tracer (should not be needed).
!
        call xcsum(sumold, fldao,ipa)
        call xcsum(sumnew, fldan,ipa)
!
        if     (sumnew.ne.0.0) then
          offset = (sumold-sumnew)/sumnew
        else
          offset = 0.0
        endif
!
!$OMP   PARALLEL DO PRIVATE(j,i,offset) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1,jj
          do i=1,ii
            if (SEA_P) then
              fld(i,j)=fld(i,j)*(1.0+offset)
!
!diag         fldan(i,j) = fld(i,j)*fcn(i,j)*scal(i,j)
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
        if     (lpipe .and. lpipe_advem) then
! ---     compare two model runs.
          call pipe_compare_sym1(fld,    ip,'ad:oset:fld ')
        endif
!
!diag   call xcsum(sumcor, fldan,ipa)
!diag   if     (mnproc.eq.1) then
!diag     write(lp,'(a,1p4e16.8)') &
!diag       'advem: ',sumold,sumnew,sumcor,offset
!diag   endif
      endif !lconserve
      return
      end subroutine advem

      subroutine advem_mpdata(fld,u,v,fco,fcn,posdef,scal,scali,dt2)
      implicit none
!
      real,    intent(in)    :: posdef,dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: u,v,scal,scali,fco,fcn
!
! LeapFrog 2nd order MPDATA.
! combined monotone scheme, for details see section 3.3 (eqs. 34 to 37)
! in smolarkiewicz and clark, 1986, j.comput.phys.,67,no 2, p. 396-438
! and smolarkiewicz and grabowski, 1989, j.comput.phys. and recently
! P.K. Smolarkiewicz and L.J. Margolin (1998): MPDATA: A finite
! difference solver for geophysical flows, J.Comput.Phys. 140 459-480.
!
! time steps are "old", "center" and "new".
!
!  fld    - scalar field, must be >0, old input but new output
!  u,v    - mass fluxes satisfying continuity equation (old to new)
!  fco    - thickness of the layer at old time step
!  fcn    - thickness of the layer at new time step
!  posdef - offset to make the field positive
!  scal   - spatial increments (squared)
!  scali  - inverse of scal
!  dt2    - temporal increment (from old to new)
!
!  on return, fld's valid halo will be 0 wide.
!
      real,    parameter :: onemu=9806.e-12 !very small layer thickness
!
      real    fcn2,fco2,flxdn,flxdp,flydn,flydp,q
      integer i,j,l,ia,ib,ja,jb,mbdy_a,margin
!
      mbdy_a = mbdy_advtyp(1)  ! = 5
!
! --- compute low-order and part of antidiffusive fluxes
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly, fmx, fmn
!
      margin = mbdy_a - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,q,jb,ja,ib,ia) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            tx1(i,j)=.5*abs(u(i,j))*(fld(i,j)-fld(i-1,j))
            if (u(i,j).ge.0.0) then
              q=fld(i-1,j)
            else
              q=fld(i  ,j)
            endif
            flx(i,j)=u(i,j)*(q+posdef)
          endif !iu
          if (SEA_V) then
            ty1(i,j)=.5*abs(v(i,j))*(fld(i,j)-fld(i,j-1))
            if (v(i,j).ge.0.0) then
              q=fld(i,j-1)
            else
              q=fld(i,j  )
            endif
            fly(i,j)=v(i,j)*(q+posdef)
          endif !iv
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))+posdef
            fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))+posdef
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(tx1, iu,'ad:11:tx1   ', &
                               ty1, iv,'ad:11:ty1   ')
        call pipe_compare_sym2(flx, iu,'ad:11:flx   ', &
                               fly, iv,'ad:11:fly   ')
        call pipe_compare_sym1(fmx, ip,'ad:11:fmx   ')
        call pipe_compare_sym1(fmn, ip,'ad:11:fmn   ')
      endif
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly, fmx, fmn
!
      margin = mbdy_a - 1
!
      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            flx(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            flx(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fly(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fly(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:22:flx   ', &
                               fly, iv,'ad:33:fly   ')
      endif
!diag if     (itests.gt.0 .and. jtests.gt.0) then
!diag   i=itests
!diag   j=jtests
!diag   write (lp,'(a,2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag              1pe9.2,0pf9.3/1pe39.2/0pf39.3)') &
!diag     'advem (1)',i+i0,j+j0, &
!diag     fld(i-1,j),u(i,j),fld(i,j-1),v(i,j), &
!diag     fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
!diag endif
!
! --- rhs: flx+, fly+, fco, fmn, fmx, fcn
! --- lhs: fld
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
!...........lo order Donor Cell step
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+ &
                         (fly(i,j+1)-fly(i,j)) )*dt2*scali(i,j)
            q=(fld(i,j)+posdef)*(fco(i,j)+onemu)-flxdiv(i,j)
            !max,min should only be active for very thin layers
            fldlo(i,j)=max( fmn(i,j), min( fmx(i,j), &
                                           q/(fcn(i,j)+onemu) ))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(fldlo,  ip,'ad:610:fldlo')
        call pipe_compare_sym1(flxdiv, ip,'ad:610:flxdv')
      endif
!
! --- finish computation of antidiffusive fluxes
!
! --- rhs: tx1, u, ty1, v, flxdiv+, fco+, fcn+
! --- lhs: flx, fly
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,fcn2,fco2) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        fco2=0.0
        fcn2=0.0
        do i=1-margin,ii+margin
          if (SEA_U) then
            fco2=fco(i,j)+fco(i-1,j)  ! inforce order on flx calc
            fcn2=fcn(i,j)+fcn(i-1,j)  ! inforce order on flx calc
            flx(i,j)=tx1(i,j)-u(i,j)*(flxdiv(i,j)+flxdiv(i-1,j)) &
               /((fco2+fcn2)+onemu)
          endif !iu
          if (SEA_V) then
            fco2=fco(i,j)+fco(i,j-1)  ! inforce order on fly calc
            fcn2=fcn(i,j)+fcn(i,j-1)  ! inforce order on fly calc
            fly(i,j)=ty1(i,j)-v(i,j)*(flxdiv(i,j)+flxdiv(i,j-1)) &
               /((fco2+fcn2)+onemu)
          endif !iv
        enddo !i
        if (fco2*fcn2.eq.1.e30) flx(1-nbdy,j)=0.0  ! prevent removal of fc*2
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx,  iu,'ad: 8:flx   ', &
                               fly,  iv,'ad: 8:fly   ')
      endif
!
! --- limit antidiffusive fluxes
! --- rp and rm used to be called flp and fln
!
! --- rhs: fmx, fmn, fldlo, fcn, flx+, fly+
! --- lhs: rp, rm
!
      margin = mbdy_a - 3
!
!$OMP PARALLEL DO PRIVATE(j,i,flxdn,flxdp,flydn,flydp) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            flxdp=min(0.0,flx(i+1,j))-max(0.0,flx(i,j))
            flxdn=max(0.0,flx(i+1,j))-min(0.0,flx(i,j))
            flydp=min(0.0,fly(i,j+1))-max(0.0,fly(i,j))
            flydn=max(0.0,fly(i,j+1))-min(0.0,fly(i,j))
            rp(i,j)=(fmx(i,j)-fldlo(i,j))*(fcn(i,j)*scal(i,j))/ &
             ((onemu-(flxdp+flydp))*dt2)
            rm(i,j)=(fldlo(i,j)-fmn(i,j))*(fcn(i,j)*scal(i,j))/ &
             ((onemu+(flxdn+flydn))*dt2)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(rp, ip,'ad:16:flp   ')
        call pipe_compare_sym1(rm, ip,'ad:16:fln   ')
      endif
!
! --- rhs: flx, fly, rp+, rm+
! --- lhs: flx, fly
!
      margin = mbdy_a - 4
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            flx(i,j)=max(0.0,flx(i,j))*min(1.0,rp(i,j),rm(i-1,j)) &
                    +min(0.0,flx(i,j))*min(1.0,rp(i-1,j),rm(i,j))
          endif !iu
          if (SEA_V) then
            fly(i,j)=max(0.0,fly(i,j))*min(1.0,rp(i,j),rm(i,j-1)) &
                    +min(0.0,fly(i,j))*min(1.0,rp(i,j-1),rm(i,j))
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:18:flx   ', &
                               fly, iv,'ad:18:fly   ')
      endif
!
!diag i=itests
!diag j=jtests
!diag write (lp,'(''advem (2)'',2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag 1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fldlo(i-1,j),u(i,j),fldlo(i,j-1), &
!diag v(i,j),fldlo(i,j),v(i,j+1),fldlo(i,j+1),u(i+1,j),fldlo(i+1,j)
!
! --- rhs: flx+, fly+, fldlo, fcn, fmx
! --- lhs: flxdiv, fld
!
      margin = mbdy_a - 5
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            fldao(i,j) = fld(i,j)*fco(i,j)*scal(i,j)
!
!...........apply antidiffusive flux correction
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+ &
                         (fly(i,j+1)-fly(i,j)) )*dt2*scali(i,j)
            !max,min should only be active for very thin layers
            fld(i,j)=max( fmn(i,j), min( fmx(i,j), &
                          fldlo(i,j)-flxdiv(i,j)/(fcn(i,j)+onemu) ))
            fld(i,j)=fld(i,j)-posdef
!
            fldan(i,j) = fld(i,j)*fcn(i,j)*scal(i,j)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(flxdiv, ip,'ad:620:flxdv')
        call pipe_compare_sym1(fld,    ip,'ad:1620:fld ')
      endif
      return
      end subroutine advem_mpdata

      subroutine advem_pcm(fld,u,v,fco,fcn,scal,scali,dt2)
      implicit none
!
      real,    intent(in)    :: dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: u,v,scal,scali,fco,fcn
!
! Piecewise Constant Method (Donor Cell, Upwind)
! Over two time steps (may require half the normal time step for stability).
!
! time steps are "old", "center" and "new".
!
!  fld   - scalar field, need not be >0, old input but new output
!  u,v   - mass fluxes satisfying continuity equation (old to new)
!  fco   - thickness of the layer at old    time step
!  fcn   - thickness of the layer at new    time step
!  scal  - spatial increments (squared)
!  scali - inverse of scal
!  dt2   - temporal increment (from old to new)
!
!  on return, fld's valid halo will be 0 wide.
!
      real,    parameter :: onemu=9806.e-12 !very small layer thickness
!
      real    q
      integer i,j,l,ia,ib,ja,jb,mbdy_a,margin
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly
!
      mbdy_a = mbdy_advtyp(0)  ! = 2
!
      margin = mbdy_a - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (u(i,j).ge.0.0) then
              q=fld(i-1,j)
            else
              q=fld(i  ,j)
            endif
            flx(i,j)=u(i,j)*q
          endif !iu
          if (SEA_V) then
            if (v(i,j).ge.0.0) then
              q=fld(i,j-1)
            else
              q=fld(i,j  )
            endif
            fly(i,j)=v(i,j)*q
          endif !iv
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
            fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:11:flx   ', &
                               fly, iv,'ad:11:fly   ')
      endif
!
      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            flx(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            flx(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fly(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fly(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:22:flx   ', &
                               fly, iv,'ad:33:fly   ')
      endif
!diag if     (itests.gt.0 .and. jtests.gt.0) then
!diag   i=itests
!diag   j=jtests
!diag   write (lp,'(a,2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag              1pe9.2,0pf9.3/1pe39.2/0pf39.3)') &
!diag     'advem (1)',i+i0,j+j0, &
!diag     fld(i-1,j),u(i,j),fld(i,j-1),v(i,j), &
!diag     fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
!diag endif
!
! --- rhs: flx+, fly+, fld, fco, fcn
! --- lhs: flxdiv, fld
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            fldao(i,j) = fld(i,j)*fco(i,j)*scal(i,j)
!
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+ &
                         (fly(i,j+1)-fly(i,j)) )*dt2*scali(i,j)
            q=fld(i,j)*(fco(i,j)+onemu)-flxdiv(i,j)
            !max,min should only be active for very thin layers
            fld(i,j)=max( fmn(i,j), min( fmx(i,j), &
                                         q/(fcn(i,j)+onemu) ))
!
            fldan(i,j) = fld( i,j)*fcn(i,j)*scal(i,j)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(fld,    ip,'ad:610:fld  ')
        call pipe_compare_sym1(flxdiv, ip,'ad:610:flxdv')
      endif
!
!diag i=itests
!diag j=jtests
!diag write (lp,'(''advem (2)'',2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag 1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fld(i-1,j),u(i,j),fld(i,j-1), &
!diag v(i,j),fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
      return
      end subroutine advem_pcm

      subroutine advem_fct2(fld,fldc,u,v,fco,fcn,scal,scali,dt2)
      implicit none
!
      real,    intent(in)    :: dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: fldc,u,v,scal,scali,fco,fcn
!
! Leapfrog 2nd order FCT
! S.T. Zalesak (1979): Fully multidimensional flux-corrected
! transport algorithms for fluids, J.Comput.Phys. 31 335-362.
!
! time steps are "old", "center" and "new".
!
!  fld   - scalar field, need not be >0, old input but new output
!  fldc  - scalar field at center time step
!  u,v   - mass fluxes satisfying continuity equation (old to new)
!  fco   - thickness of the layer at old    time step
!  fcn   - thickness of the layer at new    time step
!  scal  - spatial increments (squared)
!  scali - inverse of scal
!  dt2   - temporal increment (from old to new)
!
!  on return, fld's valid halo will be 0 wide.
!
      real,    parameter :: onemu=9806.e-12 !very small layer thickness
!
      real    flxdn,flxdp,flydn,flydp,q
      integer i,j,l,ia,ib,ja,jb,mbdy_a,margin
!
      real :: fhx, fhy, fqmax, fqmin, famax, famin
      real :: qdt2, qp, qm, fact

      mbdy_a = mbdy_advtyp(2)  ! = 5
!
! --- compute low-order and part of antidiffusive fluxes
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly, fmx, fmn
!
      margin = mbdy_a - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,q,jb,ja,ib,ia) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (u(i,j).ge.0.0) then
              q=fld(i-1,j)
            else
              q=fld(i  ,j)
            endif
            flx(i,j)=u(i,j)*q
          endif !iu
          if (SEA_V) then
            if (v(i,j).ge.0.0) then
              q=fld(i,j-1)
            else
              q=fld(i,j  )
            endif
            fly(i,j)=v(i,j)*q
          endif !iv
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
            fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
!       call pipe_compare_sym2(tx1, iu,'ad:11:tx1   ',
!    &                         ty1, iv,'ad:11:ty1   ')
!       call pipe_compare_sym2(flx, iu,'ad:11:flx   ',
!    &                         fly, iv,'ad:11:fly   ')
      endif
!
      if     (ldebug_advem .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        write (lp,'(a,2i5,3f10.5)') &
          'advem: fld, rng  ',i+i0,j+j0, &
          fld(i,j),fmx(i,j),fmn(i,j)
      endif
!
      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            flx(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            flx(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fly(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fly(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:22:flx   ', &
                               fly, iv,'ad:33:fly   ')
      endif
!
!diag if     (ldebug_advem .and. itests.gt.0 .and. jtests.gt.0) then
!diag   i=itests
!diag   j=jtests
!diag   write (lp,'(a,2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag              1pe9.2,0pf9.3/1pe39.2/0pf39.3)') &
!diag     'advem (1)',i+i0,j+j0, &
!diag     fld(i-1,j),u(i,j),fld(i,j-1),v(i,j), &
!diag     fld(i,j),v(i,j+1),fld(i,j+1),u(i+1,j),fld(i+1,j)
!diag endif
!
! --- rhs: flx+, fly+, fld, fmx, fmn, fco, fcn
! --- lhs: fldlo, fmxlo, fmnlo, flxdiv
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
!...........lo order Donor Cell step
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+ &
                         (fly(i,j+1)-fly(i,j)) )*dt2*scali(i,j)
            q=fld(i,j)*(fco(i,j)+onemu)-flxdiv(i,j)
            !max,min should only be active for very thin layers
            fldlo(i,j)=max( fmn(i,j), min( fmx(i,j), &
                                           q/(fcn(i,j)+onemu) ))
            fmxlo(i,j) = max(fld(i,j),fldc(i,j),fldlo(i,j))
            fmnlo(i,j) = min(fld(i,j),fldc(i,j),fldlo(i,j))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(fldlo,  ip,'ad:610:fldlo')
        call pipe_compare_sym1(flxdiv, ip,'ad:610:flxdv')
      endif
!
      if     (ldebug_advem .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        write (lp,'(a,2i5,3f10.5)') &
          'advem: fldlo,rng ',i+i0,j+j0, &
          fldlo(i,j),fmnlo(i,j),fmxlo(i,j)
      endif
!
!.....Leapfrog step using high order scheme
!
! --- rhs: u, v, fld+, flx, fly
! --- lhs: fax, fay
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q,jb,ja,ib,ia,fhx,fhy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            fhx=u(i,j)*0.5*(fldc(i,j)+fldc(i-1,j))  ! 2nd order in space
            fax(i,j)= fhx-flx(i,j)                  ! anti-diffusion x-flux
          endif !iu
          if (SEA_V) then
            fhy=v(i,j)*0.5*(fldc(i,j)+fldc(i,j-1))  ! 2nd order in space
            fay(i,j)= fhy-fly(i,j)                  ! anti-diffusion y-flux
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO

      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            fax(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            fax(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fay(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fay(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
!========================================================
!
! --- finish computation of antidiffusive fluxes
!
! --- rhs: fmnlo+,fmxlo+,fax+,fay+,fldlo,fcn,scal
! --- lhs: rp, rm
!
      margin = mbdy_a - 3
!
      qdt2 = 1.0/dt2
!$OMP PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb, &
!$OMP                     fqmax,fqmin,famax,famin,qp,qm) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fqmax = max(fmxlo(i,j),fmxlo(ia,j),fmxlo(ib,j), &
                                   fmxlo(i,ja),fmxlo(i,jb))
            fqmin = min(fmnlo(i,j),fmnlo(ia,j),fmnlo(ib,j), &
                                   fmnlo(i,ja),fmnlo(i,jb))
            famax = max(0.0,fax(i, j) ) - min(0.0,fax(ib,j) ) + &
                    max(0.0,fay(i, j) ) - min(0.0,fay(i, jb))
            famin = max(0.0,fax(ib,j) ) - min(0.0,fax(i, j) ) + &
                    max(0.0,fay(i, jb)) - min(0.0,fay(i, j) )
            if (famax > 0.0) then
              qp = (fqmax-fldlo(i,j)) *fcn(i,j)*scal(i,j)*qdt2
              if (qp > famax) then
                 rp(i,j) =1.0
              else
                  rp(i,j) = min(1.0, qp/famax)
              endif
            else
              rp(i,j) = 0.0
            endif
            if (famin > 0.0) then
              qm = (fldlo(i,j)-fqmin) *fcn(i,j)*scal(i,j)*qdt2
              if (qm>famin) then
                 rm(i,j) = 1.0
              else  
                 rm(i,j) = min(1.0, qm/famin)
              endif
            else
              rm(i,j) = 0.0
            endif
            fmx(i,j) = fqmax  !less restrictive maximum
            fmn(i,j) = fqmin  !less restrictive minimum
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (ldebug_advem .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        write (lp,'(a,2i5,3f10.5)') &
          'advem: fldc,rngfq',i+i0,j+j0, &
          fldc(i,j),fmn(i,j),fmx(i,j)
      endif
!
! --- rhs: rp+, rm+
! --- lhs: fax, fay
!
      margin = mbdy_a - 4
!
!$OMP PARALLEL DO PRIVATE(j,i,fact) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (fax(i,j) < 0.0) then
              fact = min(rp(i-1,j),rm(i,j))
            else
              fact = min(rp(i,j),rm(i-1,j))
            endif
            fax(i,j) = fact*fax(i,j)
          endif !iu
          if (SEA_V) then
            if (fay(i,j) < 0.0) then
              fact = min(rp(i,j-1),rm(i,j))
            else
              fact = min(rp(i,j),rm(i,j-1))
            endif
            fay(i,j) = fact*fay(i,j)
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:18:flx   ', &
                               fly, iv,'ad:18:fly   ')
      endif
!
!diag i=itests
!diag j=jtests
!diag write (lp,'(''advem (2)'',2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag 1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fldlo(i-1,j),u(i,j),fldlo(i,j-1), &
!diag v(i,j),fldlo(i,j),v(i,j+1),fldlo(i,j+1),u(i+1,j),fldlo(i+1,j)
!
! --- rhs: fax+, fay+, fld, fcn, fmx, fmn
! --- lhs: fld
!
      margin = mbdy_a - 5
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            fldao(i,j) = fld(i,j)*fco(i,j)*scal(i,j)
!
!...........apply antidiffusive flux correction
            flxdiv(i,j)=((fax(i+1,j)-fax(i,j))+ &
                         (fay(i,j+1)-fay(i,j)) )*dt2*scali(i,j)
            !max,min should only be active for very thin layers
            fld(i,j)=max( fmn(i,j), min( fmx(i,j), &
                          fldlo(i,j)-flxdiv(i,j)/(fcn(i,j)+onemu) ))
!
            fldan(i,j) = fld(i,j)*fcn(i,j)*scal(i,j)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(flxdiv, ip,'ad:620:flxdv')
        call pipe_compare_sym1(fld,    ip,'ad:1620:fld ')
      endif
!
      if     (ldebug_advem .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        write (lp,'(a,2i5,f10.5)') &
          'advem: fld       ',i+i0,j+j0, &
          fld(i,j)
      endif
      return
      end subroutine advem_fct2

      subroutine advem_fct2c(fld,fldc,u,v,fco,fcn, scal,scali,dt2)
      implicit none
!
      real,    intent(in)    :: dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: fldc,u,v,scal,scali,fco,fcn
!
! Leapfrog 2nd order FCT
! S.T. Zalesak (1979): Fully multidimensional flux-corrected
! transport algorithms for fluids, J.Comput.Phys. 31 335-362.
!
! time steps are "old", "center" and "new".
!
!  fld   - scalar field, need not be >0, old input but new output
!  fldc  - scalar field at center time step
!  u,v   - mass fluxes satisfying continuity equation (old to new)
!  fco   - thickness of the layer at old    time step
!  fcn   - thickness of the layer at new    time step
!  scal  - spatial increments (squared)
!  scali - inverse of scal
!  dt2   - temporal increment (from old to new)
!
!  on return, fld's valid halo will be 0 wide.
!
!  Remy Baraille, SHOM, September 2008.
!
      real,    parameter :: epsil=1.e-10
!
      real    flxdn,flxdp,flydn,flydp,q
      integer i,j,l,ia,ib,ja,jb,mbdy_a
      integer iter,margin
!
      real :: fhx, fhy, fqmax, fqmin, famax, famin
      real :: qdt2, qp, qm, fact
!
#if defined(RELO)
      if     (.not.allocated(uloc)) then
       allocate( &
                 uloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 vloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 hloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                dtloc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               ucumdt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               vcumdt(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               flxcum(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               flycum(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                lcalc(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))
        call mem_stat_add( 9*(idm+2*nbdy)*(jdm+2*nbdy) )
                 uloc = r_init
                 vloc = r_init
                 hloc = r_init
                dtloc = r_init
               ucumdt = r_init
               vcumdt = r_init
               flxcum = r_init
               flycum = r_init
               lcalc  = .false.
      endif !.not.allocated
#endif
!
      mbdy_a = mbdy_advtyp(2)  ! = 5
!
! --- compute low-order and part of antidiffusive fluxes
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly, fmx, fmn
!
       lcalc(:,:) = .true.
      flxcum(:,:) = 0.0
      flycum(:,:) = 0.0
       dtloc(:,:) = 0.0
        hloc(:,:) = fco(:,:)
       fldlo(:,:) = fld(:,:)
!
! --- explicitly set to zero so that arrays are zero over land
      ucumdt(:,:) = 0.0
      vcumdt(:,:) = 0.0
        uloc(:,:) = 0.0
        vloc(:,:) = 0.0
         flx(:,:) = 0.0
         fly(:,:) = 0.0
!
      do iter=1,5
!
        margin = mbdy_a
!
!$OMP   PARALLEL DO PRIVATE(j,i,q) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin, jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              q = max(u(i+1,j),0.0)-min(u(i,j),0.0) &
                + max(v(i,j+1),0.0)-min(v(i,j),0.0)
              if (q.gt.0.0) then
                 dtloc(i,j) = min( dt2, hloc(i,j) / (q*scali(i,j)))
              else
                 dtloc(i,j) = dt2
              endif
            endif !ip
          enddo !i
        enddo !j
!
        margin = mbdy_a - 1
!
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              if (ucumdt(i,j).ne.dt2) then
                if (u(i,j).ge.0) then
                  uloc(i,j) = min(dt2-ucumdt(i,j),dtloc(i-1,j))*u(i,j)
                  flx(i,j)  = fldlo(i-1,j)*uloc(i,j)
                  ucumdt(i,j) = ucumdt(i,j) +  &
                                min(dt2-ucumdt(i,j),dtloc(i-1,j))
                else !u(i,j).lt.0
                  uloc(i,j) = min(dt2-ucumdt(i,j),dtloc(i,j))*u(i,j)
                  flx(i,j)  = fldlo(i,j)*uloc(i,j)
                  ucumdt(i,j) = ucumdt(i,j) +  &
                                min(dt2-ucumdt(i,j),dtloc(i,j))
                endif !u
                flxcum(i,j) = flxcum(i,j) + flx(i,j)
              else !ucumdt==dt2
                uloc(i,j) = 0.0
                flx(i,j)  = 0.0
              endif !ucumdt
            endif !iu
            if (SEA_V) then
              if (vcumdt(i,j).ne.dt2) then
                if (v(i,j).ge.0) then
                  vloc(i,j) = min(dt2-vcumdt(i,j),dtloc(i,j-1))*v(i,j)
                  fly(i,j)  = fldlo(i,j-1)*vloc(i,j)
                  vcumdt(i,j) = vcumdt(i,j) +  &
                                min(dt2-vcumdt(i,j),dtloc(i,j-1))
                else !v.lt.0
                  vloc(i,j) = min(dt2-vcumdt(i,j),dtloc(i,j))*v(i,j)
                  fly(i,j)  = fldlo(i,j)*vloc(i,j)
                  vcumdt(i,j) = vcumdt(i,j) +  &
                                min(dt2-vcumdt(i,j),dtloc(i,j))
                endif
                flycum(i,j) = flycum(i,j) + fly(i,j)
              else !vcumdt==dt2
                vloc(i,j) = 0.0
                fly(i,j)  = 0.0
             endif !vcumdt
           endif !iv
         enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        margin = mbdy_a - 2

!$OMP   PARALLEL DO PRIVATE(j,i,qp) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin, jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              if (lcalc(i,j)) then  !lcalc from prevous iteration
                qp = hloc(i,j) - (uloc(i+1,j)-uloc(i,j)+ &
                                  vloc(i,j+1)-vloc(i,j))*scali(i,j)
!
!                                it may happen that the cfl is violated,
!                                                leading to a negative h
!                                                -----------------------
                if (qp.gt.0.0) then
                  fldlo(i,j) = ((epsil+hloc(i,j))*fldlo(i,j) -  &
                            (flx(i+1,j)-flx(i,j)+ &
                             fly(i,j+1)-fly(i,j))*scali(i,j))/(epsil+qp)
                endif !qp
                hloc(i,j)  = qp
                lcalc(i,j) = ucumdt(i+1,j).ne.dt2.or. &
                             ucumdt(i,j)  .ne.dt2.or. &
                             vcumdt(i,j+1).ne.dt2.or. &
                             vcumdt(i,j)  .ne.dt2
              endif !lcalc
            endif !ip
          enddo !i
        enddo !j
!$OMP  END PARALLEL DO
!
        call xctilr( hloc(1-nbdy,1-nbdy),1,1, mbdy_a,mbdy_a, halo_ps)
        call xctilr(fldlo(1-nbdy,1-nbdy),1,1, mbdy_a,mbdy_a, halo_ps)
!
        if     (lpipe .and. lpipe_advem) then
! ---     compare two model runs.
          call pipe_compare_sym1(fldlo,  ip,'aditer:fldlo')
          call pipe_compare_sym1(hloc,   ip,'aditer:hloc ')
          call pipe_compare_sym1(flxcum, iu,'aditer:xcum ')
          call pipe_compare_sym1(flycum, iv,'aditer:ycum ')
        endif
!
      enddo !iter
!
!.....Leapfrog step using high order scheme
!
! --- rhs: u, v, fld+, flx, fly
! --- lhs: fax, fay
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,fhx,fhy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin, jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            fhx=u(i,j)*0.5*(fldc(i,j)+fldc(i-1,j))  ! 2nd order in space
            fax(i,j)= fhx-flxcum(i,j)/dt2           ! anti-diffusion x-flux
          endif !iu
          if (SEA_V) then
            fhy=v(i,j)*0.5*(fldc(i,j)+fldc(i,j-1))  ! 2nd order in space
            fay(i,j)= fhy-flycum(i,j)/dt2           ! anti-diffusion y-flux
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!              
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.           
        call pipe_compare_sym2(fax, iu,'ad:xcum:fax ', &
                               fay, iv,'ad:ycum:fay ')
      endif                                   
!
      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            fax(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            fax(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fay(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fay(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
!              
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.           
        call pipe_compare_sym2(fax, iu,'ad:ip:0:fax ', &
                               fay, iv,'ad:ip:0:fay ')
      endif                                   
!========================================================
!
! --- finish computation of antidiffusive fluxes
!
! --- rhs: fax+,fay+,fldlo,fcn,scal
! --- lhs: rp, rm
!
      margin = mbdy_a - 3
!
      qdt2 = 1.0/dt2
!$OMP PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb, &
!$OMP                     fqmax,fqmin,famax,famin,qp,qm) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin, jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            ia=i-1; if (ip(ia,j).eq.0) ia=i
            ib=i+1; if (ip(ib,j).eq.0) ib=i
            ja=j-1; if (ip(i,ja).eq.0) ja=j
            jb=j+1; if (ip(i,jb).eq.0) jb=j
            fqmax = max(fldlo(i,j),fldlo(ia,j ),fldlo(ib,j ), &
                                   fldlo(i, ja),fldlo(i, jb))
            fqmin = min(fldlo(i,j),fldlo(ia,j ),fldlo(ib,j ), &
                                   fldlo(i, ja),fldlo(i, jb))
            famax = max(0.0,fax(i,  j)  ) - min(0.0,fax(i+1,j)  ) + &
                    max(0.0,fay(i,  j)  ) - min(0.0,fay(i,  j+1))
            famin = max(0.0,fax(i+1,j)  ) - min(0.0,fax(i,  j)  ) + &
                    max(0.0,fay(i,  j+1)) - min(0.0,fay(i,  j)  )
            if (famax > epsil) then
              qp = (fqmax-fldlo(i,j)) *hloc(i,j)*scal(i,j)*qdt2
              rp(i,j) = min(1.0, qp/famax)
            else
              rp(i,j) = 0.0
            endif
            if (famin > epsil) then
              qm = (fldlo(i,j)-fqmin) *hloc(i,j)*scal(i,j)*qdt2
              rm(i,j) = min(1.0, qm/famin)
            else
              rm(i,j) = 0.0
            endif
            fmx(i,j) = fqmax  !less restrictive maximum
            fmn(i,j) = fqmin  !less restrictive minimum
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
! --- rhs: rp+, rm+
! --- lhs: fax, fay
!
      margin = mbdy_a - 4
!
!$OMP PARALLEL DO PRIVATE(j,i,fact) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin, jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (fax(i,j) < 0.0) then
              fact = min(rp(i-1,j),rm(i,j))
            else
              fact = min(rp(i,j),rm(i-1,j))
            endif
            fax(i,j) = fact*fax(i,j)
          endif !iu
          if (SEA_V) then
            if (fay(i,j) < 0.0) then
              fact = min(rp(i,j-1),rm(i,j))
            else
              fact = min(rp(i,j),rm(i,j-1))
            endif
            fay(i,j) = fact*fay(i,j)
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(fax, iu,'ad:fact:fax ', &
                               fay, iv,'ad:fact:fay ')
      endif
!
      margin = mbdy_a - 5
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin, jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
!...........apply antidiffusive flux correction
            flxdiv(i,j)=((fax(i+1,j)-fax(i,j))+ &
                         (fay(i,j+1)-fay(i,j)) )*dt2*scali(i,j)
            if (hloc(i,j).gt.0.) then
               fld(i,j)=((epsil+hloc(i,j))*fldlo(i,j)-flxdiv(i,j))/ &
                         (epsil+hloc(i,j))
            else
               fld(i,j)=fldlo(i,j)
            endif !hloc
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(flxdiv, ip,'ad:t2c:flxdv')
        call pipe_compare_sym1(fld,    ip,'ad:fct2c:fld')
      endif
!
      return
      end subroutine advem_fct2c

      subroutine advem_fct4(fld,fldc,u,v,fco,fcn,scal,scali,dt2)
      implicit none
!
      real,    intent(in)    :: dt2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(inout) :: fld
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
               intent(in)    :: fldc,u,v,scal,scali,fco,fcn
!
! Leapfrog 4th order FCT
! S.T. Zalesak (1979): Fully multidimensional flux-corrected
! transport algorithms for fluids, J.Comput.Phys. 31 335-362.
!
! time steps are "old", "center" and "new".
!
!  fld   - scalar field, need not be >0, old input but new output
!  fldc  - scalar field at center time step
!  u,v   - mass fluxes satisfying continuity equation (old to new)
!  fco   - thickness of the layer at old    time step
!  fcn   - thickness of the layer at new    time step
!  scal  - spatial increments (squared)
!  scali - inverse of scal
!  dt2   - temporal increment (from old to new)
!
!  on return, fld's valid halo will be 0 wide.
!
      real,    parameter :: onemu=9806.e-12 !very small layer thickness
!
      real,    parameter :: ft14= 7.0/12.0,   & !4th centered inner coeff
                            ft24=-1.0/12.0   !4th centered outer coeff
!
      real    flxdn,flxdp,flydn,flydp,q
      integer i,j,l,ia,ib,ja,jb,mbdy_a,margin
!
      real :: fhx, fhy, fqmax, fqmin, famax, famin
      real :: qdt2, qp, qm, fact

      mbdy_a = mbdy_advtyp(4)  ! = 5
!
! --- compute low-order and part of antidiffusive fluxes
!
! --- rhs: u, v, fld+
! --- lhs: flx, fly, fmx, fmn
!
      margin = mbdy_a - 1
!
!$OMP PARALLEL DO PRIVATE(j,i,q,fhx,fhy,jb,ja,ib,ia) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (u(i,j).ge.0.0) then
              q=fld(i-1,j)
            else
              q=fld(i  ,j)
            endif
            flx(i,j)=u(i,j)*q
          endif !iu
          if (SEA_V) then
            if (v(i,j).ge.0.0) then
              q=fld(i,j-1)
            else
              q=fld(i,j  )
            endif
            fly(i,j)=v(i,j)*q
          endif !iv
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fmx(i,j)=max(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
            fmn(i,j)=min(fld(i,j),fld(ia,j),fld(ib,j), &
                                  fld(i,ja),fld(i,jb))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
!       call pipe_compare_sym2(tx1, iu,'ad:11:tx1   ',
!    &                         ty1, iv,'ad:11:ty1   ')
!       call pipe_compare_sym2(flx, iu,'ad:11:flx   ',
!    &                         fly, iv,'ad:11:fly   ')
      endif
!
      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            flx(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            flx(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fly(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fly(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:22:flx   ', &
                               fly, iv,'ad:33:fly   ')
      endif
!diag if     (itests.gt.0 .and. jtests.gt.0) then
!diag   i=itests
!diag   j=jtests
!diag   write (lp,'(a,2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag              1pe9.2,0pf9.3/1pe39.2/0pf39.3)') &
!diag     'advem (1)',i+i0,j+j0, &
!diag     fldc(i-1,j),u(i,j),fldc(i,j-1),v(i,j), &
!diag     fldc(i,j),v(i,j+1),fldc(i,j+1),u(i+1,j),fldc(i+1,j)
!diag endif
!
! --- rhs: flx+, fly+, fld, fmx, fmn, fco, fcn, flxdiv
! --- lhs: fldlo, fmxlo, fmnlo, flxdiv
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
!...........lo order Donor Cell step
            flxdiv(i,j)=((flx(i+1,j)-flx(i,j))+ &
                         (fly(i,j+1)-fly(i,j)) )*dt2*scali(i,j)
            q=fld(i,j)*(fco(i,j)+onemu)-flxdiv(i,j)
            !max,min should only be active for very thin layers
            fldlo(i,j)=max( fmn(i,j), min( fmx(i,j), &
                                           q/(fcn(i,j)+onemu) ))
            fmxlo(i,j) = max(fld(i,j),fldc(i,j),fldlo(i,j))
            fmnlo(i,j) = min(fld(i,j),fldc(i,j),fldlo(i,j))
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(fldlo,  ip,'ad:610:fldlo')
        call pipe_compare_sym1(flxdiv, ip,'ad:610:flxdv')
      endif
!
!.....Leapfrog step using high order scheme
!
! --- rhs: u, v, fldi++, flx, fly
! --- lhs: fax, fay
!
      margin = mbdy_a - 2
!
!$OMP PARALLEL DO PRIVATE(j,i,q,jb,ja,ib,ia,fhx,fhy) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if     (iu(i-1,j).eq.0 .or. iu(i+1,j).eq.0) then
              ! 2nd order time centered
              fhx=u(i,j)*0.5*(fldc(i,j)+fldc(i-1,j)) 
            else
              ! 4th order time centered
              fhx=u(i,j)*(ft14*(fldc(i,  j)+fldc(i-1,j))+ &
                          ft24*(fldc(i+1,j)+fldc(i-2,j)) )
            endif
            fax(i,j)= fhx-flx(i,j)  ! anti-diffusion x-flux
          endif !iu
          if (SEA_V) then
            if     (iv(i,j-1).eq.0 .or. iv(i,j+1).eq.0) then
              ! 2nd order time centered
              fhy=v(i,j)*0.5*(fldc(i,j)+fldc(i,j-1))
            else
              ! 4th order time centered
              fhy=v(i,j)*(ft14*(fldc(i,j)  +fldc(i,j-1))+ &
                          ft24*(fldc(i,j+1)+fldc(i,j-2)) )
            endif
            fay(i,j)= fhy-fly(i,j)  ! anti-diffusion y-flux
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO

      do j=1-margin,jj+margin
        do l=1,isp(j) !ok
          if     (ifp(j,l).ge. 1-margin) then
            fax(ifp(j,l)  ,j)=0.0
          endif
          if     (ilp(j,l).lt.ii+margin) then
            fax(ilp(j,l)+1,j)=0.0
          endif
        enddo !l
      enddo !j
!
      do i=1-margin,ii+margin
        do l=1,jsp(i) !ok
          if     (jfp(i,l).ge. 1-margin) then
            fay(i,jfp(i,l)  )=0.0
          endif
          if     (jlp(i,l).lt.jj+margin) then
            fay(i,jlp(i,l)+1)=0.0
          endif
        enddo !l
      enddo !i
!========================================================
!
! --- finish computation of antidiffusive fluxes
!
! --- rhs: fmnlo+,fmxlo+,fax+,fay+,fldlo,fcn,scal
! --- lhs: rp, rm
!
      margin = mbdy_a - 3
!
      qdt2 = 1.0/dt2
!$OMP PARALLEL DO PRIVATE(j,i,ia,ib,ja,jb, &
!$OMP                     fqmax,fqmin,famax,famin,qp,qm) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            ia=ipim1(i,j) !i-1 if sea and i if land
            ib=ipip1(i,j) !i+1 if sea and i if land
            ja=ipjm1(i,j) !j-1 if sea and j if land
            jb=ipjp1(i,j) !j+1 if sea and j if land
            fqmax = max(fmxlo(i,j),fmxlo(ia,j),fmxlo(ib,j), &
                                   fmxlo(i,ja),fmxlo(i,jb))
            fqmin = min(fmnlo(i,j),fmnlo(ia,j),fmnlo(ib,j), &
                                   fmnlo(i,ja),fmnlo(i,jb))
            famax = max(0.0,fax(i, j) ) - min(0.0,fax(ib,j) ) + &
                    max(0.0,fay(i, j) ) - min(0.0,fay(i, jb))
            famin = max(0.0,fax(ib,j) ) - min(0.0,fax(i, j) ) + &
                    max(0.0,fay(i, jb)) - min(0.0,fay(i, j) )
            if (famax > 0.0) then
              qp = (fqmax-fldlo(i,j)) *fcn(i,j)*scal(i,j)*qdt2
              rp(i,j) = min(1.0, qp/famax)
            else
              rp(i,j) = 0.0
            endif
            if (famin > 0.0) then
              qm = (fldlo(i,j)-fqmin) *fcn(i,j)*scal(i,j)*qdt2
              rm(i,j) = min(1.0, qm/famin)
            else
              rm(i,j) = 0.0
            endif
            fmx(i,j) = fqmax  !less restrictive maximum
            fmn(i,j) = fqmin  !less restrictive minimum
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
! --- rhs: rp+, rm+
! --- lhs: fax, fay
!
      margin = mbdy_a - 4
!
!$OMP PARALLEL DO PRIVATE(j,i,fact) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_U) then
            if (fax(i,j) < 0.0) then
              fact = min(rp(i-1,j),rm(i,j))
            else
              fact = min(rp(i,j),rm(i-1,j))
            endif
            fax(i,j) = fact*fax(i,j)
          endif !iu
          if (SEA_V) then
            if (fay(i,j) < 0.0) then
              fact = min(rp(i,j-1),rm(i,j))
            else
              fact = min(rp(i,j),rm(i,j-1))
            endif
            fay(i,j) = fact*fay(i,j)
          endif !iv
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym2(flx, iu,'ad:18:flx   ', &
                               fly, iv,'ad:18:fly   ')
      endif
!
!diag i=itests
!diag j=jtests
!diag write (lp,'(''advem (2)'',2i5,f22.3/1pe39.2/0pf21.3,1pe9.2,0pf9.3, &
!diag 1pe9.2,0pf9.3/1pe39.2/0pf39.3)') i,j,fldlo(i-1,j),u(i,j),fldlo(i,j-1), &
!diag v(i,j),fldlo(i,j),v(i,j+1),fldlo(i,j+1),u(i+1,j),fldlo(i+1,j)
!
! --- rhs: fax+, fay+, fld, fcn, fmx, fmn
! --- lhs: fld
!
      margin = mbdy_a - 5
!
!$OMP PARALLEL DO PRIVATE(j,i) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
            fldao(i,j) = fld(i,j)*fco(i,j)*scal(i,j)
!
!...........apply antidiffusive flux correction
            flxdiv(i,j)=((fax(i+1,j)-fax(i,j))+ &
                         (fay(i,j+1)-fay(i,j)) )*dt2*scali(i,j)
            !max,min should only be active for very thin layers
            fld(i,j)=max( fmn(i,j), min( fmx(i,j), &
                          fldlo(i,j)-flxdiv(i,j)/(fcn(i,j)+onemu) ))
!
            fldan(i,j) = fld( i,j)*fcn(i,j)*scal(i,j)
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
      if     (lpipe .and. lpipe_advem) then
! ---   compare two model runs.
        call pipe_compare_sym1(flxdiv, ip,'ad:620:flxdv')
        call pipe_compare_sym1(fld,    ip,'ad:1620:fld ')
      endif
      return
      end subroutine advem_fct4

      subroutine tsadvc(m,n)
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n
!
! --- -----------------------------------------------------
! --- thermodynamic variable(s): advection and diffusion.
! --- -----------------------------------------------------
! --- on entry:
! --- saln(:,:,:,n) = time step t-1
! --- saln(:,:,:,m) = time step t
!
! ---  dpo(:,:,:,n) = time step t-1
! ---  dpo(:,:,:,m) = time step t
! ---  dp( :,:,:,m) = time step t   with RA time smoothing
! ---  dp( :,:,:,n) = time step t+1
!
! --- on exit:
! --- saln(:,:,:,m) = time step t   with RA time smoothing
! --- saln(:,:,:,n) = time step t+1
! --- -----------------------------------------------------
!
      logical, parameter :: lpipe_tsadvc=.false.  !extra checking (when pipe on)
!
#if defined(RELO)
      real, save, allocatable, dimension (:,:) :: &
       sold,told,q2old,q2lold
      real, save, allocatable, dimension (:,:,:) :: &
       trold
#else
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
       sold,told,q2old,q2lold
      real, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,mxtrcr) :: &
       trold
#endif
!
      logical latemp,lath3d,ldtemp,ldth3d
      integer i,isave,j,jsave,k,ktr,l,mbdy,mdf,margin
      real sminn,smaxx,flxdiv,th3d_t, &
           factor,q,dpsold,dpsmid,dpsnew, &
           dpold,dpmid,dpnew,qdpmidn,smin,ufa,ufb,vfa,vfb
      real xmin(kdm),xmax(kdm)
      real sminny(jdm),smaxxy(jdm)
!
      character*12 text,textu,textv
!
! --- for mpdata, select posdef:
! ---   1. as a power of 2
! ---   2. so that the ratio of the standard deviation to the mean of
! ---      each field is approximately the same:
! ---        0 for -saln-,   256 for -temp-, 32 for -th3d-,
! ---        0 for -tracer-,   1 for -q2-
      real       pdzero,pdtemp,pdth3d,pdq2,eps_har
      parameter (pdzero=0.0, pdtemp=256.0, pdth3d=32.0, pdq2=1.0, &
                 eps_har=1.0e-20)
!
      real harmonc,harmonz,a,aa,b,bb
# include "stmt_fns.h"
!
! --- harmonic mean - Conservative for lateral friction
! ---               - force a and b to be non-negative for safety
      harmonz(a, b )=2.0*a*b/max((a+b),2.0*eps_har)
      harmonc(aa,bb)=harmonz(max(aa,0.0),max(bb,0.0))
#if defined(RELO)
!
      if     (.not.allocated(sold)) then
        allocate( &
                 sold(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 told(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
        call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
                 sold = r_init
                 told = r_init
        if (mxlmy) then
          allocate( &
                  q2old(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 q2lold(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) )
          call mem_stat_add( 2*(idm+2*nbdy)*(jdm+2*nbdy) )
                  q2old = r_init
                 q2lold = r_init
        endif !mxlmy
        if     (ntracr.gt.0) then
          allocate( &
                 trold(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ntracr) )
          call mem_stat_add( (idm+2*nbdy)*(jdm+2*nbdy)*ntracr )
                 trold = r_init
        endif !ntracr
      endif !.not.allocated
#endif
!
      ldebug_tsdif = .false.
      ldebug_advem = .false.
!
      itests = itest  !local copy
      jtests = jtest  !local copy
!
      if     (btrmas) then
        onetamas(:,:,n) = oneta(:,:,n)  !use dp rather than dp' for diffusion
        onetamas(:,:,m) = oneta(:,:,n)  !use dp rather than dp' for advection, note "n"
      else
        onetamas(:,:,n) = oneta(:,:,n)  !use dp rather than dp' for diffusion
        onetamas(:,:,m) = 1.0           !use                dp' for advection
      endif
!
      uflux(:,:) = 0.0  !sets land values
      vflux(:,:) = 0.0  !sets land values
!
      mbdy = mbdy_advtyp(abs(advtyp))  ! 2-8 depending on advection scheme
!
      if     (nbdy.lt.mbdy) then
        if     (mnproc.eq.1) then
          write(lp,'(/ a,i3,a /)') &
            'error: nbdy (dimensions.h) must be at least', &
            mbdy,' for the advection scheme indicated by advtyp'
        endif
        call xcstop('tsadvc')
               stop 'tsadvc'
      endif
!
      l = mbdy
! --- dp halo is up to date
      call xctilr(saln(1-nbdy,1-nbdy,1,1),1,2*kk, l,l, halo_ps)
      call xctilr(temp(1-nbdy,1-nbdy,1,1),1,2*kk, l,l, halo_ps)
      call xctilr(th3d(1-nbdy,1-nbdy,1,1),1,2*kk, l,l, halo_ps)
      call xctilr(uflx(1-nbdy,1-nbdy,1  ),1,  kk, l,l, halo_uv)
      call xctilr(vflx(1-nbdy,1-nbdy,1  ),1,  kk, l,l, halo_vv)
      do ktr= 1,ntracr
        call xctilr(tracer(1-nbdy,1-nbdy,1,1,ktr),1,2*kk, l,l, halo_ps)
      enddo !ktr
      if (mxlmy) then
        call xctilr(q2( 1-nbdy,1-nbdy,0,1),1,2*kk+4, l,l, halo_ps)
        call xctilr(q2l(1-nbdy,1-nbdy,0,1),1,2*kk+4, l,l, halo_ps)
      endif
!
      do 81 k=1,kk
!
! --- ---------------------------------------------------
! --- advection of thermodynamic variable(s) (and tracer)
! --- ---------------------------------------------------
!
! --- for isopycnic vertical coordinates:
! ---   advect -th3d- and -saln- in the mixed layer (layer 1),
! ---   advect            -saln- only in all other layers
! --- for hybrid vertical coordinates:
! ---   advect -temp- and -saln- in all layers if advflg==0,
! ---   advect -th3d- and -saln- in all layers if advflg==1
!
      latemp =  k.le.nhybrd .and. advflg.eq.0       ! advect temp
      lath3d = (k.le.nhybrd .and. advflg.eq.1) .or. &
               (k.eq.1      .and. isopyc     )      ! advect th3d
!
! --- smooth mixed-layer mass fluxes in lateral direction
      if(isopyc .and. k.eq.1) then
!
! ---   rhs: vflx+, uflx+
! ---   lhs: vflux, uflux
!
        margin = mbdy - 1
!
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_V) then
              if (iv(i-1,j).ne.0) then
                vfa=vflx(i-1,j,1)
              else 
                vfa=vflx(i,  j,1)
              endif
              if (iv(i+1,j).ne.0) then
                vfb=vflx(i+1,j,1)
              else
                vfb=vflx(i,  j,1)
              endif
              vflux(i,j)=.5*vflx(i,j,1)+.25*(vfa+vfb)
            endif !iv
            if (SEA_U) then
              if (iu(i,j-1).ne.0) then
                ufa=uflx(i,j-1,1)
              else
                ufa=uflx(i,j,  1)
              endif
              if (iu(i,j+1).ne.0) then
                ufb=uflx(i,j+1,1)
              else
                ufb=uflx(i,j,  1)
              endif
              uflux(i,j)=.5*uflx(i,j,1)+.25*(ufa+ufb)
            endif !iu
          enddo !i
        enddo !j
      endif
!
! ---   rhs: temp, saln, uflux+, vflux+, dp
! ---   lhs: told, sold, util1, util2, util3, temp, th3d
!
! ---   util1 = fco = thickness of the layer at old    time step
! ---   util2 = fcn = thickness of the layer at new    time step
!
        margin = mbdy - 1  ! util[12] at mbdy-2
!
!$OMP PARALLEL DO PRIVATE(j,i,ktr,flxdiv) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
        do i=1-margin,ii+margin
          if (SEA_P) then
!
! ---       save for time smoothing
            if     (latemp) then
              told(i,j)=temp(i,j,k,n)
            elseif (lath3d) then
              told(i,j)=th3d(i,j,k,n)
            endif
            sold(i,j)=saln(i,j,k,n)
            do ktr= 1,ntracr
              trold(i,j,ktr)=tracer(i,j,k,n,ktr)
            enddo
            if (mxlmy) then
              q2old( i,j)=q2( i,j,k,n)
              q2lold(i,j)=q2l(i,j,k,n)
            endif
!
! --- before calling 'advem', make sure (a) mass fluxes are consistent
! --- with layer thickness change, and (b) all fields are positive-definite
            if(isopyc .and. k.eq.1) then
              flxdiv=((uflux(i+1,j)  -uflux(i,j)  ) &
                     +(vflux(i,j+1)  -vflux(i,j)  ))*delt1*scp2i(i,j)
            else
              flxdiv=((uflx( i+1,j,k)-uflx( i,j,k)) &
                     +(vflx( i,j+1,k)-vflx( i,j,k)))*delt1*scp2i(i,j)
            endif
            util1(i,j)=max(onetamas(i,j,m)*dp(i,j,k,n)+flxdiv,0.0)  !old, note "m"
            util2(i,j)=max(onetamas(i,j,m)*dp(i,j,k,n),       0.0)  !new, note "m"
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      if     (lpipe .and. lpipe_tsadvc) then
! ---   compare two model runs.
        write(text,'(a10,i2)') '49:sold,k=',k
        call pipe_compare_sym1(sold,ip,text)
        write(text,'(a10,i2)') '49:told,k=',k
        call pipe_compare_sym1(told,ip,text)
        write(text,'(a10,i2)') '49:utl1,k=',k
        call pipe_compare_sym1(util1,ip,text)
        write(text,'(a10,i2)') '49:utl2,k=',k
        call pipe_compare_sym1(util2,ip,text)
        write (textu,'(a9,i3)') 'uflx   k=',k
        write (textv,'(a9,i3)') 'vflx   k=',k
        call pipe_compare_sym2(uflx(1-nbdy,1-nbdy,k),  iu,textu, &
                               vflx(1-nbdy,1-nbdy,k),  iv,textv)
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
!
! --- rhs: temp.[mn], th3d.[mn], saln.[mn], uflx, vflx, util[12]
! --- lhs: temp.n, th3d.n, saln.n
!
      if     (latemp) then
        call advem(advtyp,temp( 1-nbdy,1-nbdy,k,n), &
                          temp( 1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdtemp, scp2,scp2i,delt1, btrmas)
        call advem(advtyp,saln( 1-nbdy,1-nbdy,k,n), &
                          saln( 1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdzero, scp2,scp2i,delt1, btrmas)
      elseif (lath3d .and. hybrid) then
        call advem(advtyp,th3d( 1-nbdy,1-nbdy,k,n), &
                          th3d( 1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdth3d, scp2,scp2i,delt1, btrmas)
        call advem(advtyp,saln( 1-nbdy,1-nbdy,k,n), &
                          saln( 1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdzero, scp2,scp2i,delt1, btrmas)
      elseif (lath3d .and. isopyc) then  ! MICOM-like upper layer
        call advem(advtyp,th3d( 1-nbdy,1-nbdy,k,n), &
                          th3d( 1-nbdy,1-nbdy,k,m), &
                          uflux, &
                          vflux, &
                          util1,util2, &
                          pdth3d, scp2,scp2i,delt1, btrmas)
        call advem(advtyp,saln( 1-nbdy,1-nbdy,k,n), &
                          saln( 1-nbdy,1-nbdy,k,m), &
                          uflux, &
                          vflux, &
                          util1,util2, &
                          pdzero, scp2,scp2i,delt1, btrmas)
      else   ! exactly isopycnal layer
        call advem(advtyp,saln( 1-nbdy,1-nbdy,k,n), &
                          saln( 1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdzero, scp2,scp2i,delt1, btrmas)
      endif
      do ktr= 1,ntracr
        if     (trcflg(ktr).eq.2) then !temperature tracer
          call advem(advtyp,tracer(1-nbdy,1-nbdy,k,n,ktr), &
                            tracer(1-nbdy,1-nbdy,k,m,ktr), &
                            uflx(  1-nbdy,1-nbdy,k), &
                            vflx(  1-nbdy,1-nbdy,k), &
                            util1,util2, &
                            pdtemp, scp2,scp2i,delt1, btrmas)
        else
!         ldebug_advem = k.eq.kk !debug layer kk tracer advection
          call advem(advtyp,tracer(1-nbdy,1-nbdy,k,n,ktr), &
                            tracer(1-nbdy,1-nbdy,k,m,ktr), &
                            uflx(  1-nbdy,1-nbdy,k), &
                            vflx(  1-nbdy,1-nbdy,k), &
                            util1,util2, &
                            pdzero, scp2,scp2i,delt1, btrmas)
!         ldebug_advem = .false. !turn off debugging of advection
        endif !trcflg
      enddo !ktr
      if (mxlmy) then
        call advem(advtyp,q2(   1-nbdy,1-nbdy,k,n), &
                          q2(   1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdq2,   scp2,scp2i,delt1, btrmas)
        call advem(advtyp,q2l(  1-nbdy,1-nbdy,k,n), &
                          q2l(  1-nbdy,1-nbdy,k,m), &
                          uflx( 1-nbdy,1-nbdy,k), &
                          vflx( 1-nbdy,1-nbdy,k), &
                          util1,util2, &
                          pdq2,   scp2,scp2i,delt1, btrmas)
      endif
!
      if     (lpipe .and. lpipe_tsadvc) then
! ---   compare two model runs.
        write (text,'(a9,i3)') 'temp.n k=',k
        call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'saln.n k=',k
        call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
        write (text,'(a9,i3)') 'th3d.n k=',k
        call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
      endif
!
!diag if (itests.gt.0.0and.jtests.gt.0) &
!diag write (lp,'(i9,2i5,i3,'' th,s,dp after advection  '',2f9.3,f8.2)') &
!diag nstep,itests,jtests,k,temp(itests,jtests,k,n),saln(itests,jtests,k,n), &
!diag dp(itests,jtests,k,n)*qonem
!
      if     (mod(nstep,3).eq.0 .or. diagno) then
!$OMP   PARALLEL DO PRIVATE(j,i) &
!$OMP            SCHEDULE(STATIC,jblk) !NOCSD
        do j=1,jj
          sminny(j)= 999.  !simplifies OpenMP parallelization
          smaxxy(j)=-999.  !simplifies OpenMP parallelization
!DIR$     PREFERVECTOR
          do i=1,ii
            if (SEA_P) then
              if     (dp(i,j,k,n).gt.onemm) then
                sminny(j)=min(sminny(j),saln(i,j,k,n))
                smaxxy(j)=max(smaxxy(j),saln(i,j,k,n))
              endif
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO !NOCSD
        xmin(k) = minval(sminny(1:jj))
        xmax(k) = maxval(smaxxy(1:jj))
      endif !every 3 time steps or diagno
!
 81   continue  ! k=1,kk
!
      call pipe_comparall(m,n, 'advem,  step')
!
! --- check for negative scalar fields.
!
      if     (mod(nstep,3).eq.0 .or. diagno) then
        call xcminr(xmin(1:kk))
        call xcmaxr(xmax(1:kk))
!
        do k= 1,kk
          sminn=xmin(k)
          smaxx=xmax(k)
!
          if (sminn.lt.0.0) then
            do j=1,jj
              do i=1,ii
                if (SEA_P) then
                  if (saln(i,j,k,n).eq.sminn) then
                    write (lp,'(i9,a,2i5,i3,a,f10.2)')  &
                      nstep,' i,j,k =',i+i0,j+j0,k, &
                      ' neg. saln after advem call ', &
                      saln(i,j,k,n)
                  endif !sminn
                endif !ip
              enddo !i
            enddo !j
            call xcsync(flush_lp)
          endif
!
          if (diagno) then
            if     (mnproc.eq.1) then
            if     (sminn.le.smaxx) then
              write (lp,'(i9,i3, a,2f9.3, a,1pe9.2,a)') &
                nstep,k, &
                ' min/max of s after advection:',sminn,smaxx, &
                '   (range:',smaxx-sminn,')'
            else
              write (lp,'(i9,i3, a,a)') &
                nstep,k, &
                ' min/max of s after advection:',' N/A (thin layer)'
            endif !normal:thin layer
            call flush(lp)
            endif
          endif
        enddo !k
      endif !every 3 time steps or diagno
!
! --- --------------------------------------
! --- diffusion of thermodynamic variable(s)
! --- --------------------------------------
!
      if     (temdf2.gt.0.0) then
        mdf = 2  !Laplacian
        call xctilr(saln(   1-nbdy,1-nbdy,1,n),1,kk, mdf,mdf, halo_ps)
        call xctilr(temp(   1-nbdy,1-nbdy,1,n),1,kk, mdf,mdf, halo_ps)
        call xctilr(th3d(   1-nbdy,1-nbdy,1,n),1,kk, mdf,mdf, halo_ps)
        if (mxlmy) then
          call xctilr(q2( 1-nbdy,1-nbdy,0,n),1,kk+2, mdf,mdf, halo_ps)
          call xctilr(q2l(1-nbdy,1-nbdy,0,n),1,kk+2, mdf,mdf, halo_ps)
        endif
        do ktr= 1,ntracr
          call xctilr(tracer(1-nbdy,1-nbdy,1,n,ktr), &
                                               1,kk, mdf,mdf, halo_ps)
        enddo !ktr
!
        do k=1,kk
!
! ---     for isopycnic vertical coordinates:
! ---       diffuse -th3d- and -saln- in the mixed layer (layer 1),
! ---       diffuse            -saln- only in all other layers
! ---     for hybrid vertical coordinates:
! ---       diffuse -saln- in all layers
! ---       diffuse -temp- in all layers if temdfc < 0.0
! ---       diffuse -th3d- in all layers if temdfc < 1.0
! ---       if 0.0 < temdfc < 1.0:
! ---         combine -temp- and -th3d- diffusion in density space
!
          ldtemp =  k.le.nhybrd .and. temdfc.gt.0.0     ! diffus temp
          ldth3d = (k.le.nhybrd .and. temdfc.lt.1.0) .or. &
                   (k.eq.1      .and. isopyc     )      ! diffus th3d
          if     (ldtemp .and. ldth3d) then ! diffus temp and th3d
            call tsdff_2x(th3d(1-nbdy,1-nbdy,k,n), &
                          temp(1-nbdy,1-nbdy,k,n))
            call tsdff_1x(saln(1-nbdy,1-nbdy,k,n))
          elseif (ldtemp) then ! diffus temp
            call tsdff_2x(temp(1-nbdy,1-nbdy,k,n), &
                          saln(1-nbdy,1-nbdy,k,n))
          elseif (ldth3d) then ! diffus th3d
            call tsdff_2x(th3d(1-nbdy,1-nbdy,k,n), &
                          saln(1-nbdy,1-nbdy,k,n))
          else   ! exactly isopycnal layer
            call tsdff_1x(saln(1-nbdy,1-nbdy,k,n))
          endif
          if (mxlmy) then
            call tsdff_2x(q2(  1-nbdy,1-nbdy,k,n), &
                          q2l( 1-nbdy,1-nbdy,k,n))
          endif !mxlmy
          do ktr= 1,ntracr,2
            if     (ktr+1.le.ntracr) then
              call tsdff_2x(tracer(1-nbdy,1-nbdy,k,n,ktr), &
                            tracer(1-nbdy,1-nbdy,k,n,ktr+1))
            else
!             ldebug_tsdif = k.eq.kk !debug layer kk tracer diffusion
              call tsdff_1x(tracer(1-nbdy,1-nbdy,k,n,ktr))
!             ldebug_tsdif = .false. !turn off debugging of diffusion
            endif
          enddo !ktr
        enddo !k
!
!       non-independent thermodynamic variable
!
        margin = 0
!$OMP   PARALLEL DO PRIVATE(j,k,i,ldtemp,ldth3d,th3d_t) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do k=1,kk
            ldtemp =  k.le.nhybrd .and. temdfc.gt.0.0     ! diffus temp
            ldth3d = (k.le.nhybrd .and. temdfc.lt.1.0) .or. &
                     (k.eq.1      .and. isopyc     )      ! diffus th3d
            do i=1-margin,ii+margin
              if (SEA_P) then
                if     (ldtemp .and. ldth3d) then
! ---             combine -temp- and -th3d- diffusion in density space
                  th3d_t       =sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                  th3d(i,j,k,n)=(1.0-temdfc)*th3d(i,j,k,n) + &
                                     temdfc *th3d_t
                  temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase, &
                                       saln(i,j,k,n))
                elseif (ldtemp) then
                  th3d(i,j,k,n)=sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                elseif (ldth3d) then
                  temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase, &
                                       saln(i,j,k,n))
                else   ! exactly isopycnal layer
                  th3d(i,j,k,n)=theta(i,j,k)
                  temp(i,j,k,n)=tofsig(th3d(i,j,k,n)+thbase, &
                                       saln(i,j,k,n))
                endif !ld....
              endif !ip
            enddo !i
          enddo !k
        enddo !j
!$OMP   END PARALLEL DO
      endif !temdf2.gt.0.0
!
      do k=1,kk
        if     (lpipe .and. lpipe_tsadvc) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'util1  k=',k
          call pipe_compare_sym1(util1,ip,text)
          write (text,'(a9,i3)') 'util2  k=',k
          call pipe_compare_sym1(util2,ip,text)
          write (text,'(a9,i3)') 'temp.n k=',k
          call pipe_compare_sym1(temp(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'saln.n k=',k
          call pipe_compare_sym1(saln(1-nbdy,1-nbdy,k,n),ip,text)
          write (text,'(a9,i3)') 'th3d.n k=',k
          call pipe_compare_sym1(th3d(1-nbdy,1-nbdy,k,n),ip,text)
        endif
!
!diag   if (itests.gt.0.and.jtests.gt.0) then &
!diag     write (lp,'(i9,2i5,i3,a,2f9.3,f8.2)') &
!diag       nstep,itests+i0,jtests+j0,k, &
!diag       ' t,s,dp after isopyc.mix.', &
!diag       temp(itests,jtests,k,n),saln(itests,jtests,k,n), &
!diag       dp(itests,jtests,k,n)*qonem &
!diag     call flush(lp) &
!diag   endif
!
      enddo !k

      return
 
      contains
!
      subroutine tsdff_2x(fld1,fld2)
      real fld1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
           fld2(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
! --- Laplacian diffusion for two scalar fields
!
        margin = 1
!
!$OMP   PARALLEL DO PRIVATE(j,i,factor) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              factor=temdf2*aspux(i,j)* &
                     scuy(i,j)*harmonc(dp(i-1,j,k,n)*onetamas(i-1,j,n) &
                                      ,dp(i  ,j,k,n)*onetamas(i,  j,n))
              uflux (i,j)=factor*(fld1(i-1,j)-fld1(i,j))
              uflux2(i,j)=factor*(fld2(i-1,j)-fld2(i,j))
            endif !iu
            if (SEA_V) then
              factor=temdf2*aspvy(i,j)* &
                     scvx(i,j)*harmonc(dp(i,j-1,k,n)*onetamas(i,j-1,n) &
                                      ,dp(i,j  ,k,n)*onetamas(i,j  ,n))
              vflux (i,j)=factor*(fld1(i,j-1)-fld1(i,j))
              vflux2(i,j)=factor*(fld2(i,j-1)-fld2(i,j))
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        if     (lpipe .and. lpipe_tsadvc) then
! ---     compare two model runs.
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux, iu,textu, &
                                 vflux, iv,textv)
          write (textu,'(a9,i3)') 'uflux2 k=',k
          write (textv,'(a9,i3)') 'vflux2 k=',k
          call pipe_compare_sym2(uflux2,iu,textu, &
                                 vflux2,iv,textv)
        endif
!
! ---   rhs: dp.n, uflux+, vflux+, uflux2+, vflux2+
! ---   lhs: saln.n, temp.n, th3d.n
!
        margin = 0
!
!$OMP   PARALLEL DO PRIVATE(j,i,factor) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              factor=-delt1/(scp2(i,j)*max(dp(i,j,k,n)*onetamas(i,j,n), &
                                           eps_har))
              util1(i,j)=((uflux (i+1,j)-uflux (i,j)) &
                         +(vflux (i,j+1)-vflux (i,j)))*factor
              fld1(i,j)=fld1(i,j)+util1(i,j)
              util2(i,j)=((uflux2(i+1,j)-uflux2(i,j)) &
                         +(vflux2(i,j+1)-vflux2(i,j)))*factor
              fld2(i,j)=fld2(i,j)+util2(i,j)
!
!diag         if (i.eq.itests.and.j.eq.jtests) then
!diag           if (1.le.i .and. i.le.ii .and. &
!diag               1.le.j .and. j.le.jj      ) then &
!diag             write (lp,100) nstep,i+i0,j+j0,k,'t,s,dt,ds', &
!diag             fld1(i,j),fld2(i,j),util1(i,j),util2(i,j) &
!diag             call flush(lp) &
!diag           endif &
!diag         endif
!
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        if     (lpipe .and. lpipe_tsadvc) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'util1  k=',k
          call pipe_compare_sym1(util1,ip,text)
          write (text,'(a9,i3)') 'util2  k=',k
          call pipe_compare_sym1(util2,ip,text)
          write (text,'(a9,i3)') 'fld1.n k=',k
          call pipe_compare_sym1(fld1(1-nbdy,1-nbdy),ip,text)
          write (text,'(a9,i3)') 'fld2.n k=',k
          call pipe_compare_sym1(fld2(1-nbdy,1-nbdy),ip,text)
        endif
!
      return
      end subroutine tsdff_2x
!
      subroutine tsdff_1x(fld1)
      real fld1(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
!
! --- Laplacian  diffusion for a single scalar field
!
        margin = 1
!
!$OMP   PARALLEL DO PRIVATE(j,i,factor) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_U) then
              factor=temdf2*aspux(i,j)* &
                     scuy(i,j)*harmonc(dp(i-1,j,k,n)*onetamas(i-1,j,n) &
                                      ,dp(i  ,j,k,n)*onetamas(i,  j,n))
              uflux (i,j)=factor*(fld1(i-1,j)-fld1(i,j))
            endif !iu
            if (SEA_V) then
              factor=temdf2*aspvy(i,j)* &
                     scvx(i,j)*harmonc(dp(i,j-1,k,n)*onetamas(i,j-1,n) &
                                      ,dp(i,j  ,k,n)*onetamas(i,j  ,n))
              vflux (i,j)=factor*(fld1(i,j-1)-fld1(i,j))
            endif !iv
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        if     (lpipe .and. lpipe_tsadvc) then
! ---     compare two model runs.
          write (textu,'(a9,i3)') 'uflux  k=',k
          write (textv,'(a9,i3)') 'vflux  k=',k
          call pipe_compare_sym2(uflux, iu,textu, &
                                 vflux, iv,textv)
        endif
!
      if     (ldebug_tsdif .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        if (SEA_U) then
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: u dp ',i+i0,j+j0, &
            uflux(i,j),  dp(i-1,j,k,n)*qonem,dp(i  ,j,k,n)*qonem
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: u fld',i+i0,j+j0, &
            uflux(i,j),fld1(i-1,j    ),    fld1(i  ,j    )
        else
          write (lp,'(a,2i5,1e16.6)') &
            'tsdff: u LND',i+i0,j+j0, &
            uflux(i,j)
        endif
        if (SEA_V) then
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: v dp ',i+i0,j+j0, &
            vflux(i,j),  dp(i,j-1,k,n)*qonem,dp(i,j  ,k,n)*qonem
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: v fld',i+i0,j+j0, &
            vflux(i,j),fld1(i,j-1    ),    fld1(i,j      )
        else
          write (lp,'(a,2i5,1e16.6)') &
            'tsdff: v LND',i+i0,j+j0, &
            vflux(i,j)
        endif
        i=itests+1
        j=jtests
        if (SEA_U) then
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: u dp ',i+i0,j+j0, &
            uflux(i,j),  dp(i-1,j,k,n)*qonem,dp(i  ,j,k,n)*qonem
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: u fld',i+i0,j+j0, &
            uflux(i,j),fld1(i-1,j    ),    fld1(i  ,j    )
        else
          write (lp,'(a,2i5,1e16.6)') &
            'tsdff: u LND',i+i0,j+j0, &
            uflux(i,j)
        endif
        i=itests
        j=jtests+1
        if (SEA_V) then
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: v dp ',i+i0,j+j0, &
            vflux(i,j),  dp(i,j-1,k,n)*qonem,dp(i,j  ,k,n)*qonem
          write (lp,'(a,2i5,3e16.6)') &
            'tsdff: v fld',i+i0,j+j0, &
            vflux(i,j),fld1(i,j-1    ),    fld1(i,j      )
        else
          write (lp,'(a,2i5,1e16.6)') &
            'tsdff: u LND',i+i0,j+j0, &
            vflux(i,j)
        endif
      endif
!
! ---   rhs: dp.n, uflux+, vflux+, uflux2+, vflux2+
! ---   lhs: saln.n, temp.n, th3d.n
!
        margin = 0
!
!$OMP   PARALLEL DO PRIVATE(j,i,factor) &
!$OMP            SCHEDULE(STATIC,jblk)
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              factor=-delt1/(scp2(i,j)*max(dp(i,j,k,n)*onetamas(i,j,n), &
                                           eps_har))
              util1(i,j)=((uflux (i+1,j)-uflux (i,j)) &
                         +(vflux (i,j+1)-vflux (i,j)))*factor
              fld1(i,j)=fld1(i,j)+util1(i,j)
!
!diag         if (i.eq.itests.and.j.eq.jtests) then
!diag           if (1.le.i .and. i.le.ii .and. &
!diag               1.le.j .and. j.le.jj      ) then &
!diag             write (lp,100) nstep,i+i0,j+j0,k,'t,s,dt,ds', &
!diag             fld1(i,j),0.0,util1(i,j),0.0 &
!diag             call flush(lp) &
!diag           endif &
!diag         endif
!
            endif !ip
          enddo !i
        enddo !j
!$OMP   END PARALLEL DO
!
        if     (lpipe .and. lpipe_tsadvc) then
! ---     compare two model runs.
          write (text,'(a9,i3)') 'util1  k=',k
          call pipe_compare_sym1(util1,ip,text)
          write (text,'(a9,i3)') 'fld1.n k=',k
          call pipe_compare_sym1(fld1(1-nbdy,1-nbdy),ip,text)
        endif
!
      if     (ldebug_tsdif .and. itests.gt.0 .and. jtests.gt.0) then
        i=itests
        j=jtests
        write (lp,'(a,2i5,3e16.6)') &
          'tsdff: p fld',i+i0,j+j0, &
          fld1(i,j),util1(i,j), &
          -delt1/(scp2(i,j)*max(dp(i,j,k,n)*onetamas(i,j,n),eps_har))
      endif
!
      return
      end subroutine tsdff_1x

!-----end contains

      end subroutine tsadvc
!
      end module mod_tsadvc
!
!  Revision history:
!
!> June 1995 - eliminated setting of salinity in massless layers (loop 46)
!>             (this is now done in mxlayr.f)
!> Aug. 1995 - omitted t/s/dp time smoothin, case of abrupt mxlayr.thk.change
!> Sep. 1995 - increased temdf2 if mixed layer occupies >90% of column
!> Mar. 2000 - removed 'cushn' and added logic to assure global conservation
!> Apr. 2000 - conversion to SI units
!> Apr. 2000 - changed i/j loop nesting to j/i
!> Aug. 2000 - temp advection and diffusion only for hybrid vertical coordinate
!> Nov. 2000 - nhybrd T&S advection layers, kdm-nhybrd S advection layers
!> Nov. 2000 - T&S or th&S advection/diffusion based on advflg
!> Feb. 2001 - placed advem in a module
!> May  2002 - diffusion coefficent based on max(sc?x,sc?y)
!> Aug. 2003 - separate PCM and MPDATA versions (advtyp)
!> Aug. 2003 - added FCT2 and UTOPIA advection options.
!> Nov. 2003 - per layer diffusion routine for 1 or 2 scalar fields
!> Feb. 2008 - fixed famin,famax bugs in FCT2/4
!> Jun. 2008 - fixed q2,q2l halo update bug
!> Aug  2011 - reworked Robert-Asselin filter, RA of dp now in cnuity
!> Aug  2012 - RA filter now exactly conserves constant tracers
!> Sep  2012 - RA filter not applied if thickness < onezm
!> Apr. 2014 - replace ip with ipa for mass sums
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> May  2014 - use ipim1,ipip1,ipjm1,ipjp1 as sea-only neighbors
!> Aug. 2018 - Robert-Asselin filter now in mod_asselin
!> Aug. 2018 - btrmas added, use onetamas to simplify logic
!> Aug. 2018 - replaced itest,jtest with itests,jtests
!> Nov. 2018 - always use oneta for diffusion
