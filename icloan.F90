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
      subroutine icloan(m,n)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
      implicit none
!
      integer m,n
!
! --- 'energy loan' ice model. no advection, no dynamics. ice amount
! --- represents energy 'loaned' to water column to prevent wintertime
! --- cooling below freezing level. loan is paid back in summer.
!
! --- modified version for ice-ocean "coupling".
! --- freeze/melt energy from relaxation to the freezing temperature.
!
      integer i,j
      real    tfrz,tsur,tmxl,smxl,hfrz,paybak,borrow,hice,thkimx,t2f
      real    radfl,swfl,tdif,wind,airt,pair,rair,snsibl,emnp,dtrmui
      real    thkimxy(jdm)
!
! --- hice   = actual ice thickness (m), local variable
!
! --- thkice = average ice thickness, i.e. hice x covice (m)
! --- covice = ice coverage, i.e. cell fraction (0.0 to 1.0)
! --- temice = ice surface temperature           (degC)
! --- flxice = cell average heat  flux under ice (W/m^2)        into ocean
! --- fswice = cell average swv   flux under ice (W/m^2)        into ocean
! --- sflice = cell average salt  flux under ice (psu kg/m^2/s) into ocean
! --- wflice = cell average water flux under ice (    kg/m^2/s) into ocean
! --- wflfrz = cell average water flux under ice due to freezing/melting
!
! --- icefrq = e-folding time scale back to tfrz (no. time steps)
! --- thkfrz = maximum thickness of near-surface freezing zone (m)
! --- tfrz_0 = ice melting point (degC) at S=0psu
! --- tfrz_s = gradient of ice melting point (degC/psu)
! --- ticegr = vertical temperature gradient inside ice (deg/m)
! ---            (0.0 to get ice surface temp. from atmos. surtmp)
! --- hicemn = minimum ice thickness (m)
! --- hicemx = maximum ice thickness (m)
!
      real       tfrz_n,ticemn,ticemx,salice,rhoice,fusion,meltmx
      parameter (tfrz_n= -1.79,  & ! nominal ice melting point (degC)
                 ticemn=-50.0,   & ! minimum ice surface temperature (degC)
                 ticemx=  0.0,   & ! maximum ice surface temperature (degC)
                 salice=  4.0,   & ! salinity of ice (psu) - same as CICE
                 rhoice=917.0,   & ! density  of ice (kg/m**3)
                 fusion=334.e3,  & ! latent heat of fusion (J/kg)
                 meltmx= 33.e-7)! max. ice melting rate (m/s), 0.285 m/day
!
      real       fluxmx         !max. ice melting flux (W/m^2)
      parameter (fluxmx=meltmx*fusion*rhoice)    !~1000 W/m^2 - like CICE
!
      real       csice,csubp,pairc,rgas,tzero
      parameter (csice =0.0006,       & !ice-air sensible exchange coefficient
                 csubp =1005.7,       & !specific heat of air (j/kg/deg)
                 pairc=1013.0*100.0,  & !air pressure (mb) * 100
                 rgas =287.1,         & !gas constant (j/kg/k)
                 tzero=273.16)       !celsius to kelvin temperature offset
!
      real       sb_cst
      parameter (sb_cst=5.67e-8)     !Stefan-Boltzman constant, for lwflag=-1
!
!#     include "stmt_fns.h"
!
      dtrmui = baclin/(1.0*86400.0)  !dt*1/1days
!
! --- energy loan: add extra energy to the ocean to keep SST from dropping
! --- below tfrz in winter. return this borrowed energy to the 'energy bank'
! --- in summer.
!
! --- salt loan: analogous to energy loan.
!
!$OMP PARALLEL DO PRIVATE(j,i, &
!$OMP                     t2f,hfrz,smxl,tmxl,tfrz,borrow,paybak) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        thkimxy(j)=0.0 !simplifies OpenMP parallelization
        do i=1,ii
          if (SEA_P) then
          if     (ishlf(i,j).eq.0) then  !under an ice shelf
            flxice(i,j)=0.0
            wflice(i,j)=0.0
            sflice(i,j)=0.0
            thkice(i,j)=hicemx
             util1(i,j)=0.0
          else !standard ocean point
! ---       relax to tfrz with e-folding time of icefrq time steps
! ---       assuming the effective surface layer thickness (hfrz) 
! ---       is at most thkfrz meters
! ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
! ---       icefrq==1 is "Instant Relaxation" and an e-folding time 
! ---       of 30 days is consistent with the "Drag Law" approach
! ---         D.M. Holland (1998) On the Parameterization of Basal
! ---          Heat Flux for Sea-ice Modelling; Geophysica 34 pp 1-21
! ---       when coupling, icefrq must be at least the coupling interval
!
            hfrz   = min( thkfrz*onem, dpbl(i,j) )
            t2f    = (spcifh*hfrz)/(baclin*icefrq*g)
            smxl   = saln(i,j,1,n)
            tmxl   = temp(i,j,1,n)
            tfrz   = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
            borrow = (tfrz-tmxl)*t2f       !W/m^2 into ocean
!
! ---       limit heat flux range (for both forming and melting ice)
            borrow=max( -fluxmx, min( fluxmx, borrow ) )
!
!diag       if (i.eq.itest .and. j.eq.jtest) then
!diag         write (lp,'(i9,2i5,a,5f9.3)') &
!diag           nstep,i+i0,j+j0,'  t,tfrz,flx,hfrz,cov:', &
!diag           tmxl,tfrz,borrow,hfrz*qonem,covice(i,j)
!diag       endif
!
            if (tmxl.lt.tfrz) then
!
! ---         add energy to move tmxl towards tfrz (only if tmxl < tfrz)
! ---         include some dependance on sea ice coverage
!
              flxice(i,j)=borrow*max(covice(i,j),0.1)               !+ve
              wflice(i,j)=           -flxice(i,j)/fusion            !-ve
              thkice(i,j)=thkice(i,j)-wflice(i,j)*(baclin/rhoice)   !+ve inc.
! ---         brine rejection as ice forms
              sflice(i,j)=            wflice(i,j)*min(smxl,salice)  !-ve
            elseif (thkice(i,j).gt.0.0) then  !tmxl > tfrz
!
! ---         ice, so return the borrowed amount whenever tmxl > tfrz
! ---         only over the fraction of the cell where ice exists
!
              paybak=min( -borrow*covice(i,j),                       & !+ve
                          thkice(i,j)*(fusion*rhoice/baclin) )
              flxice(i,j)=           -paybak                        !-ve
              wflice(i,j)=           -flxice(i,j)/fusion            !+ve
              thkice(i,j)=thkice(i,j)-wflice(i,j)*(baclin/rhoice)   !-ve inc.
! ---         brine recovery from melting ice
              sflice(i,j)=            wflice(i,j)*salice            !+ve
            else !tmxl > tfrz & thkice(i,j) == 0.0
!
! ---         no ice.
!
              flxice(i,j)=0.0
              wflice(i,j)=0.0
              sflice(i,j)=0.0
!
              if (icmflg.eq.2) then
!
! ---           add extra cooling under the ice mask (tsur<=tfrz_n)
! ---           don't allow a new tsur maximum, to preserve sea ice
!
                if     (natm.eq.2) then
                  tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1) ), &
                              surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1   )
                elseif (yrflag.lt.2) then
                  tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1), &
                                   surtmp(i,j,l2), surtmp(i,j,l3) ), &
                              surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1+ &
                              surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3   )
                else
                  tsur = min( max( surtmp(i,j,l0), surtmp(i,j,l1) ), &
                              surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1   )
                endif
                if     (tsur.le.tfrz_n) then
                  surflx(i,j)=surflx(i,j)+borrow*covice(i,j)
                endif
              endif !icmflg.eq.2
            endif
!
            util1(i,j)=max(thkice(i,j)-hicemx,0.0)  !icex = ice exceeding hicemx
            thkimxy(j)=max(thkimxy(j),thkice(i,j))
          endif !ishlf:else
          wflfrz(i,j)=wflice(i,j)  !diagnostic only
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      thkimx=maxval(thkimxy(1:jj))
      call xcmaxr(thkimx)
!
! --- spread out portion of ice thicker than hicemx
      if (thkimx.gt.hicemx) then
        call psmooth(util1, 0,0, ishlf, util2)  !smooth icex
      endif
!
!$OMP PARALLEL DO PRIVATE(j,i,hice,smxl,tfrz, &
!$OMP                     radfl,tdif,wind,airt,rair,snsibl,emnp) &
!$OMP          SCHEDULE(STATIC,jblk)
      do j=1,jj
        do i=1,ii
          if (SEA_P) then
          if     (ishlf(i,j).eq.0) then  !under an ice shelf
            thkice(i,j)=hicemx
            covice(i,j)=1.0
            fswice(i,j)=0.0
            temice(i,j)=ticemn
          else !standard ocean point
            thkice(i,j)=util1(i,j)+min(thkice(i,j),hicemx) !icex_sm+rest
!
! ---       compute fractional ice coverage for energy flux calculation
            if (thkice(i,j).lt.1.e-5*hicemn) then
              covice(i,j)=0.0
            else
              covice(i,j)=min(1.0,thkice(i,j)*(1.0/hicemn))
              hice=thkice(i,j)/covice(i,j)  !minimum of hicemn
            end if
!
            if     (icmflg.eq.3) then
! ---         relax to sea ice concentration from coupler
! ---         ice thickness is therefore always then 0 and hicemn
              covice(i,j)=covice(i,j)+dtrmui*(si_c(i,j)-covice(i,j))
              thkice(i,j)=covice(i,j)*hicemn
              hice=hicemn
            endif
!
! ---       compute ice surface temperature
            if     (covice(i,j).eq.0.0) then
              temice(i,j)=ticemx
            elseif (icmflg.eq.3 .and. si_c(i,j).gt.0.0) then  !from coupler
              temice(i,j)=max( ticemn, min( ticemx, si_t(i,j) ) )
            elseif (ticegr.eq.0.0) then  !use surtmp
              if     (natm.eq.2) then
                temice(i,j)=max( ticemn, &
                                 min( ticemx, &
                                      surtmp(i,j,l0)*w0+ &
                                      surtmp(i,j,l1)*w1  ) )
              else
                temice(i,j)=max( ticemn, &
                                 min( ticemx, &
                                      surtmp(i,j,l0)*w0+ &
                                      surtmp(i,j,l1)*w1+ &
                                      surtmp(i,j,l2)*w2+ &
                                      surtmp(i,j,l3)*w3  ) )
              endif !natm
            else
              temice(i,j)=max( ticemn, ticemx-ticegr*hice )
            endif
!
! ---       atmosphere to ice surface exchange is applied to the ocean,
! ---        i.e. use the "energy-loan" approach.
! ---       don't apply to the ocean when coupling, because it has
! ---        already been applied to the ice.
!
            if     (icmflg.ne.3 .and. covice(i,j).gt.0.0) then
! ---         net radiative thermal flux (w/m**2) +ve into ocean/ice
              if     (lwflag.gt.-1) then
! ---           radflx's Qsw includes the atmos. model's surface albedo,
! ---           i.e. it already allows for ice&snow where it is observed.
                if     (natm.eq.2) then
                  radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
                else
                  radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1 &
                       +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
                endif !natm
                if    (lwflag.gt.1) then
! ---             longwave correction to radfl (Qsw+Qlw).
! ---             this will be ~zero for ticegr==0.0 (temice=surtmp)
                  if     (natm.eq.2) then
                    tdif = temice(i,j) - &
                         ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1 )
                  else
                    tdif = temice(i,j) - &
                         ( surtmp(i,j,l0)*w0+surtmp(i,j,l1)*w1 &
                          +surtmp(i,j,l2)*w2+surtmp(i,j,l3)*w3 )
                  endif !natm
                  !correction is blackbody radiation from tdif at temice
                  radfl = radfl - (4.506+0.0554*temice(i,j)) * tdif
                endif
              else  !lwflag.eq.-1
! ---           input radflx is Qlwdn
! ---           input  swflx is net Qsw with ocean albedo (~ 1-0.09)
! ---           convert  swfl to net Qsw with sea ice albedo (1-0.6)
! ---                      i.e. scale by(1-0.6)/(1-0.09)=0.44
! ---           convert radfl to net Qlw + Qsw where sea ice
                if     (natm.eq.2) then
                  radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1
                   swfl=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1
                else
                  radfl=radflx(i,j,l0)*w0+radflx(i,j,l1)*w1 &
                       +radflx(i,j,l2)*w2+radflx(i,j,l3)*w3
                   swfl=swflx (i,j,l0)*w0+swflx (i,j,l1)*w1 &
                       +swflx (i,j,l2)*w2+swflx (i,j,l3)*w3
                endif !natm
                radfl = radfl - sb_cst*(temice(i,j)+tzero)**4  &
                              + 0.44d0*swfl !!Alex
              endif  !lwflux
!
              if     (flxflg.ne.3) then
                if     (natm.eq.2) then
! ---             wind speed (m/s)
                  wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
! ---             air temperature (C)
                  airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1
                else
! ---             wind speed (m/s)
                  wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1 &
                      +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3
! ---             air temperature (C)
                  airt=airtmp(i,j,l0)*w0+airtmp(i,j,l1)*w1 &
                      +airtmp(i,j,l2)*w2+airtmp(i,j,l3)*w3
                endif !natm
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
                rair   = pair/(rgas*(tzero+airt))
                snsibl = csubp*rair*wind*csice*(temice(i,j)-airt)
              else
                snsibl = 0.0 !already in total flux (i.e. in radfl)
              endif
              flxice(i,j) = flxice(i,j) + &
                            covice(i,j)*(radfl - snsibl) !no evap
!
! ---         add a time-invarient net heat flux offset
              if     (flxoff) then
                flxice(i,j) = flxice(i,j) + covice(i,j)*offlux(i,j)
              endif
!
! ---         emnp = evaporation minus precipitation (m/sec) into atmos.
! ---         no evap (sublimation) over ice, all precip enters ocean
              if     (pcipf) then
                if     (natm.eq.2) then
                  emnp = -( precip(i,j,l0)*w0+precip(i,j,l1)*w1)
                else
                  emnp = -( precip(i,j,l0)*w0+precip(i,j,l1)*w1 &
                           +precip(i,j,l2)*w2+precip(i,j,l3)*w3)
                endif !natm
              else
                emnp =  0.0
              endif
! ---         wflice = water flux (m/s kg/m**2/sec) into ocean under ice
              wflice(i,j) = wflice(i,j) - covice(i,j)*emnp*rhoref
            endif !covice
            fswice(i,j) = 0.0 !no penetrating Qsw under ice
          endif !ishlf:else
          endif !ip
        enddo !i
      enddo !j
!$OMP END PARALLEL DO
!
      return
      end subroutine icloan
!
!
!> Revision history
!>
!> June 2000 - conversion to SI units
!> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
!> May  2003 - added option to impose an ice mask
!> June 2003 - added 8 time step e-folding time scale
!> June 2003 - limited rate of ice formation
!> June 2003 - replaced constant saldif with smxl-salice
!> Mar. 2005 - freezing point linearly dependent on salinity
!> Mar. 2005 - ice surface temperature optionally from surtmp
!> June 2006 - modified version for ice-ocean "coupling"
!> Nov. 2011 - don't apply atmosphere to ice surface exchange when "coupling"
!> May  2012 - limit brine rejection to be a non-negative salt flux
!> July 2012 - flxice and sflice now correctly represent cell average under ice
!> Nov. 2012 - weaker dependance on covice when freezing
!> Jan. 2014 - added natm
!> Apr. 2014 - added pair for time varying msl pressure (mslprf)
!> Apr. 2014 - added ice shelf logic (ishlf)
!> May  2014 - use land/sea masks (e.g. ip) to skip land
!> Aug. 2018 - added sflfrz (now wflfrz) as a diagnostic field
!> Aug. 2018 - icloan is not leapfrog (baclin, not delt1)
!> Nov. 2018 - virtual salt flux replaced with water and actual salt flux
!> Nov. 2018 - added lwflag=-1 for input radflx=Qlwdn
!> Nov. 2018 - allow for difference in ocean and sea ice albedo when lwflag=-1
