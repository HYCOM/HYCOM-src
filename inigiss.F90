      subroutine inigiss
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 2.1
      implicit none
!
! -----------------------------------------------
! --- initialize nasa-giss vertical mixing scheme
! -----------------------------------------------
!
      integer i,j,k
!
      real slq2b_00,smb_00,shb_00,ssb_00,c_y0,c_y00,deltanum, &
           deltaden,delta,rrcrn,rrcrp,theta_rcrn_deg,theta_rcrp_deg, &
           delra_r,theta_r_deg,theta_r,sm_r0,sh_r0,ss_r0,ra_r,ra_r1, &
           rit,ric,ri_r,rid_r,sm_r,sh_r,ss_r,smosloq_r,rit1,ric1, &
           ri_r1,rid_r1,slq2_r,smosloq_0r,ra_r0,rit0,ric0,ri_r0,rid_r0, &
           slq2_r0,c_y001,sisa1
!
      integer iridsign,iridstep,irid,iri,mt0s,mtm1s,idfs,idif,irisign, &
              iristep,itheta_r,jtheta_r,isailback,idifs,ira_r,ibg, &
              ipenra_r
!
      real    acosh1,xx
# include "stmt_fns.h"
      acosh1(xx) = log(xx+sqrt((xx**2)-1.0))
!
! --- initialize viscosity and diffusivity arrays
      do j=1,jdm
        do i=1,idm
          do k=1,kdm+1
            vcty(i,j,k)=diwm(i,j)
            dift(i,j,k)=diws(i,j)
            difs(i,j,k)=diws(i,j)
! --- no nonlocal forcing
            ghats(i,j,k)=0.0
          enddo
        enddo
      enddo
!
! --- dimensions for the tables in mxgiss routine
! --- the file is needed in order to preserve in the
! --- arrays for the Ri-tables
!
      pidbl=3.14159265358979312
      ntbl=251
      nextrtbl0=62
      ifexpabstable=1
      nextrtbl1=500
      nextrtbl=nextrtbl0+ifexpabstable*nextrtbl1
      nposapprox=51
      mt0=ntbl-nposapprox
      mt=mt0+nextrtbl       ! table dimensions
      mt_ra_r=nposapprox-1
      n_theta_r_oct=(((pidbl/4.)*mt_ra_r)/15.)*15
      deltheta_r=(pidbl/4.)/(n_theta_r_oct)
!
! --- set other parameters
!
      ri0=- 4.0           !C     parameter(ri0=-20.D0)
      ebase=2.71828182845904509
      ifback=5            !Temperature=Salt diffusivity model background 
                          !model swith. 
                          ! K_H,K_S (S=N/sqrt(Ri)),Ri=backfrac*Ri_Cr)
      ifsali=1            !Salinity model switch (Canuto's)
      ifepson2=2          !Background (epsilon/N^2) dimensionalization 
                          !of diffusivities switch.
                          !   cnst blw highst lvl frgr dies
      epson2_ref=.288     !reference value of dissipation/N**2
                          !Value of (epsilon/N^2)/(1 cm/sec^2) used. 
			  !See Canuto et al. JPO 2002 Sections 8&9.
		          !040126 Actual (epsilon/N^2) can vary with z,N and f .
      eps_bot0=2.e-5	  !The value of epsilon at the bottom in cgs, 
			  !St.Laurent et al. JPO2001 give epsilon = 3to9e-9 W/kg
			  !for slopes and 2to5e-9 W/kg for crests and canyons.
      scale_bot=5.e-4	  !The scale (in cm) of exponential decrease of mixing 
			  !above the bottom with height. St. Laurent et al. give
			  !150+-50 m for slopes, 500+-100 m for crests and canyons.
      eplatidepmin=7.E-2  !Gregg et al. admit their formula eq.(2) for the 
			  !latitude dependent factor L which scales turbulence 
			  !won't work at the equator where it predicts epsilon=0.
   			  !Introduce eplatidepmin, a minimum on the factor L .
      wave_30=(pi/43082.0)*acosh1(5.24e-3/(pi/43082.0))
                          !reference value at 30degN with N=5.24e-3
                          !from Garerett and Munk, as used by Gregg et. al.
      ifrafgmax=1         !Switch for limiting BackGround ra_r 
                          !to at most Foreground ra_r when Ri>0
                          !for R_r in the [R_r_crit_DoubleDiffusion,
                          !R_r_crit_SaltFingers] regime.
      ifsalback=5         !Salinity background modification switch.
                          !int.wvS=N/(Ri_i^(1/2)),Ri_icnst,
                          !   ra_r_i=cnst*ra_r_crit.(theta_r)
      ifchengcon=0        !old ocean cnsts,near-surf prof assump
      ifpolartablewrite=0 !Switch to write out polar 2D 
                          !turbulence table . 
      ifbg_theta_interp=1 !Introduce flag for use of \theta_r 
                          !arrays to interpolate background.
                          !Intrplt 2D array
                          !(slq2_r1=array for (Sl/q)^2)
                          !with (Ri,Ri_d)indices
      back_ph_0=(6.e-5)*(1.e2/(2.e0*pidbl))
                          !for ifsalback=3 case.
                          !Gargett et. al. JPO Vol.11 p.1258-71 gives 
                          !for "the deep record",
                          !\phi_0=6\times10^{-5}s^{-2}cpm^{-1}. 
                          !"cpm" is 'cycles per meter'.
                          !\phi_0=6\times10^{-5}s^{-2}(2 pidbl/100)^{-1}cm
      adjust_gargett=1.0  !Gargett et. al. favor the value, 
                          !k_0 = 0.1 cpm. But k_0=0.05-0.2 cpm 
                          !might be viable, see section 5 of their 
                          !paper. Take k_0 = 0.1 cpm * adjust_gargett, 
                          !where adjust_gargett is adjustable.
                          !Convert to radians per cm: 
                          !k_0 = 0.1 (2pi/100cm) * adjust_gargett.
                          !used for ifsalback=4 case also, but set 
                          !adjust_gargett=1 for ifsalback=4 
      back_k_0=(0.1)*(2.0)*pidbl*(1.e-2)*adjust_gargett
                          !Introduce the lengthscale 
                          !\Delta_0 \equiv pi/k_0 .
                          !The units of \Delta_0 are centimeters, 
                          !with k_0 in radians per cm.
                          !`min turb' wvnmbr (cm^-1)
      back_del_0=pidbl/back_k_0
      back_s2=1.e-14	  !back_s2 should be smaller than any normal Shear^2
      back_sm2=1.0/back_s2       !1/back_s2 (sec^2)
      ri_internal=1.0     !Parameter for ifsalback=4 case.
      backfrac = 85.e-2   !Parameter for ifback or ifsalback=5 case.
                          !ifback=5:   =cnst frac{Ri_cr};
                          !ifsalback=5:=cnst frac{ra_r_crit.(\theta_r)}
      backfact = ebase**(-1)  !Parameter for ifsalback=6 case.
      ako = 1.6           ! Kolmogorov's constant
!
      tpvot0 = 0.4        ! \tau_pv = {2 \over 5} \tau     (B.1)
                          ! "tpv/tau" = 2/5
                          ! From the printed notes Canuto
                          ! gave Armando on 980601 have:
      sgmt=0.72           !Make "sgmt" a parameter.
                          !Standard value was 0.72.      
      tptot0=(1.0/5.0)*(1.0/(1.0+(1.0/sgmt)))
                          ! \tau_p\theta over \tau
      tpcot0=tptot0       !tau_pc over \tau
      ttot0=sgmt          !tau_\theta over \tau
      tcot0=ttot0         !tau_c over \tau
      tctot0=1.0/3.0      ! tau_c\theta } over \tau
      tpvot = tpvot0
      tptot = tptot0
      tpcot = tpcot0
      ttot  = ttot0
      tcot  = tcot0
      tctot = tctot0
!
        if (mnproc.eq.1) then
        write(lp,900)
 900    format('nasa-giss mixed layer model selected'/ &
       'turbulence calculated by 040128 hycom version'/ &
       'stripped down from 030803 turb_2gi1a ncar')
        endif !1st proc
!
! --- START OF SALINITY MODEL BACKGROUND LENGTHSCALE CALCULATION SECTION.
! --- ifsali.eq. 1  therefore:   
! --- Calculate constant lengthscale for 
! --- the background for ifsalback=3,4,5
! --- \Delta_0 ={B_1 pi \over (3 Ko)^{3/2}} l_0
! --- l_0 = {(3 Ko)^{3/2} \over B_1 pi} \Delta_0
! --- "back_l_0" is the constant background 
! --- l_0 in centimeters.
!
! --- pass back B_1 from oursal2. 
      call oursal2_1a(0.,0.,slq2b_00,smb_00,shb_00,ssb_00, &
                      c_y0,c_y00,0,0)
!
      back_l_0 = (((3.*ako)**(3./2.))/(b1*pi))*back_del_0
!
        if (mnproc.eq.1) then
        write(lp,*) "Dubovikov Internal wave constants for background."
        write(lp,*) "Ratio of Background to Critical ra_r"// &
                       " [\\equiv ({Ri_T}^2 + {Ri_C}^2)^(1/2)]",backfrac
        write(lp,*) "Lengthscale, del_0/(cm) =",back_del_0
        write(lp,*) "Lengthscale, l_0/(cm) =",back_l_0
        call flush(lp)
        endif !1st proc
!
! --- Set step-size for *both* dimensions of 2D table here.
! --- ifsali.eq. 1
      dri = -ri0/real(mt0)
!
! --- BUILD SALINITY MODEL TABLES VS. "Ri = Ri_T + Ri_C" AND "Ri_d = Ri_T - Ri_C".
! --- Use separate loops for calculation of independent table variables.
!
      do iridsign=0,1
      iridstep=(-1)**iridsign 
      do irid= 0,mt*iridstep,iridstep
! --- Set Ri_d table values. (See NBP59,63=p#A27,30.)
           if(abs(irid).le.mt0) then
             ridb(irid) = real(irid)*dri
           else
             mt0s = mt0*iridstep
             mtm1s = (mt0-1)*iridstep
! --- introduction of exponential absolute val table option.
           if(ifexpabstable.eq. 0) then
               idifs = (abs(irid)-mt0)*iridstep
               ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                            ridb(mtm1s))**(idifs**2))
           else if(ifexpabstable.eq. 1) then
               idif = abs(irid)-mt0
               ridb(irid) = ridb(mt0s)*((ridb(mt0s)/ &
                            ridb(mtm1s))**(idif))
           endif
           endif
!
      enddo 
      enddo
!
      do irisign=0,1
      iristep=(-1)**irisign 
      do iri= 0,mt*iristep,iristep
! --- Set Ri table values. (See NBP59,63=p#A27,30.)
           if(abs(iri).le.mt0) then
             ribtbl(iri) = real(iri)*dri
           else
             mt0s = mt0*iristep
             mtm1s = (mt0-1)*iristep
! --- introduction of exponential absolute val table option.
           if(ifexpabstable.eq. 0) then
               idifs = (abs(iri)-mt0)*iristep
               ribtbl(iri) = ribtbl(mt0s)*((ribtbl(mt0s)/ &
                            ribtbl(mtm1s))**(idifs**2))
           else if(ifexpabstable.eq. 1) then
               idif = abs(iri)-mt0
               ribtbl(iri) = ribtbl(mt0s)*((ribtbl(mt0s)/ &
                            ribtbl(mtm1s))**(idif))
           endif
           endif
!
      enddo 
      enddo
!
! --- If using interp2d_expabs introduce ratio between adjacent Richardson
! --- numbers in nonlinear part of table.***
        rri = ribtbl(mt0)/ribtbl(mt0-1)
!
      do iridsign=0,1
      iridstep=(-1)**iridsign 
      do irid= 0,mt*iridstep,iridstep
         do irisign=0,1
         iristep=(-1)**irisign
         do iri= 0,mt*iristep,iristep
! --- Need to pass back the value of B_1 from oursal2 for use here. 
           call oursal2_1a(ribtbl(iri),ridb(irid),slq2b(iri,irid), &
                           smb(iri,irid),shb(iri,irid),ssb(iri,irid), &
                           c_y0,c_y00,iri,irid)
            if(slq2b(iri,irid).lt.0) then
              irimax(irid) = iri - 1 
              go to  15
            endif
         enddo
   15      continue
         enddo
!
      enddo
      enddo
!
! --- Add writes in salinity model case.
!diag   if (mnproc.eq.1) then
!diag   write(lp,*) "************************************************"
!diag   write(lp,*) "New Temperature-Salinity Model"
!diag   write(lp,*) "ifsali=",ifsali
!diag   write(lp,*) "ifsalback=",ifsalback
!
!diag   write(lp,*) "ifepson2=",ifepson2
!diag   if(ifepson2.GT.0) then  &
!diag write(lp,*) "epson2_ref=",epson2_ref
!diag   WRITE(lp,*) "ifdeeplat=",ifdeeplat
!diag   IF(ifdeeplat.GT.0) THEN
!diag   WRITE(*,*) "eplatidepmin=",eplatidepmin
!diag   END IF
!diag   WRITE(*,*) "ifbotenhance=",ifbotenhance
!diag   IF(ifbotenhance.EQ.1) THEN
!diag   WRITE(*,*) "eps_bot0=",eps_bot0
!diag   WRITE(*,*) "scale_bot=",scale_bot
!diag   END IF
!diag   END IF
!*****CD
!
!diag   write(lp,*)"ifrafgmax=",ifrafgmax
!diag   write(lp,*)"ifbg_theta_interp=",ifbg_theta_interp
!diag   write(lp,*)  &
!diag    "    i      ", &
!diag    "    ribtbl(i)      ","    ridb(i)     ", &
!diag    "irimax(i)  "
!diag   do i= -mt,mt
!diag     write(lp,9050) i,ribtbl(i),ridb(i),irimax(i)
!diag   enddo
!
!diag   write(lp,*) " "
!diag   write(lp,*) "irid       Ri_d        Ri(irimax)  " &
!diag 	        // "S_M        S_H        S_S        " &
!diag           // "S_M/S_H    S_S/S_H    "
!diag   do irid= -mt,mt
!diag     write(lp,9100) irid,ridb(irid),ribtbl(irimax(irid)), &
!diag            smb(irimax(irid),irid), &
!diag            shb(irimax(irid),irid), &
!diag            ssb(irimax(irid),irid), &
!diag     	 smb(irimax(irid),irid)/shb(irimax(irid),irid), &
!diag            ssb(irimax(irid),irid)/shb(irimax(irid),irid)
!diag   enddo
!diag   call flush(lp)
!diag   endif !1st proc
!
! --- CALCULATE "R_r_Critical" USING CANUTO'S 000228 ANALYTIC FORMULA
! --- FOR "R_rho_Critical". See NBp.000229-3 and 000316-4.
! --- R_rho_Canuto \equiv -Ri_C/Ri_T \equiv -R_r .
! --- In a sheet dated 000228 Canuto gave me:
! --- "R_\rho^{cr} = {1 \over \Deta} [1 {+\over-} \sqrt{1 - \Delta^2}] 
! --- \Delta \equiv {{\pi_2(1 + {15 \over 7} \pi_3)} \over
! --- {\pi_3 - \pi_2 + (15 \over 14} \pi_3^2}} ".
! --- Note that the + and - choices are reciprocals so this covers
! --- both the Salt Fingering and Double Diffusive Critical R_\rho's.
! --- From Ocean Turbulence III paper have: 
! --- \pi_{1,2,3,4,5} = 
! --- (\tau_pc,\tau_c\theta,\tau_c,\tau_p\theta,\tau_\theta)/\tau 
! --- R_r_Crit = [-1 -/+ \sqrt{1 - \Delta^2}]/Delta
! --- \Delta = {{{\tau_c\theta \over \tau} ( 1 + (15/7)*{\tau_c \over \tau})}
! --- \over {{\tau_c \over \tau} - {\tau_c\theta \over \tau} + 
! --- (15/14) {\tau_c \over \tau}^2}}
!
      deltanum = tctot*(1. + ((15./7.)*tcot))
      deltaden = tcot - tctot + ((15./14.)*(tcot**2))
      delta = deltanum/deltaden
      rrcrn = (-1. - sqrt(1. - (delta**2)))/delta
      rrcrp = (-1. + sqrt(1. - (delta**2)))/delta
      theta_rcrn = atan(rrcrn)
      theta_rcrp = atan(rrcrp)
!
! --- Make sure the right choice of arctan(R_r)=[\theta_r] is made.
! --- Arctan covers the range (-pi/2,pi/2) while 
! --- \theta_r_Crit must be in the range (-pi/4,3pi/4) (The range of Ri>0.)
!
        if(theta_rcrn.lt.-pi/4.) theta_rcrn = theta_rcrn + pi
        if(theta_rcrp.lt.-pi/4.) theta_rcrp = theta_rcrp + pi
      theta_rcrn_deg = theta_rcrn*(180./pi)
      theta_rcrp_deg = theta_rcrp*(180./pi)
!diag   if (mnproc.eq.1) then
!diag   write(lp,*) "   "
!diag   write(lp,*) "   "
!diag   write(lp,*) "   "
!diag   write(lp,*) "   "
!diag   write(lp,*) "R_r_Crit+ =",rrcrp
!diag   write(lp,*) "R_r_Crit- =",rrcrn
!diag   write(lp,*) "\\theta_r_Crit+ =",theta_rcrp
!diag   write(lp,*) "\\theta_r_Crit- =",theta_rcrn
!diag   write(lp,*) "\\theta_r_Crit+ in degrees =",theta_rcrp_deg
!diag   write(lp,*) "\\theta_r_Crit- in degrees =",theta_rcrn_deg
!diag   write(lp,*) "   "
!diag   write(lp,*) "   "
!
!diag   write(lp,*) " "
!diag   write(lp,*) " "
!diag   call flush(lp)
!diag   endif !1st proc
!
! --- Increments in radial and angular coordinates in (Ri_T,Ri_C) plane.
!
        delra_r = 1./real(mt_ra_r)
!       deltheta_r = (pi/4.)/real(n_theta_r_oct)
!
! --- Natassa
!       if (mnproc.eq.1) then
!       write(53,*)nstep,igrid,jgrid,n_theta_r_oct,deltheta_r
!       endif !1st proc
!
! --- Calculate the ratio \sigma_sa_max \equiv S_S/S_H as a function 
! --- of the angle \theta_r in Ri_T,Ri_C space,
! --- \theta_r \equiv arctan(Ri_C/Ri_T). 
! --- The range of angles where unrealizability occurs is 
! --- a subset of theta_r = -pi/4 to 3pi/4.
!
!diag   if (mnproc.eq.1) then
!diag   write(lp,*) "S_S/S_H at pre-maximum Ri as a function of" &
!diag              // "\\theta_r \\equiv Arctan(Ri_C/Ri_T)" 
!
! --- Absurd default on sisamax \equiv S_S/S_H.
!diag   write(lp,*) "Arbitrarily show the absurd value -99.999"
!diag   write(lp,*) "at angles where do not have "// &
!diag   "a maximum Ri (or radius ra_r)."
!diag   write(lp,*) " "
!diag   write(lp,*) "  \\th_r ^o  ra_r      " &
!diag           // "  Ri_T        Ri_C        Ri         Ri_d       " &
!diag           // "  S_M       S_H       S_S      S_S/S_H  "
!diag   call flush(lp)
!diag   endif !1st proc
!
! --- For Ri_T and Ri_C positive find the realizability limits  
! --- in polar coordinates in the (Ri_T,Ri_C) plane : (ra_r,theta_r).
!
        if(ifpolartablewrite.eq. 1 .and. mnproc.eq.1) then
          open(unit=uoff+98,file="turb_ra_th",status="unknown")
        endif
        do itheta_r = -n_theta_r_oct,3*n_theta_r_oct
!       do ihelp = 0,4*n_theta_r_oct
!          itheta_r=ihelp-n_theta_r_oct
         theta_r = real(itheta_r)*deltheta_r
         theta_r_deg = theta_r*(180./pi)
!
! --- Introduce jtheta_r, an angle index that begins at zero   
! --- for the purposes of letting OURSAL2 know it starts at the origin.
!
         jtheta_r = itheta_r + n_theta_r_oct
!
! --- Initialize sisamax to the impossible negative value of -99.999 to 
! --- let places where the realizability limit is not reached stand out.
         sisamax(itheta_r) = -99.999
!
! --- Initialize sm_r0,sh_r0,ss_r0 to the INCONSISTENT absurd value -9.999999.
         sm_r0 = -9.999999
         sh_r0 = -9.999999
         ss_r0 = -9.999999
!
! --- Flag ibg determines if the background value of ra_r has been calculated.
         if(ifsalback.eq. 6) ibg=0
!
! --- Flag ifunreal determines if realizability limit has been found.
         ifunreal=0
!
! --- Make the ra_r max value not too large to try to avoid numerical trouble.
!
         do ira_r = 0,(mt_ra_r**2)/4
!
            if(ira_r.le.mt_ra_r) then
                ra_r = real(ira_r)*delra_r
            else
              ra_r = ((1.+delra_r)**(ira_r - mt_ra_r)) &
                      *(real(mt_ra_r)*delra_r)
            endif
!
! --- Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
            rit = ra_r*COS(theta_r)
            ric = ra_r*SIN(theta_r)
            ri_r  = rit + ric
            rid_r = rit - ric
!
! --- Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
!
            call oursal2_1a(ri_r,rid_r,slq2_r,sm_r,sh_r,ss_r, &
                             c_y0,c_y00,ira_r,jtheta_r)
!
              if(ifpolartablewrite.eq. 1 .and. mnproc.eq.1) then
                write(uoff+98,9001)  &
                itheta_r,theta_r_deg,ira_r,ra_r,slq2_r,sm_r,sh_r,ss_r
              endif
!
! --- Calculate S_M/(S l/q) and find where it's backfact of its origin value.
            if(ifsalback.eq. 6) then
              smosloq_r = sm_r/sqrt(slq2_r)
              if(ira_r.eq. 0) smosloq_0r = smosloq_r
! --- Use radius where dimensionless K_M falls below backfact*origin value.
              if((smosloq_r.le.backfact*smosloq_0r).AND. &
                 (ibg.eq. 0)                            ) then
                ra_r1  = ra_r
                rit1   = rit
                ric1   = ric
                ri_r1  = ri_r
                rid_r1 = rid_r
                slq2_r1(itheta_r) = slq2_r
                sm_r1(itheta_r)   = sm_r
                sh_r1(itheta_r)   = sh_r
                ss_r1(itheta_r)   = ss_r
                ibg=1
              endif
            endif
!
            if(slq2_r.le.0.) then 
! --- Use value of last lattice point on this radius with "slq2" positive.
! --- Calculate the ratio of the salt and heat diffusivities there.
      sisamax(itheta_r) = ss_r0/sh_r0 
!
! --- Store in an array the maximum radius, ra_r, at this angle, theta_r,
! --- in the polar (Ri_T,Ri_C) [that is the (theta_r,ra_r)] plane.
      ra_rmax(itheta_r) = ra_r0
!
! --- Determine the background radius, ra_r, at this \theta_r.
      if(ifsalback.eq. 5) then
! --- Use a constant fraction of the maximum radius before model breakdown.
        back_ra_r(itheta_r) = backfrac*ra_rmax(itheta_r)
!
      else if(ifsalback.eq. 6) then
        back_ra_r(itheta_r) = ra_r1
      endif
!
      ifunreal = 1 
!
! --- Skip straight to write out when last point reached.
      go to 16
            endif
!
            ra_r0   = ra_r
            rit0    = rit
            ric0    = ric
            ri_r0   = ri_r
            rid_r0  = rid_r
            slq2_r0 = slq2_r
            sm_r0   = sm_r
            sh_r0   = sh_r
            ss_r0   = ss_r
!
! --- Store c_y as c_y_0 for possible use as a  guess in background calc.
            c_y_r0(itheta_r) = c_y0
!
         enddo
!
! --- Write out stability functions, the S's and sisamax.
  16    continue
!diag   if (mnproc.eq.1) then
!diag   write(lp,9150) theta_r_deg,ra_r0,rit0,ric0,ri_r0,rid_r0, &
!diag                  sm_r0,sh_r0,ss_r0,sisamax(itheta_r)
!diag   call flush(lp)
!diag   endif !1st proc
!
! --- Set background ra_r large at angles where unrealizability doesn't occur.
! --- Make the ra_r max value not too large to try to avoid numerical trouble.
        if(ifunreal.eq. 0) then
           ipenra_r = (mt_ra_r**2)/4-1
         back_ra_r(itheta_r) = ((1.+delra_r)**(ipenra_r - mt_ra_r)) &
                               *(real(mt_ra_r)*delra_r) 
      endif
!
! --- For ifsalback=5 case get value for initialization of c_y calculation. 
            if(ifsalback.eq. 5) then
              if(jtheta_r.eq. 0) then
                c_y001 = c_y0
              endif
            endif
!
      enddo

!
      if(ifpolartablewrite.eq. 1 .and. mnproc.eq.1) then
        close(uoff+98)
      endif
!
! --- Write out stability functions at background ra_r .
      if(ifsalback.GT.4) then
        do itheta_r = -n_theta_r_oct,3*n_theta_r_oct
           theta_r = real(itheta_r)*deltheta_r
           theta_r_deg = theta_r*(180./pi)
!
! --- Convert radius and angle, (ra_r,theta_r), to rectangular coordinates.
           rit1 = back_ra_r(itheta_r)*COS(theta_r)
           ric1 = back_ra_r(itheta_r)*SIN(theta_r)
           ri_r1  = rit1 + ric1
           rid_r1 = rit1 - ric1
!
! --- Calculation of turbulence functions for ifsalback=5 case.
           if(ifsalback.eq. 5) then
!
! --- Calculate turbulence functions at this radius and angle in (Ri_T,Ri_C).
             jtheta_r = itheta_r + n_theta_r_oct
!
! --- Set second table index to 1 to use last step's value except at start.
! --- Transform that "last step" value from the most recent angle step to the
! --- final realizable ra_r step at {\it this} angle in hope of more accuracy.
             call oursal2_1a(ri_r1,rid_r1,slq2_r1(itheta_r), &
                       sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
                              c_y_r0(itheta_r),c_y001,jtheta_r,1)
           endif
!
!diag        if(itheta_r.eq. -n_theta_r_oct) then
!diag          if (mnproc.eq.1) then
!diag          write(lp,*) " "
!diag          write(lp,*)  &
!diag           "Values at background ra_r=(Ri_T^2 + Ri_C^2)^(1/2)"
!diag          write(lp,*) "\\th_r ^o   ra_r       " &
!diag            // "Ri_T       Ri_C       Ri         Ri_d       " &
!diag            // "(Sl/q)^2   S_M       S_H       S_S       S_S/S_H  "
!diag          write(lp,*) " "
!diag          call flush(lp)
!diag          endif !1st proc
!diag        endif
!
           sisa1 = ss_r1(itheta_r)/sh_r1(itheta_r)
!
!          if (mnproc.eq.1) then
!          write(lp,*)
!    &       'itheta_r,theta_r_deg = ',itheta_r,theta_r_deg
!          write(lp,*)
!    &       'back_ra_r,slq2_r1    = ',
!    &        back_ra_r(itheta_r),slq2_r1(itheta_r)
!          write(lp,*)
!    &       'sm_r1,sh_r1,ss_r1    = ',
!    &        sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r)
!          call flush(lp)
!          endif !1st proc
!diag        if (mnproc.eq.1) then &
!diag      write(lp,9160) theta_r_deg,back_ra_r(itheta_r), &
!diag                      rit1,ric1,ri_r1,rid_r1,slq2_r1(itheta_r), &
!diag                   sm_r1(itheta_r),sh_r1(itheta_r),ss_r1(itheta_r), &
!diag                   sisa1
!diag        call flush(lp)
!diag        endif !1st proc
!
           if(slq2_r1(itheta_r).lt.0.) then
           if (mnproc.eq.1) then
           write(lp,*)  &
              "Negative (Sl/q)^2 in table of Background vs. \\theta_r."
           write(lp,*) "itheta_r=",itheta_r, &
      	                 "   slq2_r1(itheta_r)=",slq2_r1(itheta_r)
           write(lp,*) "Program is stopping in turb_2."
           endif !1st proc
           call xcstop('(inigiss)')
                  stop '(inigiss)'
           endif
        enddo
      endif
!
!
 9001 format(2(I8,'   ',1pe11.3),8(1pe11.3))
 9050 format(I8,'  ',2E16.4,I8,'  ')
 9100 format('  ',I8,'  ',2E12.4,3F11.6,2F11.4)
 9150 format(F11.3,5E12.4,3F10.6,F9.3)
 9160 format(F11.3,1x,6(E10.4,1x),3(F10.6,1x),F9.3)
 9200 format(I12,'    ',5E16.6)

!
      return
      end
!
      subroutine oursal2_1a(ri,rid,slq2,sm,sh,sc,c_y0,c_y00,iri,irid)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none
!
! --- Replace the numerical value of 6.25 by 1/(tpvot**2) .
! --- Version in which following OTsalche/plot000127 
! --- the timescale ratios are calculated in the 'smshsc' routine 
! --- and passed back hrough the common block /bb/
! --- to simplify the process of adjustment of timescale ratios.
! --- Submodule to calculate turbulence functions (Sl/q)^2 and S_M,S_H,S_S
! --- of Ri(=Ri_T+Ri_C) and Ri_d(=Ri_T-Ri_C) in our NCAR turbulence module.
! --- Stripped and adapted from plot981007.f.
! --- Program to generate contour and 1 variable plots vs. Ri,Ri_d based on
! --- Program to generate contour plots vs. Ri_T and Ri_C based on
! --- Program to generate plots vs. Ri_T at different Ri_C values based on
! --- .or.eC.eD PROGRAM WITH .eW VAL.e OF "p10". 'p10 = tpt*tct/(tc**2)'
! --- Program to generate K_X/((l^2) S) for Canuto based on plot980609.f:
! --- Program to generate data for plots of turbulence functions including
! --- S_{M,H,C} and Canuto's new y = (\tau_pv S)^2
! --- and n,c as functions of stability parameters in the concentration theory
! --- (structure is a 1 point closure like the generalized Mellor-Yamada, 
! --- but the constants are derived based on Dubovikov's model according
! --- to Ye Cheng). The concentration theory dimensionless parameters
! --- associated with the squares of shear, temperature contribution to 
! --- Brunt Vaisala frequency and concentration contribution to it,
! --- the new y,n,c are represented in this program by the variables
! --- c_y,c_n,c_c. 
! --- Adapted from Cheng's program mike_12.f_980528 for the Dubovikov model.

!-----------------------------------------------------------------------
!
! --- y=(tau*s)**2
! --- tau=2*e/epsilon=b1*l/q
! --- km=e*tau*sm=1/2*(b1*l)**2*s/y**(1/2)*sm
! --- kh=e*tau*sh=1/2*(b1*l)**2*s/y**(1/2)*sh
! --- ks=e*tau*ss=1/2*(b1*l)**2*s/y**(1/2)*ss
!
! --- X = {M,H,C} . 
! --- Cheng above gives K_X = (1/2)((B_1*l)^2) (S/(((\tau S)**2)^(1/2))) S_X
! --- The "old" y used above is (\tau S)^2. 
! --- The "new" y (c_y in the program) is (\tau_pv S)^2.
! --- The program variable "slq2" is (S l/q)^2 = y (B_1)^(-2), 
! --- since \tau=B_1 l/q. (S l/q)^2 = (\tau \over \tau_pv)^2 c_y (B_1)^(-2) .
! --- c_y = (S l/q)^2 * [(B_1)^2 * (\tau_pv \over \tau)^2] .
!
! --- Take \tau_pv/\tau as being calculated in the smshsc routine instead.
! --- From the printed notes Canuto gave me on 980601 have:
! ---          \tau_pv = {2 \over 5} \tau  (B.1) or parameter(tpvot = 0.4)
!
      real ri,rid,slq2,sm,sh,sc,c_y0,c_y00
      real eeps,c_yst,c_yst0,c_y,val,c_n,c_c
!DBI: all eps ==> eeps in this routine!
      integer iri,irid,iend,ier
      real       rit,ric
      common /bb/rit,ric       !rit is the temperature's part of 
                               !Ri and ric the concentration's.
      save   /bb/
!
      parameter(c_yst0 = 8.527882) !Need a guess for c_y for the solver 
                                   !for the neutral case, c_yst. Take 
                                   !c_yst = 8.527882, the approximate value 
                                   !calculated at rit=ric=0. A variable c_y00 
                                   !is intended to hold the Ri=0 value of c_y 
                                   !from the previous Ri_d row in a table the 
                                   !subroutine is being called to make and a 
                                   !variable c_y0 is intended to hold the
                                   !previous Ri value from the current Ri_d
                                   !row of that table.
        b1=16.6
!
! --- Commented excerpt from the file "sx"
!
! --- sgmt := 0.72;
!
! --- tpt := 1/(5*(1+1/sgmt))*tau;
! --- tpt  = .08372093019*tau
!
! --- tpc := 1/(5*(1+1/sgmt))*tau;
! --- tpc  = .08372093019*tau
!
! --- tt := sgmt*tau;
! --- tt   = .72*tau
!
! --- tc := sgmt*tau;
! --- tc   = .72*tau
!
! --- tct := 2/15*sgmt*tau;
! --- tct  = .09599999998*tau
!
! --- Calculate the timescale ratios in the 'smshsc' routine instead of here.
! --- Set \sigma_t0. sgmt = .72
!
! --- Calculate {\tau_C \over \tau} and {\tau_{C\theta} \over \tau}.
! --- tcot  = sgmt
! --- tctot = (2./15.)*sgmt
! --- "tpt/tau" and "tpc/tau" from the "sx" excerpt
! --- tptot = 1./(5.*(1+1/sgmt))
! --- tpcot = 1./(5.*(1+1/sgmt))
!
! --- Timescale ratios are now calculated in the 'smshsc' subroutine.
! --- Make dummy call with c_y=c_n=c_c=0 to get their values for initial use.
      call smshsc_a3(0.,0.,0.,sm,sh,sc)
!
      eeps=1.e-6                                              
      iend=300                                              
!
! --- rimax= ?
! --- rtwi finds the root of x=fct_sal(x)                     
! --- Need a guess at the root, c_yst. Use neighboring solution.
! --- Initial guess for c_yst for this value of Ri_d.
      if(iri.eq.0.and.irid.eq.0) then
      c_yst = c_yst0
      else if(iri.eq.0) then
      c_yst = c_y00
      else 
      c_yst = c_y0
      endif
!
! --- Calculate Ri_T =(Ri + Ri_d)/2 and Ri_C =(Ri - Ri_d)/2.  
       rit = (ri + rid)/2.
       ric = (ri - rid)/2.
         call rtwi(c_y,val,c_yst,eeps,sm,sh,sc,iend,ier)
!
         if(ier.ne.0) then
! --- Make error message more specific.
          if (mnproc.eq.1) then
          write(lp,*) "In oursal2 subroutine"
          write(lp,*) "c_y00=",c_y00,"	c_y0=",c_y0
          write(lp,*) "ri=",ri,"	rid=",rid
          write(lp,*) "rit=",rit,"	ric=",ric
          write(lp,*) "Initial guess for rtwi c_yst=",c_yst
!
          write(lp,*) "rtwi call problem, ier=",ier
          endif !1st proc
          call xcstop('(oursal2_1a)')
                 stop '(oursal2_1a)'
         endif
!
! --- Calculate (S l/q)^2[=program variable "slq2"] from c_y.**
! --- (S l/q)^2 = (\tau \over \tau_pv)^2 c_y (B_1)^(-2) .
! --- (S l/q)^2 = (\tau_pv \over \tau)^(-2) c_y (B_1)^(-2) .
       slq2 = c_y/((b1*tpvot)**2)
!
! --- Store value of c_y for future guesses.
       if(c_y.ge.0) then
          c_y0=c_y
       else 
! --- Turbulence model becomes unphysical for c_y negative.
! --- Realizability for negative Ri
           if(ri.lt.0) then
           if (mnproc.eq.1) then
           write(lp,*) "c_y negative at negative Ri"
           write(lp,*) "Ri=",ri," 	c_y=",c_y
           write(lp,*) "Unstable realizability limit unexpected:" 
           write(lp,*) "stopping in oursal2."
           endif !1st proc
           call xcstop('(oursal2_1a)')
                  stop '(oursal2_1a)'
           endif
       endif
!
       if(iri.eq.0) c_y00=c_y
       if((iri.eq.0).and.(irid.eq.0).and. &
           (abs(c_y - c_yst0).gt.1.e-6)) then
         if (mnproc.eq.1) then
         write(lp,*) "Inconsistency in neutral value of c_y"
         write(lp,*) "Value used =",c_yst0
         write(lp,*) "Value calculated =",c_y
         write(lp,*) "Program stopping in oursal2"
         endif !1st proc
         call xcstop('(oursal2_1a)')
                stop '(oursal2_1a)'
       endif
!
! --- From last page (#5) of "980608 AH Concentration Work" handwritten
! --- sheetsC have: 
! --- n = -{{\tau_C \tau_{C\theta}} \over {\tau_{pv}}^2 } y Ri_T
! --- c = - {{\tau_C}^2 \over {\tau_{pv}}^2} y Ri_C
! --- Decide to use the parameter "tpvot" instead of its value 2/5 \tau .
! --- n = -{{(\tau_C/\tau) (\tau_{C\theta}/\tau)} \over {\tau_{pv}/\tau}^2 }
! --- y Ri_T
! --- c = - {{\tau_C/\tau}^2 \over {\tau_{pv}/\tau}^2} y Ri_C
!
         c_n = -(tcot*tctot/(tpvot**2))*c_y*rit
         c_c = -((tcot**2)/(tpvot**2))*c_y*ric
         call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)
!
!
 1003 format(12(I8))
 1004 format(12(1pe14.5))
      end
!-----------------------------------------------------------------------
      function fct_sal(sm,sh,sc,c_y)                              
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none  
!
      real fct_sal,c_n,c_c,c_y,sm,sh,sc
!
      real       rit,ric
      common /bb/rit,ric
      save   /bb/
!
! --- Decide to use the parameter "tpvot" instead of its value 2/5 \tau .
      c_n = -((tcot*tctot)/(tpvot**2))*c_y*rit
      c_c = -((tcot**2)/(tpvot**2))*c_y*ric
      call smshsc_a3(c_y,c_n,c_c,sm,sh,sc)
!
! --- y(S_\nu - Ri_T S_h - Ri_C S_c) = 8/25 . 8/25 = 0.32 . S_\nu = sm.
! --- y = 0.32/(S_\nu - Ri_T S_h - Ri_C S_c). 
      fct_sal=(2.*(tpvot**2))/(sm-rit*sh-ric*sc)
      return                                          
      end                                            
!-----------------------------------------------------------------------
      subroutine smshsc_a3(yyy,nnn,ccc,sm,sh,sc)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none   
!
! --- .eW SUBROUTI.e WHICH calculates the "p's" from the timescale ratios.
! --- BA.eD on "smshsc2":
! --- SUBROUTI.e WHICH CALCULA.eS "p's" from "sgmt". BA.eD ON "smshsc1":
! --- .eW SUBROUTI.e WHICH U.eS .e C.eNG'S .orTRAN CO.e TO CALCULA.e CONSTANTS
! --- FROM T.e "p's" .eNT TO .e BY HIM TODAY. BA.eD ON "smshsc0". 
! --- **.or.eCT T.e VAL.e OF "p10".**
! --- p_10 = {\tau_{p \theta} \tau_{c \theta}} \over {\tau_c ^ 2}
!
! --- Replace Cheng's smsh with  smshsc, which includes concentration.
! --- The y,n,c used here are Canuto's "y,n,c" called c_y,c_n,c_c 
! --- elsewhere in this program.
      real   yyy,nnn,ccc,sm,sh,sc
      real   Nm,Nh,Nc
!
      real p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p1m,p2m
      real a0,a1,a2,a3,a4,a5
      real d0,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15
      real D
!
      integer ifrecall,ifmodelconstout
!
! --- Switch for whether(1) or not(0) to output p's a's and d's to a file.
      parameter(ifmodelconstout=0)
!
! --- Add `\tau_pv \over \tau' to the common block with timescale ratios.
!
! --- Calculate the p's.
      p1  = 0.832
      p2  = 0.545
      p3  = (5./2.)*tpcot
      p4  = (1./5.)*tpcot*(tcot**(-2))
      p5  = tpcot*tctot*(tcot**(-2))
      p6  = (1./5.)*(tcot**(-1))*(tctot**(-1))*tptot
      p7  = 5.*tctot
      p8  = (5./2.)*tptot
      p9  = ttot*tptot*((tcot*tctot)**(-1))
      p10 = tctot*tptot*(tcot**(-2))
      p11 = tpcot*(tcot**(-1))
      p1m = 1. - p1
      p2m = 1. - p2
!
!-----------------------------------------------------------------------
!results.2_1
! --- Values of a's and d's calculated from p's using Cheng's Fortran code
! --- to do so, from today's email from him, cheng990513.results.2_1 .
! --- results.2_1
!##########################
!##  Fortran code:
!##########################
      A0 = 12
      A1 = p11*(12*p9+8*p6-30*p6*p8-5*p6*(p1m+3*p2m))
      A2 = 5*(2*p4*p6*p7-p4*p9-p6*p11)*(p1m+3*p2m)+8*p6*p11+8*p4*p9- &
           16*p4*p6*p7+12*p11*p9+12*p11*p10-12*p4*p7**2*p6- &
           30*p6*p11*p8+30*p4*p6*p7*p8+30*p6*p4*p7*p3-30*p4*p9*p3
      A3 = p10*(12*p11+8*p4-30*p3*p4-5*p4*(p1m+3*p2m))
      A4 = -p6*(8-30*p8-5*p1m-15*p2m)-12*p9-12*p11
      A5 = -p4*(8-30*p3-5*p1m-15*p2m)-12*p10-12*p11
      D0 = 24
      D1 = p11*((-p6-2*p9)*p1m**2+(p6+6*p9)*p2m**2+2*p6*p8*(p1m-3*p2m))
      D2 = (2*p4*p6*p7-p4*p9-p6*p11)*(p1m**2-p2m**2)+ &
         2*(-p11*p10-p11*p9+p4*p7**2*p6)*(p1m**2-3*p2m**2)+ &
         2*(-p6*p4*p7*p3-p4*p6*p7*p8+p4*p9*p3+p6*p11*p8)*(p1m-3*p2m)
      D3 = p10*((-p4-2*p11)*p1m**2+2*p4*p3*(p1m-3*p2m)+(6*p11+p4)*p2m**2)
      D4 = -4*p6*p11*(3*p9+2*p6)
      D5 = 4*p4*p6**2*p7*(4+3*p7)-4*p4*p9*(3*p11+2*p6)- &
           4*p6*p11*(3*p9+3*p10+2*p4+2*p6)
      D6 = 4*p4**2*p6*p7*(4+3*p7)-4*p4*p9*(2*p4+3*p11)- &
           8*p4*p6*(p11+p10)-12*p10*p11*(p4+p6)
      D7 = -4*p4*p10*(2*p4+3*p11)
      D8 = (2*p9+2*p11+p6)*p1m**2-2*p6*p8*(p1m-3*p2m)- &
           (p6+6*p9+6*p11)*p2m**2
      D9 = (2*p10+p4+2*p11)*p1m**2-2*p4*p3*(p1m-3*p2m)- &
           (p4+6*p10+6*p11)*p2m**2
      D10 = 8*p6**2+4*(7*p11+3*p9)*p6+24*p11*p9
      D11 = -8*(4+3*p7)*p4*p6*p7+4*p4*(4*p6+7*p9+3*p11)+ &
             4*p6*(3*p10+7*p11)+24*p11*(p10+p9)
      D12 = 4*p10*(7*p4+6*p11)+4*p4*(2*p4+3*p11)
      D13 = 6*p2m**2-2*p1m**2
      D14 = -28*p6-24*p9-24*p11
      D15 = -24*p10-28*p4-24*p11
!results.2_1
!-----------------------------------------------------------------------
!
! --- Write out the p's.
! --- Writeout the timescale ratios as well.
        ifrecall=1
      if(ifrecall.eq.0 .and. mnproc.eq.1) then
        write(lp,*) "tau_pv/tau     =",tpvot 
        write(lp,*) "tau_ptheta/tau =",tptot
        write(lp,*) "tau_pc/tau =",tpcot
        write(lp,*) "tau_theta/tau  =",ttot
        write(lp,*) "tau_c/tau  =",tcot
        write(lp,*) "tau_ctheta/tau  =",tctot
        write(lp,*) " "
        write(lp,*) "p1 =",p1
        write(lp,*) "p2 =",p2
        write(lp,*) "p3 =",p3
        write(lp,*) "p4 =",p4
        write(lp,*) "p5 =",p5
        write(lp,*) "p6 =",p6
        write(lp,*) "p7 =",p7
        write(lp,*) "p8 =",p8
        write(lp,*) "p9 =",p9
        write(lp,*) "p10=",p10
        write(lp,*) "p11=",p11
!
! --- Write out the a's and d's as well.
        write(lp,*) "a0=",a0
        write(lp,*) "a1=",a1
        write(lp,*) "a2=",a2
        write(lp,*) "a3=",a3
        write(lp,*) "a4=",a4
        write(lp,*) "a5=",a5
        write(lp,*) "d0=",d0
        write(lp,*) "d1=",d1
        write(lp,*) "d2=",d2
        write(lp,*) "d3=",d3
        write(lp,*) "d4=",d4
        write(lp,*) "d5=",d5
        write(lp,*) "d6=",d6
        write(lp,*) "d7=",d7
        write(lp,*) "d8=",d8
        write(lp,*) "d9=",d9
        write(lp,*) "d10=",d10
        write(lp,*) "d11=",d11
        write(lp,*) "d12=",d12
        write(lp,*) "d13=",d13
        write(lp,*) "d14=",d14
        write(lp,*) "d15=",d15
!
! --- Output p#, a# and d# to the file model_constants if the switch is set.
! --- Writeout the timescale ratios as well.
          if(ifmodelconstout.eq.1 .and. mnproc.eq.1) then
            open(unit=uoff+98,file='model_constants',status='unknown')
            write(uoff+98,*) "tau_pv/tau     =",tpvot 
            write(uoff+98,*) "tau_ptheta/tau =",tptot
            write(uoff+98,*) "tau_pc/tau =",tpcot
            write(uoff+98,*) "tau_theta/tau  =",ttot
            write(uoff+98,*) "tau_c/tau  =",tcot
            write(uoff+98,*) "tau_ctheta/tau  =",tctot
            write(uoff+98,*) " "
            write(uoff+98,*) "p1 =",p1
            write(uoff+98,*) "p2 =",p2
            write(uoff+98,*) "p3 =",p3
            write(uoff+98,*) "p4 =",p4
            write(uoff+98,*) "p5 =",p5
            write(uoff+98,*) "p6 =",p6
            write(uoff+98,*) "p7 =",p7
            write(uoff+98,*) "p8 =",p8
            write(uoff+98,*) "p9 =",p9
            write(uoff+98,*) "p10=",p10
            write(uoff+98,*) "p11=",p11
            write(uoff+98,*) "a0 =",a0
            write(uoff+98,*) "a1 =",a1
            write(uoff+98,*) "a2 =",a2
            write(uoff+98,*) "a3 =",a3
            write(uoff+98,*) "a4 =",a4
            write(uoff+98,*) "a5 =",a5
            write(uoff+98,*) "d0 =",d0
            write(uoff+98,*) "d1 =",d1
            write(uoff+98,*) "d2 =",d2
            write(uoff+98,*) "d3 =",d3
            write(uoff+98,*) "d4 =",d4
            write(uoff+98,*) "d5 =",d5
            write(uoff+98,*) "d6 =",d6
            write(uoff+98,*) "d7 =",d7
            write(uoff+98,*) "d8 =",d8
            write(uoff+98,*) "d9 =",d9
            write(uoff+98,*) "d10=",d10
            write(uoff+98,*) "d11=",d11
            write(uoff+98,*) "d12=",d12
            write(uoff+98,*) "d13=",d13
            write(uoff+98,*) "d14=",d14
            write(uoff+98,*) "d15=",d15
            close(uoff+98)
          endif
!
      endif
      ifrecall = 1
!

! --- Modification of section of "sx" containing the den and nums of the "S"'s

!###############################################

         D = d0 + d1*yyy*nnn**2 + d2*yyy*nnn*ccc + d3*yyy*ccc**2 &
            + d4*nnn**3 + d5*nnn**2*ccc + d6*nnn*ccc**2 + d7*ccc**3 &
            + d8*yyy*nnn + d9*yyy*ccc + d10*nnn**2 + d11*nnn*ccc &
            + d12*ccc**2 + d13*yyy &
            + d14*nnn + d15*ccc
!########################################################################

         Nm = a0 + a1*nnn**2 + a2*nnn*ccc + a3*ccc**2 + a4*nnn + a5*ccc

!###########################################################################
         Nh  = - (30.*nnn*p6 + 30.*ccc*p4 - 60. &
               - ( 2.*p1m  + 15.*p2m**2 - 6.*p2m  - 5.*p1m**2 ) &
               * yyy ) &
               * (ccc*p4*p7 - ccc*p11 - nnn*p11 + 1.)


         Nc  =   (30.*nnn*p6 + 30.*ccc*p4 - 60. &
               - ( 2.*p1m  + 15.*p2m**2 - 6.*p2m  - 5.*p1m**2 ) &
               * yyy ) &
               * (ccc*p10 - 1. - nnn*p6*p7 + nnn*p9)
! --- Modification of section of "sx" containing Sm, Sh, Sc
!*******************************************************************************
         sm = (4./15.) * tpvot * Nm/D

         sh = (4./15.) * tptot * Nh/D

         sc = (4./15.) * tpcot * Nc/D
!*******************************************************************************
      return
 1004 format(12(1pe14.5))
      end
!-----------------------------------------------------------------------
      subroutine rtwi(xx,val,xst,eeps,sm,sh,sc,iend,ier)                      
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none
!
      real fct_sal,xx,val,xst,eeps,sm,sh,sc
      real tol,a,b,d
      integer iend,ier,i
!
      real       rit,ric
      common /bb/rit,ric
      save   /bb/
!
! --- to solve general nonlinear equations of the form x=fct_sal(x)       
! --- by means of wegsteins iteration method                         
! --- prepare iteration                                             
!
!
      ier=0                                                        
      tol=xst                                                     
      xx=fct_sal(sm,sh,sc,tol)                                                 
      a=xx-xst                                                   
      b=-a                                                     
      tol=xx                                                   
      val=xx-fct_sal(sm,sh,sc,tol)
!
! --- start iteration loop                                 
      do 6 i=1,iend                                       
!
! --- Crude fix to avoid mysterious problem which occurred 
! --- with a close but not too close guess.
      if(abs(val).lt.1.e-12) val =0.
!MHRI: The arithmetic IF statement is a deleted feature from the Fortran 2015 standard.
!MHRI      if(val) 1,7,1                                      
      if (val==0.) then
        goto 7
      else
        goto 1
      endif
!
! --- equation is not satisfied by x                    
 1    b=b/val-1.                                       
!MHRI      if(b) 2,8,2                                     
      if (b==0.) then
        goto 8
      else
        goto 2
      endif
!
! --- iteration is possible                          
 2    a=a/b                                         
      xx=xx+a                                        
      b=val                                       
      tol=xx                                      
      val=xx-fct_sal(sm,sh,sc,tol) 
!
! --- test on satisfactory accuracy            
      tol=eeps                                 
      d=abs(xx)                               
!MHRI      if(d-1.) 4,4,3                        
      if ((d-1.)<=0.) then
        goto 4
      else
        goto 3
      endif
 3    tol=tol*d                            
!MHRI 4    if(abs(a)-tol) 5,5,6                
 4    if ((abs(a)-tol)<=0.) then
        goto 5                
      else
        goto 6
      endif
!MHRI 5    if(abs(val)-10.*tol) 7,7,6         
 5    if((abs(val)-10.*tol)<=0.) then
        goto 7
      else
        goto 6
      endif
 6    continue                          
!
! --- end of iteration loop                                           
! --- no convergence after iend iteration steps. error return.       
      ier=1                                          
 7    return                                        
!
! --- error return in case of zero divisor         
 8    ier=2                                       
!
      return                                     
      end                                       
!
      subroutine interp2d_expabs(ri,rid,slq2,sm,sh,ss,m,m0,delta,rat)
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none
!
! --- Subroutine for a modular interpolation calculation.
! --- provides a faster interpolation calculation in the ifexpabstable=1 case.
!
! --- Output interpolated values: slq2,sm,sh,ss
!
        real ri,rid,slq2,sm,sh,ss,delta,rat
        real deltaridta,deltarita,deltarid,deltari,dslq2_rid,dslq2_ri, &
             dsm_rid,dsm_ri,dsh_rid,dsh_ri,dss_rid,dss_ri, &
             tabindrid,tabindri
        integer m,m0
        integer lrid0,lrid1,lri0,lri1,l
!
! --- Take values off the edge of the table back to the table on radial lines.
! --- Must use ratio of Ri's taken *before* the cut-off has taken place.
      if(ri.gt.ribtbl(m)) then
        if(abs(rid).le.ri) then
          rid = ribtbl(m)*(rid/ri)
          ri  = ribtbl(m)
        else if(rid.gt.ri) then
          ri  = ridb(m)*(ri/rid)
          rid = ridb(m)
        else if(rid.lt.-ri) then
          ri  = ridb(-m)*(ri/rid)
          rid = ridb(-m)
        end if
      else if(ri.lt.ribtbl(-m)) then
        if(abs(rid).le.-ri) then
          rid = ribtbl(-m)*(rid/ri)
          ri  = ribtbl(-m)
        else if(rid.gt.-ri) then
          ri  = ridb(m)*(ri/rid)
          rid = ridb(m)
        else if(rid.lt.ri) then
          ri  = ridb(-m)*(ri/rid)
          rid = ridb(-m)
        end if
      else if(rid.gt.ridb(m)) then
          ri  = ridb(m)*(ri/rid)
          rid = ridb(m)
      else if(rid.lt.ridb(-m)) then
          ri  = ridb(-m)*(ri/rid)
          rid = ridb(-m)
      end if
!
! --- Interpolate points within the table range.
! --- Table index ranges from -m to m with equal spacing for -m0 to m0.
      if(abs(rid).lt.ridb(m0)) then
!
! --- Find Interpolation points in the equally spaced Ri_d part of the table.
        lrid1 = int(rid/delta)+nint(sign(real(1),rid))
! --- Find Interpolation points in exponential absolute value spaced Ri_d
! --- part of the table.
      else if((abs(rid)).ge.(ridb(m))) then
! --- Special case where have a value which falls at the limit of the table.
        lrid0 = nint(sign(real(m),rid))
        lrid1 = lrid0
        GO TO 252
!
      else
        tabindrid = sign( &
           real(m0) + ((log(abs(rid)) - log(ridb(m0)))/log(rat)), &
           rid)
        lrid1 = int(tabindrid)+nint(sign(real(1),rid))
!
      endif
!
! --- It is conceivable that rounding errors may in borderline cases 
! --- throw the calculated table indices for Ri_d off by one.
! --- Check and allow moving one to either side to take care of this.
      if ((abs(ridb(lrid1))).lt.(abs(rid))) then
        lrid1 = lrid1 + sign(1,lrid1)
      else if ((abs(ridb(lrid1-sign(1,lrid1)))).gt.(abs(rid))) then
        lrid1 = lrid1 - sign(1,lrid1)
      end if
!
  250 continue
!
! --- Make lrid0 one less or greater than lrid1 according to sgn(rid).
        lrid0 = lrid1 - nint(sign(real(1),rid)) 
!
        if(rid.eq.0.0) lrid1 = 1
  252   continue
!
! --- Check that the Ri_d value falls within the interpolation interval.
      if(( rid.gt.0.0.and. &
          (rid.lt.ridb(lrid0).or.rid.gt.ridb(lrid1))).or. &
         ( rid.lt.0.0.and. &
          (rid.gt.ridb(lrid0).or.rid.lt.ridb(lrid1)))    ) then
       if (mnproc.eq.1) then
       WRITE(lp,*) "Ri_d is outside interpolation range in interp2d_e."
       WRITE(lp,*) "rid=  ",rid,"lrid0= ",lrid0,"lrid1= ",lrid1
       WRITE(lp,*) "ridb(lrid0)=  ",ridb(lrid0), &
                      "   ridb(lrid1)= ",ridb(lrid1)
       WRITE(lp,*) "Program is stopping."
       endif !1st proc
       call xcstop('(interp2d_expabs)')
              stop '(interp2d_expabs)'
      end if
!
! --- Artificially reduce Ri if it threatens to surpass Ri_max(Ri_d).
! --- This is to conform to the 1D table's realizability limit treatment. 
! --- if(ri.gt.MIN(ribtbl(irimax(lrid0)),ribtbl(irimax(lrid1)))) then
! --- ri = MIN(ribtbl(irimax(lrid0)),ribtbl(irimax(lrid1)))
! --- end if
!
! --- Set turbulence to zero if Ri threatens to surpass the realizability limit.
        if(ri.gt.MIN(ribtbl(irimax(lrid0)),ribtbl(irimax(lrid1)))) then
          slq2=0.0
          sm = 0.0
          sh = 0.0
          ss = 0.0
          return
        end if
!
! --- Table index ranges from -m to m with equal spacing for -m0 to m0.
      if(abs(ri).lt.ribtbl(m0)) then
!
! --- Find Interpolation points in the equally spaced Ri part of the table.
        lri1 = int(ri/delta)+nint(sign(real(1),ri))
!
! --- Find Interpolation points in exponential absolute value spaced Ri
! --- part of the table.
      else if((abs(ri)).ge.(ribtbl(m)))  &
         then
!
! --- Special case where have a value which falls at the limit of the table.
        lri0 = nint(sign(real(m),ri))
        lri1 = lri0
        GO TO 272
!
      else
        tabindri = sign( &
           real(m0) + ((log(abs(ri)) - log(ribtbl(m0)))/log(rat)), &
           ri)
        lri1 = int(tabindri)+nint(sign(real(1),ri))
!
  270 continue
      end if
!
! --- It is conceivable that rounding errors will in borderline cases 
! --- throw the calculated table indices for Ri off by one.
! --- Check and allow moving one to either side to take care of this.
      if((abs(ribtbl(lri1))).lt.(abs(ri))) then
        lri1 = lri1 + sign(1,lri1)
      else if((abs(ribtbl(lri1-sign(1,lri1)))).gt.(abs(ri))) then
        lri1 = lri1 - sign(1,lri1)
      end if
!
! --- Make lri0 one less or greater than lri1 according to sgn(ri).
        lri0 = lri1 - nint(sign(real(1),ri)) 
!
        if(ri.eq.0.0) lri1 = 1
  272 continue
!
! --- check that the Ri_d value falls within the interpolation interval.
      if((ri.gt.0.0.and.(ri.lt.ribtbl(lri0).or.ri.gt.ribtbl(lri1))) &
          .or.(ri.lt.0.0.and.(ri.gt.ribtbl(lri0) &
          .or.ri.lt.ribtbl(lri1)))) then
       if (mnproc.eq.1) then
       WRITE(lp,*) "Ri is outside interpolation range in interp2d_e."
       WRITE(lp,*) "ri=  ",ri,"lri0= ",lri0,"lri1= ",lri1
       WRITE(lp,*) "ribtbl(lri0)=  ",ribtbl(lri0), &
                      "   ribtbl(lri1)= ",ribtbl(lri1)
       WRITE(lp,*) "Program is stopping."
       endif !1st proc
       call xcstop('(interp2d_expabs)')
              stop '(interp2d_expabs)'
      end if
!
! --- interpolate turbulence fields and introduce table spacing variables.
      deltaridta = ridb(lrid1) - ridb(lrid0)
      deltarita  = ribtbl(lri1)  - ribtbl(lri0)
      deltarid = rid - ridb(lrid0)
      deltari  = ri - ribtbl(lri0)
!
! --- set delta field to zero in special cases falling at limit of the table. 
      if(lrid1.eq.lrid0) then
        dslq2_rid = 0.0
      else
        dslq2_rid = (slq2b(lri0,lrid1) - slq2b(lri0,lrid0))/ &
                      deltaridta
      end if
      if(lri1.eq.lri0) then
        dslq2_ri  = 0.0
      else
        dslq2_ri = (slq2b(lri1,lrid0) - slq2b(lri0,lrid0))/ &
                     deltarita
      end if
      slq2 = slq2b(lri0,lrid0) + dslq2_ri*deltari + dslq2_rid*deltarid
!
! --- sm
      if(lrid1.eq.lrid0) then
        dsm_rid   = 0.0
      else
        dsm_rid = (smb(lri0,lrid1) - smb(lri0,lrid0))/ &
                    deltaridta
      end if
      if(lri1.eq.lri0) then
        dsm_ri    = 0.0
      else
        dsm_ri = (smb(lri1,lrid0) - smb(lri0,lrid0))/ &
                   deltarita
      end if
      sm     = smb(lri0,lrid0) +  &
                  dsm_ri*deltari + dsm_rid*deltarid
!
! --- sh
      if(lrid1.eq.lrid0) then
        dsh_rid   = 0.0
      else
        dsh_rid = (shb(lri0,lrid1) - shb(lri0,lrid0))/ &
                    deltaridta
      end if
      if(lri1.eq.lri0) then
        dsh_ri    = 0.0
      else
        dsh_ri = (shb(lri1,lrid0) - shb(lri0,lrid0))/ &
                    deltarita
      end if
      sh     = shb(lri0,lrid0) +  &
                  dsh_ri*deltari + dsh_rid*deltarid
!
! --- ss
      if(lrid1.eq.lrid0) then
        dss_rid   = 0.0
      else
        dss_rid = (ssb(lri0,lrid1) - ssb(lri0,lrid0))/ &
                    deltaridta
      end if
      if(lri1.eq.lri0) then
        dss_ri    = 0.0
      else
        dss_ri = (ssb(lri1,lrid0) - ssb(lri0,lrid0))/ &
                   deltarita
      end if
!
      ss     = ssb(lri0,lrid0) +  &
                  dss_ri*deltari + dss_rid*deltarid
!
      return
      end
!
!> Revision history:
!>
