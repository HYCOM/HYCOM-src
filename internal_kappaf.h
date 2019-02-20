!-----------------------------------------------------------------------------
      real function kappaf1(t,s,r,prs,kkf)
      implicit none
!
      real    t,s,r,prs
      integer kkf
!
! --- kappaf1 used to simplify offsetting T and S,
! --- always invoke via kappaf.
!
! --- coefficients for kappa^(theta)
! --- 1=Arctic/Antarctic; 2=Atlantic; 3=Mediterranean
!
      real, parameter :: &
         sclkap=1.e-11
      real, parameter, dimension(3) :: &
        qttt = (/ -3.03869354E-05, -3.03869352E-05, -3.03869353E-05 /) &
       ,qtt  = (/  4.56625601E-03,  4.29277358E-03,  3.38116552E-03 /) &
       ,qt   = (/ -2.88801209E-01, -2.61828868E-01, -1.81335007E-01 /) &
       ,qs   = (/ -1.08670290E-01, -1.05131061E-01, -9.33336309E-02 /) &
       ,qst  = (/  7.90503772E-04,  7.71096940E-04,  1.07270585E-03 /) &
       ,qpt  = (/  1.07813750E-09,  1.00638435E-09,  7.57239852E-10 /) &
       ,qpst = (/  1.41541548E-11,  1.48598578E-11,  3.89226107E-12 /) &
       ,qptt = (/ -1.31383708E-11, -1.31383707E-11, -1.31383708E-11 /)
!
      kappaf1=(r+1000.0)* &
        (exp(sclkap*(prs-pref)* &
              ( s*( qs(kkf)+t* qst(kkf) ) + &
                t*( qt(kkf)+t*(qtt(kkf)+t*qttt(kkf))+ &
                    0.5*(prs+pref)* &
                    (qpt(kkf)+s*qpst(kkf)+t*qptt(kkf)) ) ) ) &
         -1.0)
      end function kappaf1
!
      real function kappaf(t,s,r,prs,kkf)
      implicit none
!
      real    t,s,r,prs
      integer kkf
!
! --- thermobaric compressibility coefficient (integral from prs to pref)
! ---     Sun et.al. (1999) JPO 29 pp 2719-2729.
!
! --- offset limits based on stability estimates from:
! ---     Hallberg (2005) Ocean Modelling 8 pp 279-300.
! --- t: potential temperature (degC); s: salinity (psu);
! --- r: potential density (sigma); prs: pressure; kkf: ref.state
! ---     example: kappaf(4.5,34.5,36.406,1.e7,1) = -0.12301201
! ---     example: kappaf(4.5,34.5,36.406,1.e7,2) = -0.03356404
! ---     example: kappaf(4.5,34.5,36.406,1.e7,3) =  0.05201003
!
      real, parameter, dimension(3) :: &
        toff = (/  0.0,             3.0,            13.0 /) &
       ,soff = (/ 34.5,            35.0,            38.5 /)
!
      kappaf= &
           kappaf1(max(-1.2,         t-toff(kkf) ),   & !Hallberg,T-only: -1.8,0.9
                   max(-3.0,min(1.5, s-soff(kkf))),   & !Hallberg,S-only: -4.2,2.1
                   r,prs,kkf)
      end function kappaf
!
!> Revision history
!>
!> Sep  2004 - added kkf to kappaf, select one of three reference states
!> Aug  2006 - more restrictive kappaf1 offset limits
!> Mar  2009 - modified limits in kappaf
!> Mar  2009 - more accurate kappaf, with potential density
!> Aug  2017 - kappaf and kappaf1 in separate file as internal functions
!-----------------------------------------------------------------------------
