!-----------------------------------------------------------------------------
      integer, parameter :: &
#if defined (EOS_SIG0) 
# if defined (EOS_7T)
        sigver=1  !7-term sigma-0
# elif defined (EOS_9T)
        sigver=3  !9-term sigma-0
# elif defined (EOS_12T)
        sigver=7  !12-term sigma-0
# elif defined (EOS_17T)
        sigver=5  !17-term sigma-0
# endif
#elif defined (EOS_SIG2)
# if defined (EOS_7T)
        sigver=2  !7-term sigma-2
# elif defined (EOS_9T)
        sigver=4  !9-term sigma-2
# elif defined (EOS_12T)
        sigver=8  !12-term sigma-2
# elif defined (EOS_17T)
        sigver=6  !17-term sigma-2
# endif
#endif
!
      real    sig,dsigdt,dsigds,tofsig,sofsig, &
              sigloc,dsiglocdt,dsiglocds,tofsigloc
!
      real    sofsig_a,sofsig_b,sofsig_c
      real    a0,a1,a2,cubr,cubq,cuban,cubrl,cubim
      real    c1l,c2l,c3l,c4l,c5l,c6l,c7l

      real    sig_n,sig_d,sig_q, dsigdt_n,dsigdt_d, dsigds_n,dsigds_d
      real    tofsig_a,tofsig_b,tofsig_c
      real    sigloc_n,sigloc_d,sigloc_q, &
              dsiglocdt_n,dsiglocdt_d, dsiglocds_n,dsiglocds_d
!
      real    r,s,t,pdb,prs
      integer kkf
!
      real, parameter :: &
         aone =1.0, &
         ahalf=1.0/2.0, &
         a3rd =1.0/3.0, athird =a3rd, &
         a4th =1.0/4.0, afourth=a4th
!
!

#if defined (EOS_7T)
# if defined (EOS_SIG0)
! --- coefficients for sigma-0 7-terms (based on Brydon & Sun fit)
      real, parameter :: &
         c1=-1.36471E-01,   & !const. coefficent
         c2= 4.68181E-02,   & !T      coefficent
         c3= 8.07004E-01,   & !   S   coefficent
         c4=-7.45353E-03,   & !T^2    coefficent
         c5=-2.94418E-03,   & !T  S   coefficent
         c6= 3.43570E-05,   & !T^3    coefficent
        rc6= 1.0/c6, &
         c7= 3.48658E-05,   & !T^2S   coefficent
       pref= 0.0              ! reference pressure, Pascals
# elif defined (EOS_SIG2)
! --- coefficients for sigma-2 (based on Brydon & Sun fit)
      real, parameter :: &
         c1= 9.77093E+00,   & !const. coefficent
         c2=-2.26493E-02,   & !T      coefficent
         c3= 7.89879E-01,   & !   S   coefficent
         c4=-6.43205E-03,   & !T^2    coefficent
         c5=-2.62983E-03,   & !T  S   coefficent
         c6= 2.75835E-05,   & !T^3    coefficent
        rc6= 1.0/c6, &
         c7= 3.15235E-05,   & !T^2S   coefficent
       pref= 2000.0e4         !reference pressure, Pascals
# endif
!
! --- HYCOM pressure to bar, for locally referenced equations
      real, parameter :: prs2pb=1.e-5       !Pascals to bar

#elif defined (EOS_9T)
# if defined (EOS_SIG0)
! --- sigma-theta as a function of temp (deg c) and salinity (psu)
! --- (9-term polynomial fit to T:[-2:30],S:[18:38])
!
! --- coefficients for sigma-0.
      real, parameter :: &
         c1=-4.311829E-02,   & !const. coefficent
         c2= 5.429948E-02,   & !T      coefficent
         c3= 8.011774E-01,   & !   S   coefficent
         c4=-7.641336E-03,   & !T^2    coefficent
         c5=-3.258442E-03,   & !T  S   coefficent
         c6= 3.757643E-05,   & !T^3    coefficent
        rc6=1.0/c6, &
         c7= 3.630361E-05,   & !T^2S   coefficent
         c8= 8.675546E-05,   & !   S^2 coefficent
         c9= 3.995086E-06,   & !T  S^2 coefficent
       pref= 0.0               ! reference pressure, Pascals
# elif defined (EOS_SIG2)
! --- coefficients for sigma-2.
      real, parameter :: &
         c1= 9.903308E+00,   & !const. coefficent
         c2=-1.618075E-02,   & !T      coefficent
         c3= 7.819166E-01,   & !   S   coefficent
         c4=-6.593939E-03,   & !T^2    coefficent
         c5=-2.896464E-03,   & !T  S   coefficent
         c6= 3.038697E-05,   & !T^3    coefficent
        rc6= 1.0/c6, &
         c7= 3.266933E-05,   & !T^2S   coefficent
         c8= 1.180109E-04,   & !   S^2 coefficent
         c9= 3.399511E-06,   & !T  S^2 coefficent
       pref= 2000.0e4          ! reference pressure, Pascals
# endif 
!
! --- HYCOM pressure to bar, for locally referenced equations
      real, parameter :: prs2pb=1.e-5       !Pascals to bar

#elif (EOS_12T) 
      real, parameter :: &
         c1= 1.0e-01,       & !not used, required to compile mxkrtm.f
         c2= 1.0e-02,       & !not used, required to compile mxkrtm.f
         c3= 1.0e-03,       & !not used, required to compile mxkrtm.f
         c4= 1.0e-04,       & !not used, required to compile mxkrtm.f
         c5= 1.0e-05,       & !not used, required to compile mxkrtm.f
         c6= 1.0e-06,       & !not used, required to compile mxkrtm.f
         c7= 1.0e-07          !not used, required to compile mxkrtm.f
!
! --- REFERENCE?
!
! --- coefficients for 18-term rational function sigloc().
      real, parameter :: &
       c001=-1.4627567840659594d-01,   & !num. constant    coefficent
       c002= 6.4247392832635697d-02,   & !num.    T        coefficent
       c003= 8.1213979591704621d-01,   & !num.       S     coefficent
       c004=-8.1321489441909698d-03,   & !num.    T^2      coefficent
       c005= 4.5199845091090296d-03,   & !num.    T  S     coefficent
       c006= 4.6347888132781394d-04,   & !num.       S^2   coefficent
       c007= 5.0879498675039621d-03,   & !num. P           coefficent
       c008= 1.6333913018305079d-05,   & !num. P  T        coefficent
       c009= 4.3899924880543972d-06      !num. P     S     coefficent
      real, parameter :: &
       c011= 1.0000000000000000d+00,   & !den. constant    coefficent
       c012= 1.0316374535350838d-02,   & !den.    T        coefficent
       c013= 8.9521792365142522d-04,   & !den.       S     coefficent
       c014=-2.8438341552142710d-05,   & !den.    T^2      coefficent
       c015=-1.1887778959461776d-05,   & !den.    T  S     coefficent
       c016=-4.0163964812921489d-06,   & !den.       S^2   coefficent
       c017= 1.1995545126831476d-05,   & !den. P           coefficent
       c018= 5.5234008384648383d-08,   & !den. P  T        coefficent
       c019= 8.4310335919950873d-09      !den. P     S     coefficent
      real, parameter :: &
       c004x2=c004*2.d0,               & !for dsigdt and dsiglocdt
       c014x2=c014*2.d0,               & !for dsigdt and dsiglocdt
       c006x2=c006*2.d0,               & !for dsigds and dsiglocds
       c016x2=c016*2.d0                  !for dsigds and dsiglocds
      real, parameter :: &
       sqrmin=0.d0,                    & !sqrt arg can't be negative
       sofmin=0.d0                       !salinity can't be negative
! --- reference pressure.
      real, parameter :: prs2pdb=1.d-4     !Pascals to dbar
# if defined (EOS_SIG0)
      real, parameter :: pref=   0.d0      !ref. pressure in Pascals, sigma0
# elif defined (EOS_SIG2)
      real, parameter :: pref=2000.d4      !ref. pressure in Pascals, sigma2
# endif
      real, parameter :: rpdb=pref*prs2pdb !ref. pressure in dbar
! --- coefficients for 12-term rational function sig() at rpdb.
      real, parameter :: &
       c101=c001+rpdb*c007,            & !num. constant    coefficent
       c102=c002+rpdb*c008,            & !num.    T        coefficent
       c103=c003+rpdb*c009               !num.       S     coefficent
      real, parameter :: &
       c111=c011+rpdb*c017,            & !num. constant    coefficent
       c112=c012+rpdb*c018,            & !num.    T        coefficent
       c113=c013+rpdb*c019               !num.       S     coefficent

#elif (EOS_17T)
      real, parameter :: &
         c1= 1.0e-01,       & !not used, required to compile mxkrtm.f
         c2= 1.0e-02,       & !not used, required to compile mxkrtm.f
         c3= 1.0e-03,       & !not used, required to compile mxkrtm.f
         c4= 1.0e-04,       & !not used, required to compile mxkrtm.f
         c5= 1.0e-05,       & !not used, required to compile mxkrtm.f
         c6= 1.0e-06,       & !not used, required to compile mxkrtm.f
         c7= 1.0e-07          !not used, required to compile mxkrtm.f

!
! --- Jackett, McDougall, Feistel, Wright and Griffies (2006),
! --- Algorithms for Density, Potential Temperature, Conservative
! --- Temperature, and the Freezing Temperature of Seawater, JAOT
!
! --- coefficients for 25-term rational function sigloc().
      real, parameter :: &
         c001= 9.9984085444849347d+02,      & !num. constant    coefficent
         c002= 7.3471625860981584d+00,      & !num.    T        coefficent
         c003=-5.3211231792841769d-02,      & !num.    T^2      coefficent
         c004= 3.6492439109814549d-04,      & !num.    T^3      coefficent
         c005= 2.5880571023991390d+00,      & !num.       S     coefficent
         c006= 6.7168282786692355d-03,      & !num.    T  S     coefficent
         c007= 1.9203202055760151d-03,      & !num.       S^2   coefficent
         c008= 1.0000000000000000d+00,      & !den. constant    coefficent
         c009= 7.2815210113327091d-03,      & !den.    T        coefficent
         c010=-4.4787265461983921d-05,      & !den.    T^2      coefficent
         c011= 3.3851002965802430d-07,      & !den.    T^3      coefficent
         c012= 1.3651202389758572d-10,      & !den.    T^4      coefficent
         c013= 1.7632126669040377d-03,      & !den.       S     coefficent
         c014= 8.8066583251206474d-06,      & !den.    T  S     coefficent
         c015= 1.8832689434804897d-10,      & !den.    T^3S     coefficent
         c016= 5.7463776745432097d-06,      & !den.    T  S^1.5 coefficent
         c017= 1.4716275472242334d-09         !den.    T^3S^1.5 coefficent
      real, parameter :: &
         c018= 1.1798263740430364d-02,      & !num. P           coefficent
         c019= 9.8920219266399117d-08,      & !num. P  T^2      coefficent
         c020= 4.6996642771754730d-06,      & !num. P     S     coefficent
         c021= 2.5862187075154352d-08,      & !num. P^2         coefficent
         c022= 3.2921414007960662d-12,      & !num. P^2T^2      coefficent
         c023= 6.7103246285651894d-06,      & !den. P           coefficent
         c024= 2.4461698007024582d-17,      & !den. P^2T^3      coefficent
         c025= 9.1534417604289062d-18         !den. P^3T        coefficent
! --- additional coefficients for dsiglocdt().
      real, parameter :: &
         c031= 7.3471625860981580d+00,      & !num. constant    coefficent
         c032=-1.0642246358568354d-01,      & !num.    T        coefficent
         c033= 1.0947731732944364d-03,      & !num.    T^2      coefficent
         c034= 6.7168282786692355d-03,      & !num.       S     coefficent
         c035= 7.2815210113327090d-03,      & !den. constant    coefficent
         c036=-8.9574530923967840d-05,      & !den.    T        coefficent
         c037= 1.0155300889740728d-06,      & !den.    T^2      coefficent
         c038= 5.4604809559034290d-10,      & !den.    T^3      coefficent
         c039=-8.8066583251206470d-06,      & !den.       S     coefficent
         c040= 5.6498068304414700d-10,      & !den.    T^2S     coefficent
         c041= 2.9432550944484670d-09,      & !den.    T  S^1.5 coefficent
         c042= 1.9784043853279823d-07,      & !num. P  T        coefficent
         c043= 6.5842828015921320d-12,      & !num. P^2T        coefficent
         c044= 7.3385094021073750d-17,      & !den. P^2T^2      coefficent
         c045= 9.1534417604289060d-18         !den. P^3         coefficent
! --- additional coefficients for dsiglocds().
      real, parameter :: &
         c051= 2.5880571023991390d+00,      & !num. constant    coefficent
         c052= 6.7168282786692355d-03,      & !num.    T        coefficent
         c053= 3.8406404111520300d-03,      & !num.       S     coefficent
         c054= 1.7632126669040377d-03,      & !den. constant    coefficent
         c055=-8.8066583251206470d-06,      & !den.    T        coefficent
         c056= 1.8832689434804897d-10,      & !den.    T^3      coefficent
         c057= 8.6195665118148150d-06,      & !den.       S^0.5 coefficent
         c058= 2.2074413208363504d-09,      & !den.    T^2S^0.5 coefficent
         c059= 4.6996642771754730d-06         !num. P           coefficent
!
      real, parameter :: sqrmin=0.d0          !sqrt arg can't be negative
! --- reference pressure.
      real, parameter :: prs2pdb=1.d-4        !Pascals to dbar
# if defined (EOS_SIG0)
      real, parameter :: pref=   0.d0         !ref. pressure in Pascals, sigma0
# elif defined (EOS_SIG2)
      real, parameter :: pref=2000.d4         !ref. pressure in Pascals, sigma2
# endif
      real, parameter :: rpdb=pref*prs2pdb    !ref. pressure in dbar
! --- coefficients for 17-term rational function sig() at rpdb.
      real, parameter :: &
         c101=c001+(c018-c021*rpdb)*rpdb,  & !num. constant    coefficent
         c103=c003+(c019-c022*rpdb)*rpdb,  & !num.    T^2      coefficent
         c105=c005+c020*rpdb,              & !num.       S     coefficent
         c108=c008+c023*rpdb,              & !den. constant    coefficent
         c109=c009-c025*rpdb**3,           & !den.    T        coefficent
         c111=c011-c024*rpdb**2              !den.    T^3      coefficent
! --- additional coefficients for dsigdt().
      real, parameter :: &
         c132=c032+(c042-c043*rpdb)*rpdb,  & !num.    T        coefficent
         c135=c035-c045*rpdb**3,           & !den. constant    coefficent
         c137=c037-c044*rpdb**2              !den.    T^2      coefficent
! --- additional coefficients for dsigds().
      real, parameter :: &
         c151=c051+c059*rpdb                 !num. constant    coefficent
#endif

#if defined (EOS_7T) || (EOS_9T)
!
! --- sub-coefficients for locally referenced sigma
! --- a fit towards Jackett & McDougall (1995)
      real, parameter, dimension(7) :: &
        alphap = (/ -0.1364705627213484   , 0.04681812123458564, &
                     0.80700383913187     ,-0.007453530323180844, &
                    -0.002944183249153631 , 0.00003435702568990446, &
                     0.0000348657661057688 /) &
       ,betap  = (/  0.05064226654169138  ,-0.0003571087848996894, &
                    -0.0000876148051892879, 5.252431910751829e-6, &
                     1.579762259448864e-6 ,-3.466867400295792e-8, &
                    -1.687643078774232e-8 /) &
       ,gammap = (/ -5.526396144304812e-6 , 4.885838128243163e-8, &
                     9.96026931578033e-9  ,-7.251389796582352e-10, &
                    -3.987360250058777e-11, 4.006307891935698e-12, &
                     8.26367520608008e-13 /)
!
#endif


#if defined (EOS_7T)
! --- auxiliary statements for finding root of cubic polynomial
      a0(s,r)=(c1+c3*s-r)*rc6  !constant  coefficient
      a1(s)  =(c2+c5*s  )*rc6  !linear    coefficient
      a2(s)  =(c4+c7*s  )*rc6  !quadratic coefficient
                               !cubic     coefficient is c6*rc6=1.0
      cubq(s)=a3rd*a1(s)-(a3rd*a2(s))**2
      cubr(r,s)=a3rd*(0.5*a1(s)*a2(s)-1.5*a0(s,r))-(a3rd*a2(s))**3
! --- if q**3+r**2>0, water is too dense to yield real root at given
! --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
! --- lowering sigma until a double real root is obtained.
      cuban(r,s)=a3rd*atan2(sqrt(max(0.0,-(cubq(s)**3+cubr(r,s)**2))), &
                              cubr(r,s))
      cubrl(r,s)=sqrt(-cubq(s))*cos(cuban(r,s))
      cubim(r,s)=sqrt(-cubq(s))*sin(cuban(r,s))
!
! --- -----------------
! --- equation of state
! --- -----------------
!
! --- sigma-theta as a function of temp (deg c) and salinity (psu)
! --- (friedrich-levitus, polynomial fit that is cubic in T and linear in S)
!
      sig(t,s)=(c1+c3*s+t*(c2+c5*s+t*(c4+c7*s+c6*t)))
!
! --- d(sig)/dt
      dsigdt(t,s)=(c2+c5*s+2.0*t*(c4+c7*s+1.5*c6*t))
!
! --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7))
!
! --- temp (deg c) as a function of sigma and salinity (psu)
! --- find a cubic polynominal root of t**3+a2*t**2+a1*t+a0=0
      tofsig(r,s)=-cubrl(r,s)+sqrt(3.0)*cubim(r,s)-a3rd*a2(s)
!
! --- salinity (psu) as a function of sigma and temperature (deg c)
      sofsig(r,t)=(r-c1-t*(c2+t*(c4+c6*t)))/(c3+t*(c5+c7*t))
!
! --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
! --- t: potential temperature; s: psu; prs: pressure
      c1l(prs)=alphap(1)+prs2pb*prs*(betap(1)+prs2pb*prs*gammap(1))
      c2l(prs)=alphap(2)+prs2pb*prs*(betap(2)+prs2pb*prs*gammap(2))
      c3l(prs)=alphap(3)+prs2pb*prs*(betap(3)+prs2pb*prs*gammap(3))
      c4l(prs)=alphap(4)+prs2pb*prs*(betap(4)+prs2pb*prs*gammap(4))
      c5l(prs)=alphap(5)+prs2pb*prs*(betap(5)+prs2pb*prs*gammap(5))
      c6l(prs)=alphap(6)+prs2pb*prs*(betap(6)+prs2pb*prs*gammap(6))
      c7l(prs)=alphap(7)+prs2pb*prs*(betap(7)+prs2pb*prs*gammap(7))
      sigloc(t,s,prs)=c1l(prs)+c3l(prs)*s+ &
             t*(c2l(prs)+c5l(prs)*s+t*(c4l(prs)+c7l(prs)*s+c6l(prs)*t))
      dsiglocdt(t,s,prs)=(c2l(prs)+c5l(prs)*s+ &
             2.0*t*(c4l(prs)+c7l(prs)*s+1.5*c6l(prs)*t))
      dsiglocds(t,s,prs)=(c3l(prs)+t*(c5l(prs)+t*c7l(prs)))

#elif defined (EOS_9T)
! --- auxiliary statements for finding root of cubic polynomial
      a0(s,r)=(c1+s*(c3+s*c8)-r)*rc6  !constant  coefficient
      a1(s)  =(c2+s*(c5+s*c9)  )*rc6  !linear    coefficient
      a2(s)  =(c4+s* c7        )*rc6  !quadratic coefficient
                                      !cubic     coefficient is c6*rc6=1.0
      cubq(s)  =a3rd*     a1(s)                   -(a3rd*a2(s))**2
      cubr(r,s)=a3rd*(0.5*a1(s)*a2(s)-1.5*a0(s,r))-(a3rd*a2(s))**3
! --- if q**3+r**2>0, water is too dense to yield real root at given
! --- salinitiy. setting q**3+r**2=0 in that case is equivalent to
! --- lowering sigma until a double real root is obtained.
      cuban(r,s)=a3rd*atan2(sqrt(max(0.0,-(cubq(s)**3+cubr(r,s)**2))), &
                            cubr(r,s))
      cubrl(r,s)=sqrt(-cubq(s))*cos(cuban(r,s))
      cubim(r,s)=sqrt(-cubq(s))*sin(cuban(r,s))
!
! --- -----------------
! --- equation of state
! --- -----------------
!
! --- sigma-theta as a function of temp (deg c) and salinity (psu)
! --- (polynomial fit that is cubic in T and quadratic in S)
!
      sig(t,s)=(c1+s*(c3+s* c8)+ &
                   t*(c2+s*(c5+s*c9)+t*(c4+s*c7+t*c6)))
!
! --- d(sig)/dt
      dsigdt(t,s)=(c2+s*(c5+s*c9)+2.0*t*(c4+s*c7+1.5*t*c6))
!
! --- d(sig)/ds
      dsigds(t,s)=(c3+t*(c5+t*c7)+2.0*s* c8)
!
! --- temp (deg c) as a function of sigma and salinity (psu)
! --- find a cubic polynominal root of t**3+a2*t**2+a1*t+a0=0
      tofsig(r,s)=-cubrl(r,s)+sqrt(3.0)*cubim(r,s)-a3rd*a2(s)
!
! --- salinity (psu) as a function of sigma and temperature (deg c)
! --- find a quadratic polynominal root of a*s**2+b*s+c=0
      sofsig_a(r,t)=(c8+t* c9)                !quadratic coefficient
      sofsig_b(r,t)=(c3+t*(c5+t* c7))         !linear    coefficient
      sofsig_c(r,t)=(c1+t*(c2+t*(c4+t*c6))-r) !constant  coefficient
      sofsig(r,t)=(2.0*sofsig_c(r,t))/ &
                  (-sofsig_b(r,t) &
                   -sign(sqrt(sofsig_b(r,t)**2- &
                              4.0*sofsig_a(r,t)*sofsig_c(r,t)), &
                         sofsig_b(r,t)))
!
! --- locally referenced sigma, a fit towards Jackett & McDougall (1995)
! --- t: potential temperature; s: psu; prs: pressure
      c1l(prs)=alphap(1)+prs2pb*prs*(betap(1)+prs2pb*prs*gammap(1))
      c2l(prs)=alphap(2)+prs2pb*prs*(betap(2)+prs2pb*prs*gammap(2))
      c3l(prs)=alphap(3)+prs2pb*prs*(betap(3)+prs2pb*prs*gammap(3))
      c4l(prs)=alphap(4)+prs2pb*prs*(betap(4)+prs2pb*prs*gammap(4))
      c5l(prs)=alphap(5)+prs2pb*prs*(betap(5)+prs2pb*prs*gammap(5))
      c6l(prs)=alphap(6)+prs2pb*prs*(betap(6)+prs2pb*prs*gammap(6))
      c7l(prs)=alphap(7)+prs2pb*prs*(betap(7)+prs2pb*prs*gammap(7))
      sigloc(t,s,prs)=c1l(prs)+c3l(prs)*s+ &
             t*(c2l(prs)+c5l(prs)*s+t*(c4l(prs)+c7l(prs)*s+c6l(prs)*t))
      dsiglocdt(t,s,prs)=(c2l(prs)+c5l(prs)*s+ &
             2.0*t*(c4l(prs)+c7l(prs)*s+1.5*c6l(prs)*t))
      dsiglocds(t,s,prs)=(c3l(prs)+t*(c5l(prs)+t*c7l(prs)))

#elif (EOS_12T)
!
! --- -----------------
! --- equation of state
! --- -----------------
!
! --- sigma at rpdb (dbar) as a function of temp (deg c) and salinity (psu)
!
      sig_n(t,s) = c101+(c102+c004*t+c005*s)*t  + &
                        (c103+       c006*s)*s
      sig_d(t,s) = c111+(c112+c014*t+c015*s)*t  + &
                        (c113       +c016*s)*s
      sig_q(t,s) = aone/sig_d(t,s)
      sig(  t,s) = sig_n(t,s)*sig_q(t,s)
!
! --- d(sig)/dt
      dsigdt_n(t,s) = c102+c004x2*t+c005*s
      dsigdt_d(t,s) = c112+c014x2*t+c015*s
      dsigdt(  t,s) = (dsigdt_n(t,s)- &
                       dsigdt_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
!
! --- d(sig)/ds
      dsigds_n(t,s) = c103+c005*t+c006x2*s
      dsigds_d(t,s) = c113+c015*t+c016x2*s
      dsigds(  t,s) = (dsigds_n(t,s)- &
                       dsigds_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
!
! --- temp (deg c) as a function of sigma and salinity (psu)
! --- find a quadratic polynominal root of a*t**2+b*t+c=0
      tofsig_a(r,s)=(   c004 - &
                     r* c014   )                  !quadratic coefficient
      tofsig_b(r,s)=(  (c102+      c005*s) - &
                     r*(c112+      c015*s)  )     !linear    coefficient
      tofsig_c(r,s)=(  (c101+(c103+c006*s)*s) - &
                     r*(c111+(c113+c016*s)*s)  )  !constant  coefficient
      tofsig(r,s)=( -tofsig_b(r,s) &
                    -sqrt(max(sqrmin, &
                                  tofsig_b(r,s)**2 - &
                              4.0*tofsig_a(r,s)*tofsig_c(r,s))) ) / &
                  (2.0*tofsig_a(r,s))
!
! --- salinity (psu) as a function of sigma and temperature (deg c)
! --- find a quadratic polynominal root of a*s**2+b*s+c=0
      sofsig_a(r,t)=(   c006 - &
                     r* c016   )                  !quadratic coefficient
      sofsig_b(r,t)=(  (c103+      c005*t) - &
                     r*(c113+      c015*t)  )     !linear    coefficient
      sofsig_c(r,t)=(  (c101+(c102+c004*t)*t) - &
                     r*(c111+(c112+c014*t)*t)  )  !constant  coefficient
      sofsig(r,s)=max(sofmin, &
                      ( -sofsig_b(r,s) &
                        +sqrt(max(sqrmin, &
                                      sofsig_b(r,s)**2 - &
                                  4.0*sofsig_a(r,s)*sofsig_c(r,s))) ) / &
                      (2.0*sofsig_a(r,s)) )
!
! --- locally referenced sigma, using the 18-term equation of state.
! --- t: potential temperature; s: psu; prs: pressure
!
      sigloc_n(t,s,pdb) = c001+(c002+c004*t+c005*s)*t  + &
                               (c003+       c006*s)*s  + &
                               (c007+c008*t+c009*s)*pdb
      sigloc_d(t,s,pdb) = c011+(c012+c014*t+c015*s)*t  + &
                               (c013       +c016*s)*s  + &
                               (c017+c018*t+c019*s)*pdb
      sigloc_q(t,s,pdb) = aone/sigloc_d(t,s,pdb)
      sigloc(  t,s,prs)=sigloc_n(t,s,prs*prs2pdb)* &
                        sigloc_q(t,s,prs*prs2pdb)
!
! --- d(sig)/dt
      dsiglocdt_n(t,s,pdb) = c002+c004x2*t+c005*s+c008*pdb
      dsiglocdt_d(t,s,pdb) = c012+c014x2*t+c015*s+c018*pdb
      dsiglocdt(  t,s,prs)=(dsiglocdt_n(t,s,prs*prs2pdb)- &
                            dsiglocdt_d(t,s,prs*prs2pdb)* &
                               sigloc_n(t,s,prs*prs2pdb)* &
                               sigloc_q(t,s,prs*prs2pdb) ) * &
                               sigloc_q(t,s,prs*prs2pdb)
!
! --- d(sig)/ds
      dsiglocds_n(t,s,pdb) = c003+c005*t+c006x2*s+c009*pdb
      dsiglocds_d(t,s,pdb) = c013+c015*t+c016x2*s+c019*pdb
      dsiglocds(  t,s,prs)=(dsiglocds_n(t,s,prs*prs2pdb)- &
                            dsiglocds_d(t,s,prs*prs2pdb)* &
                               sigloc_n(t,s,prs*prs2pdb)* &
                               sigloc_q(t,s,prs*prs2pdb) ) * &
                               sigloc_q(t,s,prs*prs2pdb)
!

#elif defined (EOS_17T)
!
! --- -----------------
! --- equation of state
! --- -----------------
!
! --- sigma at rpdb (dbar) as a function of pot.temp (deg c) and salinity (psu)
!
      sig_n(t,s) = c101 + t*(c002+t*(c103+t*c004)) + &
                          s*(c105-t*c006+s*c007)
      sig_d(t,s) = c108 + t*(c109+t*(c010+t*(c111+t*c012))) + &
                          s*(c013-t*(c014+t*t*c015) + &
        sqrt(max(sqrmin,s))*(c016+t*t*c017))
      sig_q(t,s) = aone/sig_d(t,s)
      sig(  t,s) = sig_n(t,s)*sig_q(t,s) - 1000.0
!
! --- d(sig)/dt
      dsigdt_n(t,s) = c031 + t*(c132+t*c033) - &
                             s* c034
      dsigdt_d(t,s) = c135 + t*(c036+t*(c137+t*c038)) + &
                             s*(c039-t*t*c040+ &
                   sqrt(max(sqrmin,s))*t*c041 )
      dsigdt(  t,s) = (dsigdt_n(t,s)- &
                       dsigdt_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
!
! --- d(sig)/ds
! --- additional coefficients for dsigds().
      dsigds_n(t,s) = c151 - t*c052 + s*c053
      dsigds_d(t,s) = c054 +       t*(c055-t*t*c056) + &
                 sqrt(max(sqrmin,s))*(c057+t*t*c058)
      dsigds(  t,s) = (dsigds_n(t,s)- &
                       dsigds_d(t,s)*sig_n(t,s)*sig_q(t,s))*sig_q(t,s)
!
! --- temp (deg c) as a function of sigma and salinity (psu)
! --- NOT AVAILABLE AS AN EXPRESSION - DO NOT USE
      tofsig(r,s)=99.0
!
! --- salinity (psu) as a function of sigma and temperature (deg c)
! --- NOT AVAILABLE AS AN EXPRESSION - DO NOT USE
      sofsig(r,t)=99.0
!
! --- locally referenced sigma, using the 25-term equation of state.
! --- t: potential temperature (degC); s: salinity (psu); prs: pressure (dbar)
      sigloc_n(t,s,pdb) =      c001 + &
                            t*(c002 + &
                            t*(c003 + &
                            t* c004  )) + &
                            s*(c005 - &
                            t* c006 + &
                            s* c007  ) + &
                          pdb*(c018 + &
                          t*t* c019 + &
                            s* c020 - &
                          pdb*(c021 + &
                          t*t* c022  ))
      sigloc_d(t,s,pdb) =      c008 + &
                            t*(c009 + &
                            t*(c010 + &
                            t*(c011 + &
                            t* c012  ))) + &
                            s*(c013 - &
                            t*(c014 + &
                          t*t* c015  ) + &
          sqrt(max(sqrmin,s))*(c016 + &
                          t*t* c017  )) + &
                          pdb*(c023 - &
                   pdb*t*(t*t* c024 + &
                          pdb* c025  ))
      sigloc_q(t,s,pdb) = aone/sigloc_d(t,s,pdb)
      sigloc(t,s,prs)=sigloc_n(t,s,prs*prs2pdb)* &
                      sigloc_q(t,s,prs*prs2pdb) - 1000.0 
!
! --- d(sig)/dt
      dsiglocdt_n(t,s,pdb) =   c031 + &
                            t*(c032 + &
                            t* c033  ) - &
                            s* c034 + &
                        pdb*t*(c042 - &
                          pdb* c043  )
      dsiglocdt_d(t,s,pdb) =   c035 + &
                            t*(c036 + &
                            t*(c037 + &
                            t* c038  )) + &
                            s*(c039 - &
                          t*t* c040 + &
        sqrt(max(sqrmin,s))*t* c041  ) - &
                 pdb*pdb*(t*t* c044 + &
                          pdb* c045  )
      dsiglocdt(t,s,prs)=(dsiglocdt_n(t,s,prs*prs2pdb)- &
                          dsiglocdt_d(t,s,prs*prs2pdb)* &
                             sigloc_n(t,s,prs*prs2pdb)* &
                             sigloc_q(t,s,prs*prs2pdb) ) * &
                             sigloc_q(t,s,prs*prs2pdb)
!
! --- d(sig)/ds
      dsiglocds_n(t,s,pdb) =   c051 - &
                            t* c052 + &
                            s* c053 + &
                          pdb* c059
      dsiglocds_d(t,s,pdb) =   c054 + &
                            t*(c055 - &
                          t*t* c056  ) + &
          sqrt(max(sqrmin,s))*(c057 + &
                          t*t* c058  )
      dsiglocds(t,s,prs)=(dsiglocds_n(t,s,prs*prs2pdb)- &
                          dsiglocds_d(t,s,prs*prs2pdb)* &
                             sigloc_n(t,s,prs*prs2pdb)* &
                             sigloc_q(t,s,prs*prs2pdb) ) * &
                             sigloc_q(t,s,prs*prs2pdb)
!

#endif


!
!> Revision history
!>
!> May  2000 - conversion to SI units
!> Jul  2000 - removed rarely used functions, constants via parameter
!> Jan  2002 - removed geometery functions
!> Dec  2002 - new thermobaricity fit with toff=0.0,soff=34.0
!> Jun  2003 - removed sigma4
!> Jun  2003 - added locally referenced sigma
!> Sep  2004 - added kkf to kappaf, select one of three reference states
!> Aug  2006 - more restrictive kappaf1 offset limits
!> Sep  2006 - 9-term polynominal fit to T:[-2:30],S:[18:38]
!> May  2007 - added sigver
!> Mar  2009 - modified limits in kappaf
!> Mar  2009 - more accurate kappaf, with potential density
!> Oct  2010 - 17-term rational function equation of state
!> Aug  2017 - moved kappaf to internal_kappaf.h
!> Jan  2019 - define MACRO to choose EOS at compile time (EOS_SIGX; EOS_XXT)
!-----------------------------------------------------------------------------
