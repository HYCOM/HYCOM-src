      module mod_cb_arrays_nuopc_glue
      
        ! The purpose of this module is to make select common block entities
        ! available as public module variables. Doing this through an extra
        ! Fortran module is cleaner than simply including the entyre common
        ! blocks header in the actual glue module. Further, the common blocks
        ! header is in fixed format, and that means that which ever unit pulls
        ! it in via include also needs to be in fixed format... not something
        ! that would be nice for the actual hycom_nuopc_glue module.
      
        use mod_xc
        use mod_cb_arrays       ! HYCOM saved arrays

        implicit none
    
        private
  
        ! grid c-grid staggering
        public plon     ! lon at the p pts
        public plat     ! lat at the p pts
        public qlon     ! lon at the p pts
        public qlat     ! lat at the p pts
        public scp2     ! grid cell area at the p pts
        public scux     ! grid cell area at the p pts
        public scvy     ! grid cell area at the p pts

        ! import from ATM
        public cpl_taux,    imp_taux, imp_taue ! xstress [Pa], eastward stress [Pa]
        public cpl_tauy,    imp_tauy, imp_taun ! ystress [Pa], northward stress [Pa]
        public cpl_wndspd,  imp_wndspd    ! wind speed [m s-1]
        public cpl_ustara,  imp_ustara    ! friction speed [m s-1]
        public cpl_airtmp,  imp_airtmp    ! air temperature [K]
        public cpl_vapmix,  imp_vapmix    ! specific humidity [kg kg-1]
        public cpl_swflx,   imp_swflx     ! shortwave net flux [W m-2]
        public cpl_lwmdnflx,imp_lwdflx    ! longwave down flux [W m-2]
        public cpl_lwmupflx,imp_lwuflx    ! longwave up   flux [W m-2]
        public cpl_precip,  imp_precip    ! precipitation [m s-1]
        public cpl_surtmp,  imp_surtmp    ! air surface temperature [C]
        public cpl_seatmp,  imp_seatmp    ! sea surface temperature [C]
        public cpl_latflx,  imp_latflx    ! latent heat flux   [W m-2]
        public cpl_sensflx, imp_sensflx   ! sensible heat flux [W m-2]
        public cpl_orivers, imp_orivers   ! ocean runoffs [kg m-2 s-1]
        public cpl_irivers, imp_irivers   ! ice   runoffs [kg m-2 s-1]
        
        
        ! import from SEA-ICE
        public cpl_sic,     sic_import    ! Sea Ice Concentration
        public cpl_sitx,    sitx_import   ! Sea Ice X-Stress
        public cpl_sity,    sity_import   ! Sea Ice Y-Stress
        public cpl_siqs,    siqs_import   ! Solar Heat Flux thru Ice to Ocean
        public cpl_sifh,    sifh_import   ! Ice Freezing/Melting Heat Flux
        public cpl_sifs,    sifs_import   ! Ice Freezing/Melting Salt Flux
        public cpl_sifw,    sifw_import   ! Ice Net Water Flux
        public cpl_sit,     sit_import    ! Sea Ice Temperature
        public cpl_sih,     sih_import    ! Sea Ice Thickness
        public cpl_siu,     siu_import    ! Sea Ice X-Velocity
        public cpl_siv,     siv_import    ! Sea Ice Y-Velocity

        ! native exports
        public temp           ! temp at the p pts
        public u, v           ! velocity components
        public ubavg, vbavg   ! barotropic velocities
        public srfhgt         ! sea surface height, g*ssh(m)
        public dpbl           ! turbulent boundary layer depth
        public saln           ! salinity
        public pang           ! ANGLET
        public dhde           ! eastward ssh slope
        public dhdn           ! northward ssh slope
        public umxl           ! eastward u velocity
        public vmxl           ! northward v velocity
        public tml            ! export T average over one coupling sequence
        public sml            ! export S average over one coupling sequence
        public sshm           ! export SSH average over one coupling sequence
        public frzh           ! export of freezing and melting potential 
        
        ! native imports
        public covice         ! ice coverage (rel.units)
        public thkice         ! grid-cell avg. ice thknss (m)
        public temice         ! ice surface temperature
        public fswice         ! swv  flux under ice
        public flxice         ! heat flux under ice
        public sflice         ! salt flux under ice
        public wflice         ! water flux under ice
        public si_c           ! ice concentration   on p-grid from coupler
        public si_h           ! ice thickness       on p-grid from coupler
        public si_t           ! ice temperature     on p-grid from coupler
        public si_u           ! ice u-velocity      on p-grid from coupler
        public si_v           ! ice v-velocity      on p-grid from coupler
        public si_tx          ! x-stress  under ice on p-grid from coupler
        public si_ty          ! y-stress  under ice on p-grid from coupler
        
        ! ocn coupling frequency
        public ocn_cpl_frq    ! in days in CESM
        public icefrq         ! in time step in HYCOM
        
        ! scalars and flags
        public g, onem, qonem, thkfrz, baclin, spcifh, pcp_fact
        public tfrz_0, tfrz_s
        public iceflg, icmflg, sstflg,itest,jtest,ltripolar
        
      end module
