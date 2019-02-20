!-----------------------------------------------------------------------------
! --- START OF REGION AND TILING SPECIFIC PARAMETERS
! --- Static memory version, not used when macro /* RELO */ is set (dynamic memory)
! --- See: README.dimensions and README.OpenMP for more details.
!
! --- itdm  = total grid dimension in i direction
! --- jtdm  = total grid dimension in j direction
! --- kdm   =       grid dimension in k direction
      integer    itdm,jtdm,kdm
      parameter (itdm=500,jtdm=382,kdm=41)  ! GLBT0.72
!
! --- iqr   = maximum number of tiles in i direction
! --- jqr   = maximum number of tiles in j direction
      integer    iqr,jqr
      parameter (iqr=10,jqr=10)  ! multiple tiles (TYPE=ompi or mpi or shmem)
!
! --- idm   = maximum single tile grid dimension in i direction
! --- jdm   = maximum single tile grid dimension in j direction
      integer    idm,jdm
!!!!!!parameter (idm=itdm,jdm=jtdm)  ! always works if enough memory
      parameter (idm= 250,jdm= 191)  ! NMPI=4,8,16,24,32,40,47,64
!
! --- mxthrd= maximum number of OpenMP threads
      integer    mxthrd
      parameter (mxthrd=1)  ! NOMP=0,1
!
! --- kkwall= grid dimension in k direction for wall relax arrays
! --- kknest= grid dimension in k direction for nest relax arrays
      integer    kkwall,kknest
      parameter (kkwall=  1)  ! must be 1 or kdm
      parameter (kknest=  1)  ! must be 1 or kdm
!
! --- kkmy25= grid dimension in k direction for M-Y 2.5 arrays
      integer    kkmy25
      parameter (kkmy25= -1)  ! must be -1 or kdm
!
! --- nlgiss= size of lookup table for GISS
      integer    nlgiss
      parameter (nlgiss=  1)  ! must be 1 (no GISS) or 762
!
! --- mxtrcr= maximum number of tracers
      integer    mxtrcr
      parameter (mxtrcr=1)
!
! --- natm  = number of saved atmospheric fields
      integer    natm
      parameter (natm=2)      ! must be 2 (high freq.) or 4 (monthly)
!
! --- max_nsteps_batrop = maximum barotropic steps per baroclinic time step
      integer    max_nsteps_batrop
      parameter (max_nsteps_batrop = 128)
!
! ---   END OF REGION AND TILING SPECIFIC PARAMETERS
!-----------------------------------------------------------------------------
