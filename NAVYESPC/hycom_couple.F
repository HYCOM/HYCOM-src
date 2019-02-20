      module hycom_couple

      use mod_xc  ! HYCOM communication interface
      use mod_cb_arrays


      implicit none

      integer idim_size,jdim_size
      integer, dimension(:,:,:), allocatable :: deBList

      real, dimension(:,:), allocatable :: lon_e
      real, dimension(:,:), allocatable :: lat_e
      real, dimension(:,:), allocatable :: mask_e

      real, dimension(:,:), allocatable :: lon_q
      real, dimension(:,:), allocatable :: lat_q


      contains

c===================================================

      subroutine hycom_couple_init(nPets,rc) 

      implicit none


      real,   allocatable, dimension(:,:) :: tmx
      real*4, allocatable, dimension(:,:) :: alon,alat
      real*4, allocatable, dimension(:,:) :: qlon,qlat

      real,   allocatable, dimension(:,:) :: tmp_e


      integer i,j,rc
      integer nPets

      if(mnproc.eq.1) print *,"hycom_couple_init called..."
      rc=0
c---------------------
c     grid size
      idim_size=itdm
      jdim_size=jtdm

c-----------------------
c     deBlockList

c     directly from HYCOM 

# ifdef ESPC_COUPLE
      IF (.not.allocated(deBList)) THEN
          allocate ( deBList(2,2,nPets ) )
      END IF


c     should be something like the following
c      do ij=1,nPets
c        deBList(1,1,ij)=1+i0
c        deBList(2,1,ij)=1+j0
c        deBList(1,2,ij)=ii+i0
c        deBList(2,2,ij)=jj+j0
c      enddo

      do i=1,nPets
          deBList(1,1,i)=deBlockList(1,1,i)
          deBList(2,1,i)=deBlockList(2,1,i)
          deBList(1,2,i)=deBlockList(1,2,i)
          deBList(2,2,i)=deBlockList(2,2,i)
      enddo


      if(mnproc.eq.1) then
        print *,'itdm,jtdm=',itdm,jtdm
        print *,'hycom,deBList BL11 BL21 BL12 BL22'
        do i =1,nPets
          write(*,555)i,deBList(1,1,i),deBList(2,1,i),
     &      deBList(1,2,i),deBList(2,2,i), 
     &      deBList(1,2,i)-deBList(1,1,i)+1, 
     &      deBList(2,2,i)-deBList(2,1,i)+1

        enddo
      endif
 555  format(I4,4I8,3x,2I8)

# endif

c-----------------------
c     lat/lon/mask

      if(mnproc.eq.1) then
        allocate(lon_e(itdm,jtdm))
        allocate(lat_e(itdm,jtdm))
        allocate(mask_e(itdm,jtdm))
        allocate(lon_q(itdm,jtdm))
        allocate(lat_q(itdm,jtdm))
      else
        allocate(lon_e(1,1))
        allocate(lat_e(1,1))
        allocate(mask_e(1,1))
        allocate(lon_q(1,1))
        allocate(lat_q(1,1))
      endif
      allocate(tmp_e(itdm,jtdm))

      if(mnproc.eq.1) print *,'readHycomLatLon check0...' 

c     for plon -> lon_e
c     for plat -> lat_e
c     get and return lons, lats to all nodes

c      vland=-99999.9 !data void marker

      allocate(tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))



c lon
cx      lon_e(:,:)=9999.
cx      call xcaget(lon_e,plon,1)
      call xcaget(tmp_e,plon,1)
      if(mnproc.eq.1) lon_e(:,:)=tmp_e(:,:)


      if(mnproc.eq.1) then
        do j=1,jtdm
        do i=1,itdm
          if(lon_e(i,j) .ge. 360) then
            lon_e(i,j)=lon_e(i,j)-360.
          else
            lon_e(i,j)=lon_e(i,j)
          endif
        enddo
        enddo
      endif

c lat
cx      lat_e(:,:)=9999.
cx      call xcaget(lat_e,plat,1)
      call xcaget(tmp_e,plat,1)
      if(mnproc.eq.1) lat_e(:,:)=tmp_e(:,:)



c sea/land mask
      tmx(:,:)=0.
      do j= 1,jj
        do i= 1,ii
            tmx(i,j) = ishlf(i,j)
cx            tmx(i,j) = ip(i,j)
        enddo
      enddo
      mask_e(:,:)=0.
cx      call xcaget(mask_e,tmx,1)
      call xcaget(tmp_e,tmx,1)
      if(mnproc.eq.1) mask_e(:,:)=tmp_e(:,:)



      if(mnproc.eq.1) then
        allocate(alon(itdm,jtdm))
        allocate(alat(itdm,jtdm))
        allocate(qlon(itdm,jtdm))
        allocate(qlat(itdm,jtdm))
      else
        allocate(alon(1,1))
        allocate(alat(1,1))
        allocate(qlon(1,1))
        allocate(qlat(1,1))
      endif

      if(mnproc.eq.1) then
        print *,'readHycomLatLon check...'

c       read hycom regional.grid.a
c        call readHycomLatLon(alat,alon,itdm,jtdm)
        call readHycomLatLon(alat,alon,qlat,qlon,itdm,jtdm)

        do j=1,jtdm
        do i=1,itdm
          if(lon_e(i,j).eq.0 .and. lat_e(i,j).eq.0) then
            lon_e(i,j)=alon(i,j)
            lat_e(i,j)=alat(i,j)
          endif
          lon_q(i,j)=qlon(i,j)
          lat_q(i,j)=qlat(i,j)
        enddo
        enddo

      endif




c      vland=0.0 !restore to default

      if(allocated(tmx)) deallocate(tmx)

      if(allocated(alon)) deallocate(alon)
      if(allocated(alat)) deallocate(alat)
      if(allocated(qlon)) deallocate(qlon)
      if(allocated(qlat)) deallocate(qlat)

      if(allocated(tmp_e)) deallocate(tmp_e)

      if(mnproc.eq.1)  print *,"hycom_couple_init, end..."

 
      return
      end subroutine hycom_couple_init




c===================================================
      subroutine set_hycom_import_flag(k,fieldName)

c      use ocn_couple_impexp
c      use mod_cb_arrays

      implicit none

      integer k 
      character(len=30) fieldName

      if(mnproc.eq.1) print *,
     &    "set_hycom_import_flag start...,k,name=",k,fieldName

      if(k.eq.1) then
        cpl_taux=.false.
        cpl_tauy=.false.
        cpl_u10=.false.
        cpl_v10=.false.
        cpl_wndspd=.false.
        cpl_ustara=.false.
        cpl_airtmp=.false.
        cpl_vapmix=.false.
        cpl_swflx_net=.false.
        cpl_lwflx_net=.false.
        cpl_swflx_net2down=.false.
        cpl_lwflx_net2down=.false.
        cpl_swflxd=.false.
        cpl_lwflxd=.false.
        cpl_mslprs=.false.
        cpl_precip=.false.
        cpl_surtmp=.false.
        cpl_seatmp=.false.
        cpl_sbhflx=.false.
        cpl_lthflx=.false.
        cpl_sic=.false.
        cpl_sitx=.false.
        cpl_sity=.false.
        cpl_siqs=.false.
        cpl_sifh=.false.
        cpl_sifs=.false.
        cpl_sifw=.false.
        cpl_sit=.false.
        cpl_sih=.false.
        cpl_siu=.false.
        cpl_siv=.false.
      endif



      if(fieldName .eq. 'taux10' ) then
        cpl_taux=.true.

      else if(fieldName .eq. 'tauy10' ) then
        cpl_tauy=.true.
        if     (.not.cpl_taux) then
          if(mnproc.eq.1) print *,"error - tauy before taux"
          call xcstop('(set_hycom_import_flag)')
                 stop '(set_hycom_import_flag)'
        endif !error

      else if(fieldName .eq. 'u10' ) then
        cpl_u10=.true.

      else if(fieldName .eq. 'v10' ) then
        cpl_v10=.true.
        if     (.not.cpl_u10) then
          if(mnproc.eq.1) print *,"error - v10 before u10"
          call xcstop('(set_hycom_import_flag)')
                 stop '(set_hycom_import_flag)'
        endif !error

      else if(fieldName .eq. 'wndspd10' ) then
        cpl_wndspd=.true.

      else if(fieldName .eq. 'ustara10' ) then
        cpl_ustara=.true.

      else if(fieldName .eq. 'airtmp' ) then
        cpl_airtmp=.true.

      elseif(fieldName .eq. 'airhum' ) then
        cpl_vapmix=.true.

      else if(fieldName .eq. 'swflx_net' ) then
        cpl_swflx_net=.true.

      else if(fieldName .eq. 'lwflx_net' ) then
        cpl_lwflx_net=.true.

      else if(fieldName .eq. 'swflx_net2down' ) then
        cpl_swflx_net2down=.true.

      else if(fieldName .eq. 'lwflx_net2down' ) then
        cpl_lwflx_net2down=.true.

      else if(fieldName .eq. 'swflxd' ) then
        cpl_swflxd=.true.

      else if(fieldName .eq. 'lwflxd' ) then
        cpl_lwflxd=.true.
      else if(fieldName .eq. 'mslprs' ) then
       cpl_mslprs=.true.
      else if(fieldName .eq. 'prcp' ) then
        cpl_precip=.true.

      else if(fieldName .eq. 'gt' ) then
        cpl_surtmp=.true.
        cpl_seatmp=.true.

      else if(fieldName .eq. 'sbhflx' ) then
        cpl_sbhflx=.true.

      else if(fieldName .eq. 'lthflx' ) then
        cpl_lthflx=.true.

      else if(fieldName .eq. 'sic' ) then
c       import ice concentration
        cpl_sic=.true.

      else if(fieldName .eq. 'sitx' ) then
c       import ice x-stress
        cpl_sitx=.true.

      else if(fieldName .eq. 'sity' ) then
c       import ice y-stress
        cpl_sity=.true.

      else if(fieldName .eq. 'siqs' ) then
c       import solar thru grid cell ave.
        cpl_siqs=.true.

      else if(fieldName .eq. 'sifh' ) then
c       import freeze, melt, H. Flux
        cpl_sifh=.true.

      else if(fieldName .eq. 'sifs' ) then
c       import salt flux
        cpl_sifs=.true.

      else if(fieldName .eq. 'sifw' ) then
c       import water flux
        cpl_sifw=.true.

      else if(fieldName .eq. 'sit_sfc' ) then
c       import sea ice temperature
        cpl_sit=.true.

      else if(fieldName .eq. 'sih' ) then
c       import sea ice thickness
        cpl_sih=.true.

      else if(fieldName .eq. 'siu' ) then
c       import sea ice x-velocity
        cpl_siu=.true.

      else if(fieldName .eq. 'siv' ) then
c       import sea ice y-velocity
        cpl_siv=.true.

      endif  !if fieldName


      if(mnproc.eq.1) print *,"import_hycom end..."

      return
      end subroutine set_hycom_import_flag


c====================================================
      subroutine export_from_hycom_deb(tlb,tub,expData,
     &  fieldName,show_minmax)
c      use mod_xc  ! HYCOM communication interface
c      use mod_cb_arrays

      implicit none
c
c      integer           k
c      real              mgrid(ii,jj)
      integer tlb(2),tub(2)
      real expData(tlb(1):tub(1),tlb(2):tub(2))

      character(len=30) fieldName
      real, allocatable, dimension(:,:) :: ocn_msk
      real, allocatable, dimension(:,:) :: field_tmp
      real, allocatable, dimension(:,:) :: tmx

      integer i,j,jja
      logical show_minmax

!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))


c
      if(mnproc.eq.1) print *,"export_from_hycom_deb start..."
c      print *,"idm,jdm,nbdy,ii,jj=",mnproc,idm,jdm,nbdy,ii,jj

      call export_from_hycom_tiled(util2,fieldName)  !can't use util1

#if defined(ARCTIC)
c --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif

      if(fieldName .eq. 'sst' ) then
        do j=1,jj
        do i= 1,ii
#ifndef ESPC_NOCANONICAL_CONVERT
c           canonical unit conversion: sst (C) -> (K)
            util2(i,j) = util2(i,j)+273.15d0
#endif

        enddo
        enddo
      endif

      expData(:,:)=0.
      do j=1,jja
      do i=1,ii
c        mgrid(i,j)=util2(i,j)
        expData(i+i0,j+j0)=util2(i,j)
      enddo
      enddo


      if(show_minmax) then
        if(mnproc.eq.1) then
          allocate(ocn_msk(itdm,jtdm))
          allocate(field_tmp(itdm,jtdm))
        else
          allocate(ocn_msk(1,1))
          allocate(field_tmp(1,1))
        endif
        allocate(tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))


c       sea/land mask
        tmx(:,:)=0.
        do j= 1,jja
        do i= 1,ii
            tmx(i,j) = ishlf(i,j)
cx            tmx(i,j) = ip(i,j)
        enddo
        enddo

        call xcaget(ocn_msk,tmx,1)
c        call xcsync(no_flush)

        tmx(:,:)=0.
        do j= 1,jja
        do i= 1,ii
            tmx(i,j) = expData(i+i0,j+j0)
        enddo
        enddo

        call xcaget(field_tmp,tmx,1)
c        call xcsync(no_flush)

        if(mnproc.eq.1) then

          write(*,992)trim(fieldName),
     &    maxval(field_tmp,MASK=ocn_msk.eq.1 ),
     &    minval(field_tmp,MASK=ocn_msk.eq.1 ),
     &    sum(field_tmp,MASK=ocn_msk.eq.1 )/
     &    count(ocn_msk.eq.1 )
 992      format('export_from_hycom_deb,max,min,mean=',A10,3E23.15)

c          write(*,992) trim(fieldName),maxval(mgrid),minval(mgrid)
c 992      format('export_from_hycom,max,min=',A10,2E12.4)
           print *,"export_from_hycom_deb end..."

        endif

c       test check pang
        call xcaget(ocn_msk,pang,1)
        if(mnproc.eq.1) then
          print *,'export_from_hycom pang, min,max=',
     &     minval(ocn_msk),maxval(ocn_msk)
        endif

        if(allocated(ocn_msk)) deallocate(ocn_msk)
        if(allocated(field_tmp)) deallocate(field_tmp)
        if(allocated(tmx)) deallocate(tmx)

      endif

      if(mnproc.eq.1) print *,"export_from_hycom_deb end..."

      return
      end subroutine export_from_hycom_deb

c==================================================
      subroutine import_to_hycom_deb(tlb,tub,impData,
     &  fieldName,show_minmax,data_init_flag)

c      use mod_xc  ! HYCOM communication interface
c      use ocn_couple_impexp
c      use mod_cb_arrays

      implicit none
c      include 'common_blocks.h'
c
      character(len=30) fieldName
c      integer           k
      integer tlb(2),tub(2)
      real impData(tlb(1):tub(1),tlb(2):tub(2))
c
      integer i,j,mcnt
      real    uij,vij

      real, allocatable, dimension(:,:) :: ocn_msk
      real, allocatable, dimension(:,:) :: field_tmp
      real, allocatable, dimension(:,:) :: tmx

      real, parameter :: sstmin = -1.8
      real, parameter :: sstmax = 35.0

      integer ierr
      logical show_minmax
      integer jja

      real    albw,degtorad
      logical data_init_flag


!   (1+i0,ii+i0) could be the subset of (tlb(1),tub(1))
!   (1+j0,jja+j0) == (tlb(2),tub(2))

c      if(mnproc.eq.1) print *,"import_to_hycom_deb start...,k,name=",k,
c     &   fieldName
      
c      if(k.eq.1 .and. mnproc.eq.1) print *,"w0,w1..=",
c     &    w0,w1,w2,w3

#if defined(ARCTIC)
c --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif


      if(show_minmax) then
cjc-01262014
        if(mnproc.eq.1) then
          allocate(ocn_msk(itdm,jtdm))
          allocate(field_tmp(itdm,jtdm))
        else
          allocate(ocn_msk(1,1))
          allocate(field_tmp(1,1))
        endif

        allocate(tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy))

c       sea/land mask
        tmx(:,:)=0.
        do j= 1,jja
        do i= 1,ii
          tmx(i,j) = ishlf(i,j)
cx          tmx(i,j) = ip(i,j)
        enddo
        enddo

        call xcaget(ocn_msk,tmx,1)
c        call xcsync(no_flush)

        tmx(:,:)=0.
        do j= 1,jja
        do i= 1,ii
!          tmx(i,j) = mgrid(i,j)
          tmx(i,j) = impData(i+i0,j+j0)
        enddo
        enddo

        call xcaget(field_tmp,tmx,1)
c       call mpi_barrier(mpi_comm_hycom,ierr)
c       call xcsync(no_flush)

        if(mnproc.eq.1) then

          write(*,992)trim(fieldName),
     &    maxval(field_tmp,MASK=ocn_msk.eq.1 ),
     &    minval(field_tmp,MASK=ocn_msk.eq.1 ),
     &    sum(field_tmp,MASK=ocn_msk.eq.1 )/
     &    count(ocn_msk.eq.1 )
c 992      format('import_to_hycom_deb,max,min,mean=',A10,3E12.4)
 992      format('import_to_hycom_deb,max,min,mean=',A10,3E23.15)

        endif

        if(allocated(ocn_msk)) deallocate(ocn_msk)
        if(allocated(field_tmp)) deallocate(field_tmp)
        if(allocated(tmx)) deallocate(tmx)

      endif




c==>  import from atm
c=================================================

      if(fieldName .eq. 'taux10' ) then

c       import xstress: Pa


        do j=1,jja
        do i=1,ii
!          taux(i,j,l0)=mgrid(i,j)
          if(ishlf(i,j).eq.1) then
            taux(i,j,l0)=impData(i+i0,j+j0)
          else
            taux(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(taux(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
#endif
        call xctilr(taux(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)

c---
c=================================================

      else if(fieldName .eq. 'tauy10' ) then

c       import ystress: Pa

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            tauy(i,j,l0)=impData(i+i0,j+j0)
          else
            tauy(i,j,l0)=0.
          endif
        enddo
        enddo

        do j=1,jja
        do i=1,ii

            uij = taux(i,j,l0)
            vij = tauy(i,j,l0)
c           rotate to (x,y)ward
            taux(i,j,l0)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
            tauy(i,j,l0)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
        enddo !i
        enddo !j

#if defined(ARCTIC)
        call xctila(taux(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
        call xctila(tauy(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
#endif
        call xctilr(taux(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)
        call xctilr(tauy(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)

c---
c=================================================

      else if(fieldName .eq. 'u10' ) then

c       import u wind at 10m height: ms-1


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            taux(i,j,l0)=impData(i+i0,j+j0)
          else
            taux(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(taux(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
#endif
        call xctilr(taux(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)

c---
c=================================================

      else if(fieldName .eq. 'v10' ) then

c       import v wind at 10m height: ms-1

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            tauy(i,j,l0)=impData(i+i0,j+j0)
          else
            tauy(i,j,l0)=0.
          endif
        enddo
        enddo

        do j=1,jja
        do i=1,ii
          uij = taux(i,j,l0)
          vij = tauy(i,j,l0)
c         rotate to (x,y)ward
          taux(i,j,l0)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
          tauy(i,j,l0)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
        enddo !i
        enddo !j

#if defined(ARCTIC)
        call xctila(taux(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
        call xctila(tauy(1-nbdy,1-nbdy,l0),1,1,            halo_pv)
#endif
        call xctilr(taux(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)
        call xctilr(tauy(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_pv)

c---
c=================================================

      else if(fieldName .eq. 'wndspd10' ) then

c       import wind speed: m s-1 

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            wndspd(i,j,l0)=impData(i+i0,j+j0)
          else
            wndspd(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      wndspd(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'ustara10' ) then

c       import friction speed: m s-1

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            ustara(i,j,l0)=impData(i+i0,j+j0)
          else
            ustara(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      ustara(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'airtmp' ) then
c        cpl_airtmp=.true.

c       import air temperature
c       canonical unit conversion: airtmp (K) -> (C)
        do j=1,jja
        do i=1,ii
          impData(i+i0,j+j0)=impData(i+i0,j+j0)-273.15
        enddo
        enddo

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            airtmp(i,j,l0)=impData(i+i0,j+j0)
          else
            airtmp(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      airtmp(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================


      elseif(fieldName .eq. 'airhum' ) then

c       import specific humidity: kg kg-1 

c       convert from specific humidity to mixing ratio

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
!!Alex      flxflg.eq.4 => mixing ratio
            vapmix(i,j,l0)=impData(i+i0,j+j0)/(1.-impData(i+i0,j+j0))
!!Alex      flxflg.eq.5 => specific humidity
            if (flxflg.eq.5) vapmix(i,j,l0)=impData(i+i0,j+j0)
          else
            vapmix(i,j,l0)=0.01
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      vapmix(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================


      else if(fieldName .eq. 'swflx_net' ) then

c       import sw flux: w m-2

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            swflx(i,j,l0)=impData(i+i0,j+j0)
          else
            swflx(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      swflx(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'swflx_net2down'
     &      .or. fieldName .eq. 'swflxd' ) then

c       import downward sw flux: w m-2


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            swflx(i,j,l0)=impData(i+i0,j+j0)
          else
            swflx(i,j,l0)=0.
          endif
        enddo
        enddo



        if     (albflg.ne.0) then  !swflx is Qswdn
c ---     use the same method as on forfun.F
c ---     convert swflx to net shortwave into the ocean
c ---     shortwave through sea ice is handled separately
          if     (albflg.eq.1) then
            do j=1,jja
            do i=1,ii
c              if(ishlf(i,j).eq.1) then
                swflx(i,j,l0) = swflx(i,j,l0)*(1.0-0.09)  !NAVGEM albedo
c              else
c                swflx(i,j,l0) = 0.
c              endif
            enddo
            enddo
          else   !albflg.eq.2
            degtorad = 4.d0*atan(1.d0)/180.d0
            do j=1,jja
            do i=1,ii

c ---           latitudinally-varying ocean albedo (Large and Yeager, 2009)
c ---           5.8% at the equator and 8% at the poles
                albw = ( 0.069 - 0.011*cos(2.0*degtorad*plat(i,j) ) )
c                if(ishlf(i,j).eq.1) then
                  swflx(i,j,l0) = swflx(i,j,l0)*(1.0-albw)
c                else
c                  swflx(i,j,l0) = 0.
c                endif
            enddo
            enddo
          endif !albflg
        endif

#if defined(ARCTIC)
        call xctila(      swflx(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'lwflx_net' ) then

c       import lw flux: w m-2
c       canonical unit conversion: lwflx_net (upward) -> (downward)
        do j=1,jja
        do i=1,ii
          impData(i+i0,j+j0)=impData(i+i0,j+j0)*(-1.)
        enddo
        enddo

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            radflx(i,j,l0)=impData(i+i0,j+j0)
          else
            radflx(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      radflx(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'lwflx_net2down'
     &      .or. fieldName .eq. 'lwflxd' ) then

c       import downward lw flux: w m-2
c       +ve into ocean


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            radflx(i,j,l0)=impData(i+i0,j+j0)
          else
            radflx(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      radflx(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================


      else if(fieldName .eq. 'prcp' ) then

c       import precip: m s-1
c       canonical unit conversion: prcp (kg_m-2_s-1)-> m_s-1 
        do j=1,jja
        do i=1,ii
          impData(i+i0,j+j0)=impData(i+i0,j+j0)*(0.001)

        enddo
        enddo


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            precip(i,j,l0)=impData(i+i0,j+j0)
          else
            precip(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      precip(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

c---
c=================================================

      else if(fieldName .eq. 'gt' ) then
c       canonical unit conversion: gt (K) -> (C)
        do j=1,jja
        do i=1,ii
          impData(i+i0,j+j0)=impData(i+i0,j+j0)-273.15
        enddo
        enddo


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            surtmp(i,j,l0)=impData(i+i0,j+j0)
          else
            surtmp(i,j,l0)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      surtmp(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif


        if     (sstflg.ne.3) then  !use atmos. sst as "truth"
          do j=1,jja
          do i=1,ii

              seatmp(i,j,l0) = max( sstmin, 
     &            min(surtmp(i,j,l0), sstmax ) )
          enddo
          enddo
#if defined(ARCTIC)
          call xctila(      seatmp(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

        endif

c---
c=================================================

      else if(fieldName .eq. 'mslprs' ) then
c       import sea level pressure anomaly: Pa


        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            mslprs(i,j,l0)=impData(i+i0,j+j0)
          else
            mslprs(i,j,l0)=1000.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      mslprs(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif
        call xctilr(mslprs(1-nbdy,1-nbdy,l0),1,1, nbdy,nbdy, halo_ps)


c---
c=================================================

c==>    import from sea ice
      else if(fieldName .eq. 'sic' ) then
c       import ice concentration

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sic_import(i,j)=impData(i+i0,j+j0)
          else
            sic_import(i,j)=0.
          endif
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      sic_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sitx' ) then
c       import ice x-stress

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sitx_import(i,j)=impData(i+i0,j+j0)
          else
            sitx_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sitx_import(1-nbdy,1-nbdy),1,1,halo_pv)
#endif

      else if(fieldName .eq. 'sity' ) then
c       import ice y-stress

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sity_import(i,j)=impData(i+i0,j+j0)
          else
            sity_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sity_import(1-nbdy,1-nbdy),1,1,halo_pv)
#endif

      else if(fieldName .eq. 'siqs' ) then
c       import solar thru grid cell ave.
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            siqs_import(i,j)=impData(i+i0,j+j0)
          else
            siqs_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      siqs_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sifh' ) then
c       import freeze, melt, H. Flux
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sifh_import(i,j)=impData(i+i0,j+j0)
          else
            sifh_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sifh_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sifs' ) then
c       import salt flux
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sifs_import(i,j)=impData(i+i0,j+j0)
          else
            sifs_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sifs_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sifw' ) then
c       import water flux
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sifw_import(i,j)=impData(i+i0,j+j0)
          else
            sifw_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sifw_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sit_sfc' ) then
c       import sea ice temperature
#ifndef ESPC_NOCANONICAL_CONVERT
c       canonical unit conversion: sit_sfc (K) -> (C)
        do j=1,jja
        do i=1,ii
          impData(i+i0,j+j0)=impData(i+i0,j+j0)-273.15d0
        enddo
        enddo
#endif

        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sit_import(i,j)=impData(i+i0,j+j0)
          else
            sit_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sit_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'sih' ) then
c       import sea ice thickness
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            sih_import(i,j)=impData(i+i0,j+j0)
          else
            sih_import(i,j)=0.
          endif
        enddo
        enddo
#if defined(ARCTIC)
        call xctila(      sih_import(1-nbdy,1-nbdy),1,1,halo_ps)
#endif

      else if(fieldName .eq. 'siu' ) then
c       import sea ice x-velocity
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            siu_import(i,j)=impData(i+i0,j+j0)
          else
            siu_import(i,j)=0.
          endif
        enddo
        enddo

      else if(fieldName .eq. 'siv' ) then
c       import sea ice y-velocity
        do j=1,jja
        do i=1,ii
          if(ishlf(i,j).eq.1) then
            siv_import(i,j)=impData(i+i0,j+j0)
          else
            siv_import(i,j)=0.
          endif
        enddo
        enddo

#ifndef ESPC_NOCANONICAL_CONVERT
c ---   rotate to (x,y)ward

        do j=1,jja
        do i=1,ii
          uij = siu_import(i,j)
          vij = siv_import(i,j)
c         rotate to (x,y)ward
          siu_import(i,j)=cos(pang(i,j))*uij + sin(pang(i,j))*vij
          siv_import(i,j)=cos(pang(i,j))*vij - sin(pang(i,j))*uij
        enddo !i
        enddo !j
#endif


#if defined(ARCTIC)
        call xctila(      siu_import(1-nbdy,1-nbdy),1,1,halo_pv)
        call xctila(      siv_import(1-nbdy,1-nbdy),1,1,halo_pv)
#endif


      endif  !if fieldName



      if(mnproc.eq.1) print *,"import_hycom_deb end..."

      return
      end subroutine import_to_hycom_deb
c================================
      subroutine ocn_import_forcing()

!      include 'common_blocks.h'
!!Alex      use mod_cb_arrays

      integer i,j,m,n,jja

#if defined(ARCTIC)
! --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif

!# ifdef ESPC_ATM
# if defined(ESPC_ATM) || defined(ESPC_NAVGEM) || defined(ESPC_DATA_ATM)

      if (lwflag.eq.2) then
!       radflx is defined as net lwflx+swflx, +ve into ocean
!       the imported radflx is actually net lw flux, +ve into ocean

        do j=1,jja
        do i= 1,ii
          radflx(i,j,l0)=radflx(i,j,l0)+swflx(i,j,l0)
        enddo
        enddo

#if defined(ARCTIC)
        call xctila(      radflx(1-nbdy,1-nbdy,l0),1,1,halo_ps)
#endif

      endif

# endif

      if(.not.cpl_sic  .or. .not.cpl_sitx .or. .not.cpl_sity .or. 
     &    .not.cpl_siqs .or. .not.cpl_sifh .or. 
     &    .not.cpl_sifs .or. .not.cpl_sifw .or. 
     &    .not.cpl_sit  .or. .not.cpl_sih  .or.
     &    .not.cpl_siu  .or. .not. cpl_siv) then

         if(mnproc.eq.1) print *, 
     &   'warning... no feedback from CICE to HYCOM(ocn_import_forcing)'
         return

      endif
      do j=1,jja
      do i= 1,ii


      if(ishlf(i,j).eq.1) then  !standard ocean point

        if     (iceflg.ge.2 .and. icmflg.ne.3) then
          covice(i,j) = sic_import(i,j) !Sea Ice Concentration
          si_c(i,j) = sic_import(i,j) !Sea Ice Concentration
          if     (covice(i,j).gt.0.0) then
             si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
             si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
            fswice(i,j) =  siqs_import(i,j) !Solar Heat Flux thru Ice to Ocean
            flxice(i,j) =  fswice(i,j) + 
     &                     sifh_import(i,j) !Ice Freezing/Melting Heat Flux
            sflice(i,j) =  sifs_import(i,j)*1.e3 !Ice Salt Flux
            wflice(i,j) =  sifw_import(i,j) !Ice freshwater Flux
            temice(i,j) =   sit_import(i,j) !Sea Ice Temperature
              si_t(i,j) =   sit_import(i,j) !Sea Ice Temperature
            thkice(i,j) =   sih_import(i,j) !Sea Ice Thickness
              si_h(i,j) =   sih_import(i,j) !Sea Ice Thickness
              si_u(i,j) =   siu_import(i,j) !Sea Ice X-Velocity
              si_v(i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
          else
             si_tx(i,j) = 0.0
             si_ty(i,j) = 0.0
            fswice(i,j) = 0.0
            flxice(i,j) = 0.0
            sflice(i,j) = 0.0
            wflice(i,j) = 0.0
            temice(i,j) = 0.0
              si_t(i,j) = 0.0
            thkice(i,j) = 0.0
              si_h(i,j) = 0.0
              si_u(i,j) = 0.0
              si_v(i,j) = 0.0
          endif !covice

        elseif (iceflg.ge.2 .and. icmflg.eq.3) then
           si_c(i,j) =  sic_import(i,j) !Sea Ice Concentration
           if     (si_c(i,j).gt.0.0) then
             si_tx(i,j) = -sitx_import(i,j) !Sea Ice X-Stress into ocean
             si_ty(i,j) = -sity_import(i,j) !Sea Ice Y-Stress into ocean
              si_h(i,j) =   sih_import(i,j) !Sea Ice Thickness
              si_t(i,j) =   sit_import(i,j) !Sea Ice Temperature
              si_u(i,j) =   siu_import(i,j) !Sea Ice X-Velocity
              si_v(i,j) =   siv_import(i,j) !Sea Ice Y-Velocity
           else
             si_tx(i,j) = 0.0
             si_ty(i,j) = 0.0
              si_h(i,j) = 0.0
              si_t(i,j) = 0.0
              si_u(i,j) = 0.0
              si_v(i,j) = 0.0
           endif !covice
        endif !iceflg>=2 (icmflg)

      endif ! if(ishlf
      enddo
      enddo


#if defined(ARCTIC)

! --- update last active row of array
!jcx     call xctila( sic_import,1,1,halo_ps)  !Sea Ice Concentration
!jcx      call xctila(sitx_import,1,1,halo_pv)  !Sea Ice X-Stress
!jcx      call xctila(sity_import,1,1,halo_pv)  !Sea Ice Y-Stress
!jcx      call xctila(siqs_import,1,1,halo_ps)  !Solar Heat Flux thru Ice to Ocean
!jcx      call xctila(sifh_import,1,1,halo_ps)  !Ice Freezing/Melting Heat Flux
!jcx      call xctila(sifs_import,1,1,halo_ps)  !Ice Freezing/Melting Salt Flux
!jcx      call xctila(sifw_import,1,1,halo_ps)  !Ice Net Water Flux
!jcx      call xctila( sit_import,1,1,halo_ps)  !Sea Ice Temperature
!jcx      call xctila( sih_import,1,1,halo_ps)  !Sea Ice Thickness
!jcx      call xctila( siu_import,1,1,halo_pv)  !Sea Ice X-Velocity
!jcx      call xctila( siv_import,1,1,halo_pv)  !Sea Ice Y-Velocity

      if     (iceflg.ge.2 .and. icmflg.ne.3) then
        call xctila(covice,1,1,halo_ps)  !Sea Ice Concentration
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(fswice,1,1,halo_ps)  !Solar Heat Flux thru Ice to Ocean
        call xctila(flxice,1,1,halo_ps)  !Ice Freezing/Melting Heat Flux
        call xctila(sflice,1,1,halo_ps)  !Ice Salt Flux
        call xctila(wflice,1,1,halo_ps)  !Ice Freshwater Flux
        call xctila(temice,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(thkice,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      elseif (iceflg.ge.2 .and. icmflg.eq.3) then
        call xctila(  si_c,1,1,halo_ps)  !Sea Ice Concentration
        call xctila( si_tx,1,1,halo_pv)  !Sea Ice X-Stress into ocean
        call xctila( si_ty,1,1,halo_pv)  !Sea Ice Y-Stress into ocean
        call xctila(  si_h,1,1,halo_ps)  !Sea Ice Thickness
        call xctila(  si_t,1,1,halo_ps)  !Sea Ice Temperature
        call xctila(  si_u,1,1,halo_pv)  !Sea Ice X-Velocity
        call xctila(  si_v,1,1,halo_pv)  !Sea Ice Y-Velocity
      endif


#endif
! --- Smooth Sea Ice velocity fields
      call psmooth(si_u,0,0,ishlf,util1)
      call psmooth(si_v,0,0,ishlf,util1)

#if defined(ARCTIC)
      call xctila(si_u,1,1,halo_pv)
      call xctila(si_v,1,1,halo_pv)
#endif
!      call xctilr(si_u,1,1, nbdy,nbdy, halo_pv)
!      call xctilr(si_v,1,1, nbdy,nbdy, halo_pv)

! --- copy back from si_ to _import for archive_ice
      do j= 1,jja
      do i= 1,ii
          if     (si_c(i,j).gt.0.0) then
            siu_import(i,j) = si_u(i,j) !Sea Ice X-Velocity
            siv_import(i,j) = si_v(i,j) !Sea Ice Y-Velocity
          endif !si_c
      enddo !i
      enddo !j

      end subroutine ocn_import_forcing

      subroutine hycom_couple_final()

      if(allocated(deBList)) deallocate(deBList)
      if(allocated(lon_e)) deallocate(lon_e)
      if(allocated(lat_e)) deallocate(lat_e)
      if(allocated(mask_e)) deallocate(mask_e)
      if(allocated(lon_q)) deallocate(lon_q)
      if(allocated(lat_q)) deallocate(lat_q)

      return
      end subroutine hycom_couple_final



c========================================================

      end module hycom_couple

c========================================================
