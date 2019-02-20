c==================================================
      subroutine export_from_hycom_tiled(tmx, fieldName)
      use mod_xc  ! HYCOM communication interface
      use mod_cb_arrays

      implicit none
c      include 'common_blocks.h'
c
      real              tmx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      character(len=30) fieldName
c
      integer i,j,k
      real    hfrz,tfrz,t2f,ssfi,tmxl,smxl,umxl,vmxl
      real    dp1,usur1,vsur1,psur1,dp2,usur2,vsur2,psur2,thksur
c
      integer jja
c

#if defined(ARCTIC)
c --- Arctic (tripole) domain, top row is replicated (ignore it)
      jja = min( jj, jtdm-1-j0 )
#else
      jja = jj
#endif

      tmx(:,:)=0.

c==>  export to atm,ice
      if(fieldName .eq. 'sst' ) then
        do j=1,jj
        do i= 1,ii
            tmx(i,j) = 0.5*(temp(i,j,1,1)+temp(i,j,1,2))
        enddo
        enddo
c==>  export to ice
      else if(fieldName .eq. 'ssu' .or.
     &        fieldName .eq. 'ssv'     ) then
        call xctilr(u(    1-nbdy,1-nbdy,1,1),1,2*kk, 1,1, halo_uv)
        call xctilr(ubavg(1-nbdy,1-nbdy,  1),1,   2, 1,1, halo_uv)
        call xctilr(v(    1-nbdy,1-nbdy,1,1),1,2*kk, 1,1, halo_vv)
        call xctilr(vbavg(1-nbdy,1-nbdy,  1),1,   2, 1,1, halo_vv)


        do j=1,jja
        do i= 1,ii
          if     (ishlf(i,j).eq.1) then
c ---       average currents over top thkcdw meters
            thksur = onem*min( thkcdw, depths(i,j) )
            usur1  = 0.0
            vsur1  = 0.0
            psur1  = 0.0
            usur2  = 0.0
            vsur2  = 0.0
            psur2  = 0.0
            do k= 1,kk
              dp1   = min( dp(i,j,k,1), max( 0.0, thksur-psur1 ) )
              usur1 = usur1 + dp1*(u(i,j,k,1)+u(i+1,j,k,1))
              vsur1 = vsur1 + dp1*(v(i,j,k,1)+v(i,j+1,k,1))
#if defined(STOKES)
              usur1 = usur1 + dp1*(usd(i,j,k)+usd(i+1,j,k))
              vsur1 = vsur1 + dp1*(vsd(i,j,k)+vsd(i,j+1,k))
#endif
              psur1 = psur1 + dp1
              dp2   = min( dp(i,j,k,2), max( 0.0, thksur-psur2 ) )
              usur2 = usur2 + dp2*(u(i,j,k,2)+u(i+1,j,k,2))
              vsur2 = vsur2 + dp2*(v(i,j,k,2)+v(i,j+1,k,2))
#if defined(STOKES)
              usur2 = usur2 + dp2*(usd(i,j,k)+usd(i+1,j,k))
              vsur2 = vsur2 + dp2*(vsd(i,j,k)+vsd(i,j+1,k))
#endif
              psur2 = psur2 + dp2
              if     (min(psur1,psur2).ge.thksur) then
                exit
              endif
            enddo !k
            umxl  = 0.25*( usur1/psur1 + ubavg(i,  j,1) +
     &                                   ubavg(i+1,j,1) +
     &                     usur2/psur2 + ubavg(i,  j,2) +
     &                                   ubavg(i+1,j,2)  )
            vmxl  = 0.25*( vsur1/psur1 + vbavg(i,j,  1) +
     &                                   vbavg(i,j+1,1) +
     &                     vsur2/psur2 + vbavg(i,j,  2) +
     &                                   vbavg(i,j+1,2)  )
           if     (fieldName .eq. 'ssu') then
#ifndef ESPC_NOCANONICAL_CONVERT
c            rotate to e-ward
             tmx(i,j)=cos(pang(i,j))*umxl - sin(pang(i,j))*vmxl
#else
             tmx(i,j)=umxl
#endif
           else !'ssv'
#ifndef ESPC_NOCANONICAL_CONVERT
c            rotate to n-ward
             tmx(i,j)=cos(pang(i,j))*vmxl + sin(pang(i,j))*umxl
#else
             tmx(i,j)=vmxl
#endif
           endif
          endif !ishlf:else
        enddo !i
        enddo !j
c ---   Smooth surface ocean velocity fields
#if defined(ARCTIC)
        call xctila( tmx,1,1,halo_pv)
#endif
        call psmooth(tmx,0,0,ishlf,util1)

      else if(fieldName .eq. 'sss' ) then
c       sea surface salinity  (ppt)
        do j=1,jj
        do i= 1,ii
          tmx(i,j) = 0.5*(saln(i,j,1,1)+saln(i,j,1,2))
        enddo
        enddo
      else if(fieldName .eq. 'ssh' ) then
c       sea surface height in m
        do j=1,jj
        do i= 1,ii
          tmx(i,j) = 1./g*srfhgt(i,j)
        enddo
        enddo
      else if(fieldName .eq. 'ssfi' ) then
c       Oceanic Heat Flux Available to Sea Ice
        do j=1,jj
        do i= 1,ii
c ---       quantities for available freeze/melt heat flux
c ---       relax to tfrz with e-folding time of icefrq time steps
c ---       assuming the effective surface layer thickness is hfrz
c ---       multiply by dpbl(i,j)/hfrz to get the actual e-folding time
            hfrz = min( thkfrz*onem, dpbl(i,j) )
            t2f  = (spcifh*hfrz)/(baclin*icefrq*g)
c ---       average both available time steps, to avoid time splitting.
            smxl = 0.5*(saln(i,j,1,1)+saln(i,j,1,2))
            tmxl = 0.5*(temp(i,j,1,1)+temp(i,j,1,2))
            tfrz = tfrz_0 + smxl*tfrz_s  !salinity dependent freezing point
            ssfi = (tfrz-tmxl)*t2f       !W/m^2 into ocean
c
            tmx(i,j) = max(-1000.0,min(1000.0,ssfi))
        enddo
        enddo
      else if(fieldName .eq. 'mlt' ) then
c       Ocean Mixed Layer Thickness
        do j=1,jj
        do i= 1,ii
          tmx(i,j) = dpbl(i,j)*qonem
        enddo
        enddo
      else if(fieldName .eq. 'sbhflx' ) then
c      sensible heat flux
        do j=1,jj
        do i= 1,ii
            tmx(i,j) = exp_sbhflx(i,j)
        enddo
        enddo
      else if(fieldName .eq. 'lthflx' ) then
c      latent heat flux
        do j=1,jj
        do i= 1,ii
            tmx(i,j) = exp_lthflx(i,j)
        enddo
        enddo

      else if(fieldName .eq. 'dummy_ocn' ) then
c       dummy_ocn
        do j=1,jj
        do i= 1,ii
c          tmx(i,j)=plat(i,j)
            tmx(i,j) = 1.
        enddo
        enddo

      else
c       unknown export fieldName
        do j=1,jj
        do i= 1,ii
          tmx(i,j) = 0.
        enddo
        enddo
      endif

      return
      end subroutine export_from_hycom_tiled

