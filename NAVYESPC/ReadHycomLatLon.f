c===================================================
c      subroutine readHycomLatLon(plat,plon,itdm,jtdm)
      subroutine readHycomLatLon(plat,plon,qlat,qlon,itdm,jtdm)

      implicit none

      REAL*4 plat(itdm,jtdm),plon(itdm,jtdm)
      REAL*4 qlat(itdm,jtdm),qlon(itdm,jtdm)
      REAL*4,  allocatable :: pad(:)

      integer npad
      integer ios,nrecl
      integer i,j
      CHARACTER*240 cfilea
      integer itdm,jtdm

      cfilea = 'regional.grid.a'

      npad = 4096 - MOD(itdm*jtdm,4096)
      if(npad.eq.4096) npad=0

      allocate(pad(npad))

      INQUIRE( IOLENGTH=nrecl) plon,pad
c      print *,"itdm,jtdm=",itdm,jtdm
c      print *,"npad=",npad
c      print *,"nrecl=",nrecl


      open(unit=11,file=cfilea, form='unformatted', status='old',
     *         access='direct', recl=nrecl, iostat=ios)


      IF (ios.ne.0) THEN
        print *,"error in reading regional.grid.a"
        call exit(1)
      endif

      read(11,rec=1,iostat=ios) plon
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, plon"
        call exit(2)
      endif

      read(11,rec=2,iostat=ios) plat
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, plat"
        call exit(3)
      endif
     
      read(11,rec=3,iostat=ios) qlon
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, qlon"
        call exit(2)
      endif

      read(11,rec=4,iostat=ios) qlat
      if (ios.ne.0) then
        print *,"error in reading regional.grid.a, qlat"
        call exit(3)
      endif
 
      do j=1,jtdm
      do i=1,itdm
        if(plon(i,j).ge.360) plon(i,j)=plon(i,j)-360.
        if(qlon(i,j).ge.360) qlon(i,j)=qlon(i,j)-360.
      enddo
      enddo

      print *,'readHycomLatLon,plat, min,max=',
     &    minval(plat),maxval(plat)
      print *,'readHycomLatLon,plon, min,max=',
     &    minval(plon),maxval(plon)


c      print *,"**** readHycomLatLon, lat_hycom ***"
c      do j=1,jtdm
c        print *, "j=", j
c        write(*,12)(plat(i,j),i=1,itdm)
c      enddo

c      print *,"**** readHycomLatLon, lon_hycom ***"
c      do j=1,jtdm
c        print *, "j=", j
c        write(*,12)(plon(i,j),i=1,itdm)
c      enddo
c 12   format(10F12.5/(10F12.5))


      if(allocated(pad)) deallocate(pad)

      return
      end subroutine readHycomLatLon


