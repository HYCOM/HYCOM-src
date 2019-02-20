      subroutine inikpp
      use mod_xc         ! HYCOM communication interface
      use mod_cb_arrays  ! HYCOM saved arrays
!
! --- hycom version 1.0
      implicit none
!
      integer    nzehat,nustar
      parameter (nzehat=890,nustar=192)
!
      real, dimension (0:nzehat+1,0:nustar+1) :: &
       wmt             & ! momentum velocity scale table
      ,wst               ! scalar   velocity scale table
      common/kppltr/ wmt,wst
      save  /kppltr/
!
! -------------------------------------------------------------------
! --- initialize large, mc williams, doney kpp vertical mixing scheme
! -------------------------------------------------------------------
!
      integer i,j
      real zehat,zeta,am,cm,c22,zetam,as,c33,zetas,usta
!
      data am,cm,c22,zetam/1.257,8.380,16.0,-0.2/
      data as,c33,zetas/-28.86,16.0,-1.0/
!
# include "stmt_fns.h"
!
! --- 'vonk'      = von karman constant
! --- 'zmin,zmax' = zehat limits for velocity scale lookup table, m**3/s**3
! --- 'umin,umax' = ustar limits for velocity scale lookup table
! --- 'epsilon'   = vertical coordinate scale factor
!
      vonk   =  0.4
      zmin   = -0.4e-6
      zmax   =  0.0
      umin   =  0.0
      umax   =  0.16
      epsilon=  0.1
!
! --- construct the velocity-scale lookup tables
!
      deltaz = (zmax-zmin)/(nzehat+1)
      deltau = (umax-umin)/(nustar+1)
!
      do i=0,nzehat+1
        zehat=deltaz*i+zmin
        do j=0,nustar+1
          usta=deltau*j+umin
          zeta=zehat/(usta**3+epsil)
          if (zehat.ge.0.) then
            wmt(i,j)=vonk*usta/(1.+c11*zeta)
            wst(i,j)=wmt(i,j)
          else
            if (zeta.gt.zetam) then
              wmt(i,j)=vonk*usta*(1.-c22*zeta)**afourth
            else
              wmt(i,j)=vonk*(am*usta**3-cm*zehat)**athird
            endif
            if (zeta.gt.zetas) then
              wst(i,j)=vonk*usta*(1.-c33*zeta)**ahalf
            else
              wst(i,j)=vonk*(as*usta**3-cs*zehat)**athird
            endif
          endif
        enddo
      enddo
!
! --- set derived constants
      vtc=sqrt(.2/cs/epsilon)/vonk**2/ricr
      cg=cstar*vonk*(cs*vonk*epsilon)**athird
      dp0enh=2.0*dp00
!
      return
      end
!
!
!> Revision history:
!>
!> May  2001 - increased nustar and umax by a factor of 4
