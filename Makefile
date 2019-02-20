#
# --- HYCOM 2.2 makefile 
#
# --- Stand-alone HYCOM, or HYCOM ESMF component, or HYCOM+CICE.
#
# --- Tunable parameters in ../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F90 .o

.F90:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

include config/$(ARCH)_$(TYPE)

MODS =   mod_dimensions.o mod_xc.o mod_za.o mod_cb_arrays.o mod_pipe.o \
         mod_incupd.o \
         mod_floats.o mod_stokes.o mod_tides.o mod_mean.o mod_archiv.o \
         mod_tsadvc.o mod_momtum.o mod_barotp.o mod_asselin.o mod_restart.o\
         mod_hycom.o

MODD =   mod_dimensions.o mod_xc.o mod_za.o mod_cb_arrays.o mod_pipe.o \
         mod_incupd.o \
         mod_floats.o mod_stokes.o mod_tides.o mod_mean.o mod_archiv.o \
         mod_tsadvc.o mod_momtum.o mod_barotp.o mod_asselin.o mod_restart.o\
         mod_hycom_dummy.o

OBJS =	                   bigrid.o blkdat.o  cnuity.o convec.o \
	diapfl.o dpthuv.o  dpudpv.o forfun.o  geopar.o hybgen.o \
	icloan.o inicon.o inigiss.o inikpp.o   inimy.o latbdy.o \
	matinv.o mxkprf.o  mxkrt.o  mxkrtm.o  mxpwp.o \
	overtn.o poflat.o  prtmsk.o  psmoo.o  \
	thermf.o trcupd.o  \
       machine.o  wtime.o machi_c.o  isnan.o s8gefs.o

hycom:	$(MODS) $(OBJS) hycom.o
	$(LD)  $(LDFLAGS) -o hycom  hycom.o $(MODS) $(OBJS) $(EXTRALIBS)

esmf:	$(MODS) $(OBJS)
	@echo "--- ESMF hycom component has been built ---"

hycom_cice:	$(MODS) $(OBJS) mod_OICPL.o hycom_cice.o
	$(LD)  $(LDFLAGS) -o hycom_cice \
                             hycom_cice.o mod_OICPL.o \
                             $(MODS) $(OBJS) \
                             ${CICE_DIR}/esmf/compile/*.o \
                             $(EXTRALIBS)

dummy_cice:	$(MODS) $(OBJS) mod_OICPL.o dummy_cice.o
	$(LD)  $(LDFLAGS) -o dummy_cice \
                             dummy_cice.o mod_OICPL.o \
                             $(MODD) $(OBJS) \
                             ${CICE_DIR}/esmf/compile/*.o \
                             $(EXTRALIBS)

hycom.o:        hycom.F90       mod_hycom.o
hycom_cice.o:   hycom_cice.F90  mod_hycom.o       mod_OICPL.o
dummy_cice.o:   dummy_cice.F90  mod_hycom_dummy.o mod_OICPL.o

bigrid.o:  bigrid.F90  mod_xc.o 
blkdat.o:  blkdat.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_incupd.o \
	                                                            mod_floats.o \
	                                                            mod_tides.o \
	                                                            mod_stokes.o
cnuity.o:  cnuity.F90  mod_xc.o mod_cb_arrays.o                     mod_pipe.o \
	                                                            mod_stokes.o
convec.o:  convec.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h
diapfl.o:  diapfl.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_stokes.o
dpthuv.o:  dpthuv.F90  mod_xc.o mod_cb_arrays.o
dpudpv.o:  dpudpv.F90  mod_xc.o 
forfun.o:  forfun.F90  mod_xc.o mod_cb_arrays.o            mod_za.o
geopar.o:  geopar.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h mod_za.o
hybgen.o:  hybgen.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o
icloan.o:  icloan.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h
inicon.o:  inicon.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
                                                                    mod_restart.o
inigiss.o: inigiss.F90 mod_xc.o mod_cb_arrays.o stmt_fns.h
inikpp.o:  inikpp.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h
inimy.o:   inimy.F90   mod_xc.o mod_cb_arrays.o stmt_fns.h
isnan.o:   isnan.F90
latbdy.o:  latbdy.F90  mod_xc.o mod_cb_arrays.o                     mod_tides.o
machine.o: machine.F90
machi_c.o: machi_c.c
matinv.o:  matinv.F90  mod_xc.o mod_cb_arrays.o
mxkprf.o:  mxkprf.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
	                                                            mod_stokes.o
mxkrt.o:   mxkrt.F90   mod_xc.o mod_cb_arrays.o stmt_fns.h
mxkrtm.o:  mxkrtm.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h
mxpwp.o:   mxpwp.F90   mod_xc.o mod_cb_arrays.o stmt_fns.h
overtn.o:  overtn.F90  mod_xc.o mod_cb_arrays.o
poflat.o:  poflat.F90
prtmsk.o:  prtmsk.F90
psmoo.o:   psmoo.F90   mod_xc.o 
s8gefs.o:  s8gefs.F90
thermf.o:  thermf.F90  mod_xc.o mod_cb_arrays.o stmt_fns.h
trcupd.o:  trcupd.F90  mod_xc.o mod_cb_arrays.o                     mod_pipe.o
wtime.o:   wtime.F90
mod_hycom.o: \
        mod_hycom.F90  mod_xc.o mod_cb_arrays.o            mod_za.o mod_pipe.o \
	                                                            mod_incupd.o \
	                                                            mod_mean.o \
	                                                            mod_floats.o \
	                                                            mod_momtum.o \
	                                                            mod_tsadvc.o \
	                                                            mod_barotp.o \
	                                                            mod_asselin.o \
                                                                    mod_restart.o \
	                                                            mod_stokes.o
mod_hycom_dummy.o: \
        mod_hycom_dummy.F90  mod_xc.o mod_cb_arrays.o      mod_za.o mod_pipe.o \
	                                                            mod_incupd.o \
	                                                            mod_mean.o \
	                                                            mod_floats.o \
	                                                            mod_momtum.o \
	                                                            mod_tsadvc.o \
	                                                            mod_barotp.o \
	                                                            mod_asselin.o \
                                                                    mod_restart.o \
	                                                            mod_stokes.o
mod_asselin.o: \
	mod_asselin.F90 mod_xc.o mod_cb_arrays.o stmt_fns.h         mod_pipe.o
mod_barotp.o: \
	mod_barotp.F90 mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
	                                                            mod_tides.o \
	                                                            mod_stokes.o
mod_cb_arrays.o: \
        mod_cb_arrays.F90 mod_dimensions.o
mod_momtum.o: \
	mod_momtum.F90 mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o \
	                                                            mod_tides.o \
	                                                            mod_stokes.o
mod_restart.o: \
	mod_restart.F90 mod_xc.o mod_cb_arrays.o           mod_za.o mod_tides.o
mod_tsadvc.o: \
	mod_tsadvc.F90 mod_xc.o mod_cb_arrays.o stmt_fns.h          mod_pipe.o
mod_incupd.o: \
        mod_incupd.F90 mod_xc.o mod_cb_arrays.o            mod_za.o
mod_floats.o: \
        mod_floats.F90 mod_xc.o mod_cb_arrays.o            mod_za.o mod_pipe.o \
	                                                            mod_stokes.o
mod_pipe.o: \
        mod_pipe.F90   mod_xc.o mod_cb_arrays.o                     mod_stokes.o
mod_stokes.o: \
        mod_stokes.F90 mod_xc.o mod_cb_arrays.o            mod_za.o
mod_tides.o: \
        mod_tides.F90  mod_xc.o mod_cb_arrays.o            mod_za.o
mod_mean.o: \
        mod_mean.F90   mod_xc.o mod_cb_arrays.o            mod_za.o mod_stokes.o
mod_archiv.o: \
        mod_archiv.F90 mod_xc.o mod_cb_arrays.o            mod_za.o mod_stokes.o

mod_dimensions.o:   mod_dimensions.F90 dimensions.h
mod_xc.o: mod_xc.F90  mod_dimensions.o mod_xc_sm.h mod_xc_mp.h
mod_za.o: mod_za.F90  mod_xc.o         mod_za_sm.h mod_za_mp.h mod_za_mp1.h mod_za_zt.h

mod_OICPL.o: mod_OICPL.F90
