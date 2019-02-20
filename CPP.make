#
# --- HYCOM 2.0 surce code makefile 
#
# --- Tunable parameters in ../config/$(ARCH)_$(TYPE)
#

.SUFFIXES: 
.SUFFIXES: .c .F .f .o

.F:
	@echo "Must have an explicit rule for" $*
.f:
	@echo "Must have an explicit rule for" $*
.c:
	@echo "Must have an explicit rule for" $*

include ../config/$(ARCH)_$(TYPE)

.F.f:
	$(RM) $<.f $<.C
	sed -e 's? */// *?/ / /?g' -e 's? *// *?/ /?g' $< >  $<.C
	$(CPP) $(CPPFLAGS) $<.C | sed -e '/^ *$$/d' > $<.f
#	-\mv $<.f $*.f
	$(RM) $<.C

default: hycom.f

hycom.f:        hycom.F
hycom_cice.f:   hycom_cice.F
dummy_cice.f:   dummy_cice.F

barotp.f:  barotp.F
blkdat.f:  blkdat.F
cnuity.f:  cnuity.F
convec.f:  convec.F
diapfl.f:  diapfl.F
dpthuv.f:  dpthuv.F
dpudpv.f:  dpudpv.F
geopar.f:  geopar.F
hybgen.f:  hybgen.F
icloan.f:  icloan.F
inicon.f:  inicon.F
isnan.f:   isnan.F
latbdy.f:  latbdy.F
machine.f: machine.F
mxkprf.f:  mxkprf.F
mxkrt.f:   mxkrt.F
mxkrtm.f:  mxkrtm.F
mxpwp.f:   mxpwp.F
overtn.f:  overtn.F
thermf.f:  thermf.F
trcupd.f:  trcupd.F
wtime.f:   wtime.F
mod_hycom.f:       mod_hycom.F
mod_hycom_dummy.f: mod_hycom_dummy.F
mod_cb_arrays.f:   mod_cb_arrays.F
mod_dimensions.f:  mod_dimensions.F
mod_momtum.f: mod_momtum.F
mod_tsadvc.f: mod_tsadvc.F
mod_incupd.f: mod_incupd.F
mod_floats.f: mod_floats.F
mod_pipe.f:   mod_pipe.F
mod_stokes.f: mod_stokes.F
mod_tides.f:  mod_tides.F
mod_mean.f:   mod_mean.F
mod_archiv.f: mod_archiv.F
mod_xc.f:     mod_xc.F
mod_za.f:     mod_za.F
mod_OICPL.f:  mod_OICPL.F
