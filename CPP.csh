#!/bin/csh
#
set echo
#
# --- run HYCOM through CPP.
#
#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./CPP.com >& CPP.log
#
# --- Build a source code set for this ARCH and TYPE only.
# --- assumes dimensions.h is correct for $TYPE.
#
#setenv ARCH intel-pgi-relo
setenv ARCH intelsse-impi-sm-relo
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
else
# foreach f ( *.F )
# foreach f ( barotp.F blkdat.F cnuity.F mod_archiv.F mod_floats.F mod_hycom.F mod_mean.F mod_momtum.F mod_pipe.F mxkprf.F thermf.F )
  foreach f ( barotp.F mod_momtum.F mxkprf.F )
    make -f CPP.make $f:r.f ARCH=$ARCH TYPE=$TYPE
  end
endif
