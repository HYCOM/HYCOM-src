#!/bin/csh
#
set echo
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
# --- make hycom with TYPE from this directory's name (src_*_$TYPE).
# --- assumes dimensions.h is correct for $TYPE.
#
# --- set ARCH to the correct value for this machine.
# --- ARCH that start with A are for ARCTIC patch regions
#
#module swap compiler compiler/intel/12.1.3
#module swap mpi      mpi/intel/ibmpe
#module list
#setenv ARCH Aintelsse-pe-sm-relo
#
#module swap compiler compiler/intel/12.1.3
#module swap mpi      mpi/intel/impi/4.1.3
#module list
#setenv ARCH Aintelsse-impi-sm-SD-relo
#setenv ARCH Aintelsse-impi-sm-relo
#
unset echo
#module switch PrgEnv-cray PrgEnv-intel
module unload cray-libsci
module switch intel intel/15.0.0.090
module switch cray-mpich cray-mpich/7.0.3
module load   craype-hugepages2M
module list
set echo
setenv ARCH xc40-intel-relo
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
echo "ARCH = " $ARCH "  TYPE = " $TYPE
#
if (! -e ./config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
endif

# CPP flags for compilations
# Equation Of State
setenv OCN_SIG  -DEOS_SIG2 ## Sigma-2
#setenv OCN_SIG -DEOS_SIG0 ## Sigma-0

#setenv OCN_EOS -DEOS_7T  ## EOS  7-term
setenv OCN_EOS -DEOS_9T  ## EOS  9-term
#setenv OCN_EOS -DEOS_12T ## EOS 12-term
#setenv OCN_EOS -DEOS_17T ## EOS 17-term

# Optional CPP flags
# Global or regional
setenv OCN_GLB -DARCTIC ## global tripolar simulation
#setenv OCN_GLB ""

# Thermobaricity correction centered
#setenv OCN_KAPP -DKAPPAF_CENTERED
setenv OCN_KAPP ""

# Miscellaneous CPP flags (-DSTOKES -DOCEANS2 etc...)
# -DSTOKES  : Stokes drift
# -DOCEANS2 : master and slave HYCOM in same executable
setenv OCN_MISC ""

# CPP_EXTRAS
setenv CPP_EXTRAS "${OCN_SIG} ${OCN_EOS} ${OCN_GLB} ${OCN_KAPP} ${OCN_MISC}"

#
# --- some machines require gmake
#
#gmake ARCH=$ARCH TYPE=$TYPE hycom
make hycom ARCH=$ARCH TYPE=$TYPE
#
if ( $ARCH == "Asp5" || $ARCH == "sp5") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
if ( $ARCH == "Asp6" || $ARCH == "sp6") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
if ( $ARCH == "Asp6-nofl" || $ARCH == "sp6-nofl") then
  ldedit -bdatapsize=64K -bstackpsize=64K hycom
endif
