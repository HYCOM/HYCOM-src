#
set echo
#
setenv C ~abozec/hycom/GOMd0.08/src_2.2.110_relo_mpi`
#
echo "*****     *****     *****     *****     *****     *****     *****"
foreach f ( Makefile Makefile.NUOPC Make.csh Make_cice.csh )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -bw $f $C/$f
end
foreach f ( ALT_CODE/*.h *.h *.c )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f $C/$f
end
#allow for possible switch from .f90 to .F90 or .F90 to .f90
foreach f ( *.f90 *.F90 )
  echo "*****     *****     *****     *****     *****     *****     *****"
  diff -ibw $f $C/$f:r.[Ff]90
end
