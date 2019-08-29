README.src/README.src.newmachine:

The Makefile sources ./config/$(ARCH)_$(TYPE) where ARCH defines exactly 
what machine architecture to target and TYPE is the parallelization 
strategy and precision (one, omp, mpi, ompi).  

The make process is automated by the script Make.csh, which should be used instead 
of directly invoking the make command.  It keys on the source code directory name 
which should end with _${TYPE}, where ${TYPE} is the parallelization type (one, omp,
mpi, ompi).  The script Make.csh should be edited to define ${ARCH} appropriately 
for the machine, several options are provided as examples.  In order for this to work, 
the file ./config/${ARCH}_${TYPE} must exist and must contain the machine-specific 
parts of Makefile (see README.src.config).

In addition the Make.csh script selects some CPP flags for compilations:
    OCN_SIG  -- Sigma-0 or Sigma-2
    OCN_EOS  -- Equation Of State (7,9,12,17 term)
    OCN_GLB  -- for global tripolar cases
    OCN_KAPP -- for centered thermobaricity correction
    OCN_MISC -- other CPP flags (-DSTOKES -DOCEANS2 etc...)

The executable is then created by the command:

    ./Make.csh >& Make.log

The script Make_global.csh is just Make.csh configured for the typical global tripole
case.  The script Make_cice.csh is for HYCOM+CICEv4 using the very old ESMFv4 coupler.
