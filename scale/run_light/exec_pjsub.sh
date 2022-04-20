#!/bin/sh 
#
#PJM -g hp150019
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM -L "rscgrp=small"
#PJM -L "node=12"
#PJM -L "elapse=00:20:00"
#PJM --mpi "max-proc-per-node=4"
#PJM -j
#PJM -s
#
export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:/vol0004/apps/oss/spack-v0.17.0/opt/spack/linux-rhel8-a64fx/fj-4.7.0/netcdf-c-4.8.1-p4tpurxsxgclg6lfb3m7gwsmjqkuh55c/lib:/vol0004/apps/oss/spack-v0.17.0/opt/spack/linux-rhel8-a64fx/fj-4.7.0/netcdf-fortran-4.5.3-6gjvknqptlidwabdy7a5nyfxe7eap47b/lib:/vol0004/apps/oss/spack-v0.17.0/opt/spack/linux-rhel8-a64fx/fj-4.7.0/parallel-netcdf-1.12.2-v3xpd4kbmrw5dtedaqf2veexcdvaqmhw/lib:/vol0004/apps/oss/spack-v0.17.0/opt/spack/linux-rhel8-a64fx/fj-4.7.0/hdf5-1.10.7-e55vbz2u6w6tzq74f6tp6t36b2uza6dh/lib:$LD_LIBRARY_PATH

echo "scale-rm_init_ens"
mpiexec -std-proc log/scale_init/NOUT -n 48 ./scale-rm_init_ens config/scale-rm_init_ens_20220101000000.conf
#echo "scale-rm_ens"
#mpiexec -std-proc log/scale/NOUT -n 48 ./scale-rm_ens config/scale-rm_ens_20220101000000.conf
#echo "letkf"
#mpiexec -std-proc log/letkf/NOUT -n 48 ./letkf config/letkf_20220101060000.conf
echo "done."
