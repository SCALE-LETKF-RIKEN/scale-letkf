#!/bin/bash
#PJM -L "rscgrp=small"
#PJM -L "node=1"
#PJM -L "elapse=1800"
#PJM --mpi "max-proc-per-node=4"
#PJM -s
#PJM -j

export PARALLEL=12
export OMP_NUM_THREADS=12
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD



PYTHON_DIR=/home/apps/oss/TensorFlow-2.1.0
PYTHON=python3
WORK_PATH=`pwd`
module switch lang/tcsds-1.2.33
export PATH=${PYTHON_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${PYTHON_DIR}/lib:$LD_LIBRARY_PATH
ldd `dirname \`which python3\``/../lib/python3.8/site-packages/tensorflow_core/libtensorflow_framework.so.2
ldd `dirname \`which python3\``/../lib/python3.8/site-packages/horovod/tensorflow/mpi_lib.cpython-38-aarch64-linux-gnu.so
#PROC=4
#MPI="mpirun -np $PROC --mca mpi_print_stats 1"

#mpiexec -std-proc log/NOUT -n 16 ./scale-rm_single run_mean.conf


$PYTHON ./init_perturb.py init_bdy_20200825-180000.000


