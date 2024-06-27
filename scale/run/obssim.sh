#!/bin/bash 
# This script launchs obssim
#===============================================================================
# Configuration

. ./config.main || exit $?
. ./config.fcst || exit $?

. src/func_datetime.sh || exit $?
. src/func_distribute.sh || exit $?
. src/func_util.sh || exit $?

. src/func_common_static.sh || exit $?
#. src/func_${job}_static.sh || exit $?

echo "[$(datetime_now)] Start obssim $@"

#===============================================================================
# Setting
TSTART=8
TEND=13
TINTERVAL=${FCSTOUT}

#-------------------------------------------------------------------------------
# Create and clean up the temporary directory
TMP=${TMP}_obssim
safe_init_tmpdir $TMP || exit $?

# copy config files
mkdir -p $TMP/config
cp $SCRP_DIR/config.nml.* ${TMP}/config/
cp $SCRP_DIR/config.[c,f,r]* ${TMP}/config/
cp $SCRP_DIR/config.main.${PRESET} ${TMP}/config/

cp $SCRP_DIR/../obs/obssim ${TMP}/

#-------------------------------------------------------------------------------
# copy data files
mkdir -p ${TMP}/dat
cp -r ${SCALEDIR}/data/rad ${TMP}/dat/

if [ -n "${RTTOV_COEF_PATH}" ] && [ -n "${RTTOV_COEF_FILE}" ] ; then
  mkdir -p ${TMP}/dat/rttov
  cp -r ${RTTOV_COEF_PATH}/${RTTOV_COEF_FILE} ${TMP}/dat/rttov/
  cp -r ${RTTOV_COEF_PATH}/${RTTOV_COEF_FILE_CLD} ${TMP}/dat/rttov/
fi

#-------------------------------------------------------------------------------
# create directory for log files
logd=$OUTDIR/$STIME/log/obssim
mkdir -p $logd

#-------------------------------------------------------------------------------
# prepare config file
config=$TMP/config/obssim.conf
cat << EOF > $config
&PARAM_LOG
 LOG_LEVEL = 3,
/

&PARAM_ENSEMBLE
  MEMBER = ${MEMBER},
  MEMBER_RUN = 1,
  DET_RUN = .false.,
  EFSO_RUN = .false.,
/

&PARAM_PROCESS
  PPN = ${PPN},
  !--MEM_NODES--
  MEM_NODES = $((SCALE_NP_X*SCALE_NP_Y/PPN)),
  !--NUM_DOMAIN--
  NUM_DOMAIN = 1,
  !--PRC_DOMAINS--
  PRC_DOMAINS = $((SCALE_NP_X*SCALE_NP_Y)), 
/

&PARAM_OBSSIM
  OBSSIM_IN_TYPE = 'history',
  OBSSIM_HISTORY_IN_BASENAME = '${OUTDIR}/${STIME}/fcst/mean/history',
  OBSSIM_TIME_START = ${TSTART},
  OBSSIM_TIME_END   = ${TEND},
  OBSSIM_TIME_INTERVAL_SEC = ${TINTERVAL},
  OBSSIM_HIM = T,
  OBSSIM_NC_OUT_BASENAME = '${OUTDIR}/${STIME}/fcst/mean/himawari',
/
EOF

    cat $TMP/config/config.nml.letkf | \
        sed -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
            -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
            -e "/!--SLOT_START--/a SLOT_START = $slot_s," \
            -e "/!--SLOT_END--/a SLOT_END = $slot_e," \
            -e "/!--SLOT_BASE--/a SLOT_BASE = $slot_b," \
            -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = ${LTIMESLOT}.D0," \
            -e "/!--LETKF_TOPOGRAPHY_IN_BASENAME--/a LETKF_TOPOGRAPHY_IN_BASENAME = \"${TOPO_PATH}/topo/topo\"," \
            -e "/!--HIM_RTTOV_THREADS--/a HIM_RTTOV_THREADS = ${THREADS}," \
            -e "/!--HIM_NOWDATE--/a HIM_NOWDATE = ${STIME:0:4}, ${STIME:4:2}, ${STIME:6:2}, ${STIME:8:2}, ${STIME:10:2}, ${STIME:12:2}," \
            -e "/!--RTTOV_COEF_PATH--/a RTTOV_COEF_PATH = \"${TMP}/dat/rttov\"," \
            -e "/!--RTTOV_COEF_FILE--/a RTTOV_COEF_FILE = \"${RTTOV_COEF_FILE##*/}\"," \
            -e "/!--RTTOV_COEF_FILE_CLD--/a RTTOV_COEF_FILE_CLD = \"${RTTOV_COEF_FILE_CLD##*/}\"," \
            -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME =  \"${OUTDIR[$d]}/$STIME/log/obssim/obssim_LOG\"," \
    >> ${config}

    cat $TMP/config/config.nml.scale | \
        sed -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = ${STIME:0:4}, ${STIME:4:2}, ${STIME:6:2}, ${STIME:8:2}, ${STIME:10:2}, ${STIME:12:2}," \
            -e "/!--TIME_DURATION--/a TIME_DURATION = ${FCSTLEN}.d0," \
            -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMP}/dat/rad/cira.nc\"," \
            -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMP}/dat/rad/MIPAS\"," \
    >> ${config}

#-------------------------------------------------------------------------------
# Creat a job script

jobscrp="$TMP/obssim_job.sh"
TIME_LIMIT="00:30:00"

echo "[$(datetime_now)] Create a job script '$jobscrp'"

# FUGAKU
NNODES_USE=$((SCALE_NP_X*SCALE_NP_Y/PPN))
if [ "$PRESET" = 'FUGAKU' ]; then

  if (( NNODES_USE > 384 )) ; then
    RSCGRP="large"
  else
    RSCGRP="small"
  fi

  TPROC=$((NNODES_USE*PPN))

  CVOLUME=$(pwd | cut -d "/" -f 2) # current volume (e.g., /vol0X0Y or /vol000X)
  NUM_VOLUME=${CVOLUME:4:1} # get number of current volume 

  if [ "$NUM_VOLUME" = "0" ] ; then
    VOLUMES="/"${CVOLUME}
  else
    VOLUMES="/vol000${NUM_VOLUME}"
  fi

  if [ $VOLUMES != "/vol0004" ] ;then
    VOLUMES="${VOLUMES}:/vol0004" # spack
  fi

cat > $jobscrp << EOF
#!/bin/sh 
#
#PJM -g ${GROUP} 
#PJM -x PJM_LLIO_GFSCACHE=${VOLUMES}
#PJM -L "rscgrp=${RSCGRP}"
#PJM -L "node=$(((TPROC+PPN-1)/PPN))"
#PJM -L "elapse=${TIME_LIMIT}"
#PJM --mpi "max-proc-per-node=${PPN}"
#PJM -j
#PJM -s
EOF

  if (( BDY_TMP == 1 )) && (( BDY_ENS == 1 )); then
    echo "#PJM --llio localtmp-size=${LLIO_TMP_SIZE}Gi" >> $jobscrp
  fi

cat >> $jobscrp << EOF
#
#
export PARALLEL=${THREADS}
export OMP_NUM_THREADS=\${PARALLEL}
export FORT90L=-Wl,-T
export PLE_MPI_STD_EMPTYFILE=off
export OMP_WAIT_POLICY=active
export FLIB_BARRIER=HARD

EOF

  if (( USE_LLIO_BIN == 1 )); then
    for i in $(seq $nsteps) ; do
      echo "llio_transfer ${stepexecbin[$i]}" >> $jobscrp 
    done
    echo "" >> $jobscrp
  fi

  if (( USE_LLIO_DAT == 1 )); then
    echo "/home/system/tool/dir_transfer -l ./ ${TMP}/dat" >> $jobscrp
    echo "" >> $jobscrp
  fi

  if (( USE_SPACK > 0 )); then

cat << EOF >>  $jobscrp 
SPACK_FJVER=${SPACK_FJVER}
. /vol0004/apps/oss/spack/share/spack/setup-env.sh
spack load netcdf-c%fj@\${SPACK_FJVER}
spack load netcdf-fortran%fj@\${SPACK_FJVER}
spack load parallel-netcdf%fj@\${SPACK_FJVER}

export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:\$LD_LIBRARY_PATH

EOF

  else

    if [ -z "$SCALE_NETCDF_C" ] || [ -z "$SCALE_NETCDF_F" ] || [ -z "$SCALE_PNETCDF" ] || [ -z "$SCALE_HDF" ] ; then
      echo "[Error] Export SCALE environmental parameters (e.g., SCALE_NETCDF_C)"
      exit 1
    fi

cat << EOF >>  $jobscrp 
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:\$LD_LIBRARY_PATH

EOF

  fi

cat << EOF >>  $jobscrp 
mpiexec -std-proc ${logd}/%/200r/NOUT -n $((NNODES_USE*PPN)) $TMP/obssim ${config} || exit $?
EOF

  if (( USE_LLIO_BIN == 1 )); then
    for i in $(seq $nsteps) ; do
      echo "llio_transfer --purge ${stepexecbin[$i]}" >> $jobscrp 
    done
  fi

  if (( USE_LLIO_DAT == 1 )); then
    echo "/home/system/tool/dir_transfer -p -l ./ ${TMP}/dat" >> $jobscrp
    echo "" >> $jobscrp
  fi

  echo "[$(datetime_now)] Run ${job} job on PJM"
  echo

  job_submit_PJM $jobscrp
  echo
  
  job_end_check_PJM $jobid
  res=$?
fi

#-------------------------------------------------------------------------------

echo "[$(datetime_now)] Finish obssim $@"
