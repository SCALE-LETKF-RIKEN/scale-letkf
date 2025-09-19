#!/bin/bash -l

wkdir="$( cd "$( dirname "$0" )" && pwd )"
datadir=$wkdir/../../../../../../external
rawdir=$datadir/ncepobs_gdas

cd ${wkdir}
now=`date -u +'%Y-%m-%d %H:%M:%S'`

if [ "$1" == "" ] ;then
  echo "input STIME"
  exit 1
fi

STIME=$1

GET_TIME=`date -u -d "${STIME:0:4}-${STIME:4:2}-${STIME:6:2} ${STIME:8:2}" +'%Y-%m-%d %H'`
YYYY=`date -u -d "$GET_TIME" +'%Y'`
MM=`date -u -d "$GET_TIME" +'%m'`
DD=`date -u -d "$GET_TIME" +'%d'`
HH=`date -u -d "$GET_TIME" +'%H'`
YYYYMMDDHH="$YYYY$MM$DD$HH"

echo "$now [TRY ] $YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log

mkdir -p ${rawdir}/$YYYYMMDDHH
cd ${rawdir}/$YYYYMMDDHH

allget=1
if [ ! -s "prepbufr.$YYYYMMDDHH" ]; then
  rm -f gdas.t${HH}z.prepbufr.nr
  echo "$now [GET ] start" >> ${wkdir}/get_ncep_obs.log
  wget --max-redirect 0 --cache=off --timeout=20 --no-check-certificate https://nomads.ncep.noaa.gov/pub/data/nccf/com/obsproc/prod/gdas.${YYYY}${MM}${DD}/gdas.t${HH}z.prepbufr.nr --referer="https://nomads.ncep.noaa.gov"
  if [ -s "gdas.t${HH}z.prepbufr.nr" ]; then
    mv -f gdas.t${HH}z.prepbufr.nr prepbufr.$YYYYMMDDHH
    now=`date -u +'%Y-%m-%d %H:%M:%S'`
    echo "$now [GET ] prepbufr.$YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log
  else
    allget=0
  fi
fi

cd ${wkdir}

if [ "$allget" -eq 1 ]; then
  now=`date -u +'%Y-%m-%d %H:%M:%S'`
  echo "$now [CONV] $YYYYMMDDHH: dec_prepbufr" >> ${wkdir}/get_ncep_obs.log

  CVOLUME=$(realpath $(pwd) | cut -d "/" -f 2) # current volume (e.g., /vol0X0Y or /vol000X)
  NUM_VOLUME=${CVOLUME:4:1} # get number of current volume 
  if [ "$NUM_VOLUME" = "0" ] ; then
    VOLUMES="/"${CVOLUME}
  else
    VOLUMES="/vol000${NUM_VOLUME}"
  fi
  if [ $VOLUMES != "/vol0004" ] ;then
    VOLUMES="${VOLUMES}:/vol0004" # spack
  fi

  SRCBIN=$wkdir/../../../obs/dec_prepbufr
  ln -sf $SRCBIN $wkdir/
  ln -sf ${rawdir}/$YYYYMMDDHH/prepbufr.$YYYYMMDDHH $wkdir/prepbufr.in

  if [ -z "$SCALE_NETCDF_C" ] || [ -z "$SCALE_NETCDF_F" ] || [ -z "$SCALE_PNETCDF" ] || [ -z "$SCALE_HDF" ] ; then
    echo "[Error] Export SCALE environmental parameters (e.g., SCALE_NETCDF_C)"
    exit 1
  fi
  export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib:${SCALE_NETCDF_C}/lib:${SCALE_NETCDF_F}/lib:${SCALE_PNETCDF}/lib:${SCALE_HDF}/lib:$LD_LIBRARY_PATH
  echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" > exec_dec_prepbufr.sh 
  echo "$wkdir/dec_prepbufr" >> exec_dec_prepbufr.sh

  echo "Submitting an interactive PJM JOB for dec_prepbufr"
  pjsub --interact -g ${GROUP} -L "node=1" --mpi "max-proc-per-node=1" -x PJM_LLIO_GFSCACHE=${VOLUMES} -L "elapse=600" -j -s --sparam "wait-time=300" "./exec_dec_prepbufr.sh"

  if [ -f $wkdir/fort.90 ]; then
    mkdir -p $datadir/ncepobs_gdas_letkf/${YYYYMMDDHH}
    mv fort.90 $datadir/ncepobs_gdas_letkf/${YYYYMMDDHH}/obs_${YYYYMMDDHH}0000.dat
    rm $wkdir/exec_dec_prepbufr.sh 
    now=`date -u +'%Y-%m-%d %H:%M:%S'`
    echo "$now [DONE] $YYYYMMDDHH" >> ${wkdir}/get_ncep_obs.log
    rm $wkdir/exec_dec_prepbufr.sh.*.out
    rm $wkdir/exec_dec_prepbufr.sh.*.stats
  else
    echo "dec_prepbufr stopped abnormally."
    exit 1
  fi

fi
