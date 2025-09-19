#!/bin/bash -l

cd $(dirname $0)
wkdir=$(pwd)
datadir=$wkdir/../../../../../../external/ncepgfs_grads
rawdir=$wkdir/../../../../../../external/ncepgfs

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

### Get archive
day_realtime=7
day_archive_newest=14
day_archive_oldest=$((175+day_archive_newest))
nows=$(date -ud "$now" +%s)
gets=$(date -ud "$GET_TIME" +%s)

if (( $((nows-gets)) < $((day_realtime*86400)) )) ;then
  get_from_archive=0
elif (( $((nows-gets)) <= $((day_archive_oldest*86400)) && $((nows-gets)) >= $((day_archive_newest*86400)) )) ;then
  get_from_archive=1
else
  echo "$GET_TIME data not found in the database." 
  exit 
fi

echo "$now [TRY ] $YYYYMMDDHH" >> ${wkdir}/get_ncep_gfs.log

mkdir -p ${rawdir}/$YYYYMMDDHH
cd ${rawdir}/$YYYYMMDDHH

allget=1

t=0
while ((t <= 18)); do

  tf=`printf '%03d' $t`
  TIME_fcst=`date -ud "${t} hour $YYYY-$MM-$DD $HH" +'%Y-%m-%d %H:%M:%S'`
  YYYYMMDDHHMMSS_fcst=`date -ud "$TIME_fcst" +'%Y%m%d%H%M%S'`

  if [ ! -s "gfs.$YYYYMMDDHHMMSS_fcst" ]; then
    rm -f gfs.t${HH}z.pgrb2f${tf}
    if ((get_from_archive==1));then
      rawsrc="https://www.ncei.noaa.gov/data/global-forecast-system/access/grid-004-0.5-degree/forecast/${YYYY}${MM}/${YYYY}${MM}${DD}/gfs_4_${YYYY}${MM}${DD}_${HH}00_${tf}.grb2" ### half-year archive
    else
      rawsrc="https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${YYYY}${MM}${DD}/${HH}/atmos/gfs.t${HH}z.pgrb2.0p50.f${tf}" ### real time
    fi
    rawfile=$(basename $rawsrc)

    itry=0
    while [ $itry -le 10 ] ; do
      wget -N --cache=off --timeout=20 --no-check-certificate $rawsrc 
      [ $? == 0 ] && break
      itry=`expr $itry + 1`
      echo "retry " $itry " ..."
    done

    if [ -s $rawfile ]; then
      mv -f $rawfile gfs.$YYYYMMDDHHMMSS_fcst
      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [GET ] $YYYYMMDDHH -> gfs.$YYYYMMDDHHMMSS_fcst" >> ${wkdir}/get_ncep_gfs.log
    else
      allget=0
    fi
  fi

  if ((allget == 1)); then
      now=`date -u +'%Y-%m-%d %H:%M:%S'`
      echo "$now [CONV] $YYYYMMDDHH -> gfs.$YYYYMMDDHHMMSS_fcst - GrADS for SCALE" >> ${wkdir}/get_ncep_gfs.log
      bash $wkdir/convert/convert.sh "$TIME_fcst" "$TIME_fcst" "$rawdir/${YYYYMMDDHH}/gfs" "$datadir/${YYYYMMDDHH}" \
       > ${wkdir}/convert_grads_scale.log 2>&1
  fi

t=$((t+3))
done

echo "done." 

