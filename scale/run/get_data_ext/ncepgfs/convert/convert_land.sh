#!/bin/bash

tstart="$1"
tend="$2"
gribfile_prefix="$3"
outdir="$4"

#tstart='2016-08-10 06'
#tend='2016-08-10 06'
#gribfile_prefix='/data7/gylien/realtime/ncepgfs/2016081000/gfs'
#outdir='/data7/gylien/realtime/ncepgfs_grads/2016081000'

gribfile_dir="$(dirname ${gribfile_prefix})"
gribfile_base="$(basename ${gribfile_prefix})"
tint=21600

#----

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p ${outdir}


tinit=`echo $gribfile_dir | grep -o '[0-9]\{10\}'`
stimefgrads=$(date -ud "${tinit:0:4}-${tinit:4:2}-${tinit:6:2} ${tinit:8:2}" '+%H:%MZ%d%b%Y')

cat ${wkdir}/ctl/land.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/land.ctl

time="$tstart"
while (($(date -ud "$time" '+%s') <= $(date -ud "$tend" '+%s'))); do

  timef=$(date -ud "$time" '+%Y-%m-%d %H')
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  echo "[$timef]"

  gfile=${gribfile_prefix}.${timef2}
  ofile_land=${outdir}/land_${timef2}.grd

  tfile="tmp.grb2"
  rm -f ${tfile} ${ofile_land}
  #--Land data
  wgrib2 ${gfile} -match ":LAND:"  -match ":surface:"             -grib ${tfile}
  wgrib2 ${gfile} -match ":TMP:"   -match ":surface:"     -append -grib ${tfile}
#  wgrib2 ${gfile} -match ":TMP:"   -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":TSOIL:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":SOILW:" -match "below ground:" -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_land} -g2clib 0

  rm -f ${tfile}

time=$(date -ud "$tint second $time" '+%Y-%m-%d %H')
done
