#!/bin/bash

cd $(dirname $0)
wkdir=$(pwd)

#wgrib2dir="/share/hp150019/bin"
#export PATH=$wgrib2dir:$PATH

tstart="$1"
gfile_atm="$2"
gfile_sfc="$3"
outdir="$4"

ofile_sfc_base=${outdir}/sfc
ofile_atm_base=${outdir}/atm

#----
stimefgrads=$(date -ud "$tstart" '+%H:%MZ%d%b%Y')

echo cat ${wkdir}/ctl/atm.ctl sed "s/<STIME>/${stimefgrads}/g" ${outdir}/atm.ctl

cat ${wkdir}/ctl/atm.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/atm.ctl
cat ${wkdir}/ctl/sfc.ctl | sed "s/<STIME>/${stimefgrads}/g" > ${outdir}/sfc.ctl

wkdir="$( cd "$( dirname "$0" )" && pwd )"
cd ${wkdir}

mkdir -p ${outdir}

time="$tstart"

  timef=$(date -ud "$time" '+%Y-%m-%d %H')
  echo "[$timef]"

  gfile=${gfile_sfc}
  
  tstamp=$(date -ud "$time" '+%Y%m%d%H%M%S')
  tfile="tmp.grb2.${tstamp}"

  #--Surface data
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  ofile_sfc=${ofile_sfc_base}_${timef2}.grd
  rm -f ${tfile} ${ofile_sfc}
  wgrib2 ${gfile} -match ":anl:" -match ":PRMSL:" -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":PRES:"  -match ":surface:"           -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":UGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":VGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":TMP:"   -match ":1.5 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":RH:"    -match ":1.5 m above ground:"  -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_sfc} -g2clib 0
  rm $tfile

  for ft in 3 6 9 12 15 ; do
  timef2=$(date -ud "$ft hour $time" '+%Y%m%d%H%M%S')
  ofile_sfc=${ofile_sfc_base}_${timef2}.grd
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":PRMSL:" -match ":mean sea level:"            -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":PRES:"  -match ":surface:"           -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":UGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":VGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":TMP:"   -match ":1.5 m above ground:"  -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":RH:"    -match ":1.5 m above ground:"  -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_sfc} -g2clib 0
  rm $tfile
  done

  gfile=${gfile_atm}
  #--Upper data
  timef2=$(date -ud "$time" '+%Y%m%d%H%M%S')
  ofile_atm=${ofile_atm_base}_${timef2}.grd
  wgrib2 ${gfile} -match ":anl:" -match ":HGT:"  | wgrib2 -i ${gfile} -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":UGRD:" | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":VGRD:" | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":TMP:"  | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":anl:" -match ":RH:"   | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_atm} -g2clib 0
  rm $tfile

  for ft in 3 6 9 12 15 ; do
  timef2=$(date -ud "$ft hour $time" '+%Y%m%d%H%M%S')
  ofile_atm=${ofile_atm_base}_${timef2}.grd
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":HGT:"  | wgrib2 -i ${gfile} -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":UGRD:" | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":VGRD:" | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":TMP:"  | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${gfile} -match ":$ft hour fcst:" -match ":RH:"   | wgrib2 -i ${gfile} -append -grib ${tfile}
  wgrib2 ${tfile} -no_header -ieee ${ofile_atm} -g2clib 0
  rm $tfile
  done

