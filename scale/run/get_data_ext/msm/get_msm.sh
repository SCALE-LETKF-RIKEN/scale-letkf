#!/bin/sh 

SRCDIR="http://database.rish.kyoto-u.ac.jp/arch/jmadata/data/gpv/original"

cd $(dirname $0)
wkdir=$(pwd)
scrpdir=$(cd $wkdir/../.. ; pwd)

MSM_DIR="$scrpdir/../../../../external/jma/msm"

TIMEf=$1

if [ "$1" == "" ]; then
  echo "input TIMEf"
  exit 1
fi 

if [ "${TIMEf:9:1}" == "" ] ;then
  echo "specify YYYYMMDDHH."
  exit 1
fi

TIMEf=${TIMEf:0:10}0000

yyyy=${TIMEf:0:4}
mm=${TIMEf:4:2}
dd=${TIMEf:6:2}
hh=${TIMEf:8:2}

outdir=$MSM_DIR/${TIMEf:0:10}

TIME="$yyyy-$mm-$dd $hh"

fhs="00-15"

gfile_atm="Z__C_RJTD_${TIMEf}_MSM_GPV_Rjp_L-pall_FH${fhs}_grib2.bin"
wget --cache=off $SRCDIR/$yyyy/$mm/$dd/$gfile_atm

gfile_sfc="Z__C_RJTD_${TIMEf}_MSM_GPV_Rjp_Lsurf_FH${fhs}_grib2.bin"
wget --cache=off $SRCDIR/$yyyy/$mm/$dd/$gfile_sfc

echo "convert..."
./convert.sh "$TIME" $gfile_atm $gfile_sfc $outdir &> convert.log

#rm $gfile_atm $gfile_sfc

echo "done."

