#!/bin/sh

#mem=$1

filebase="init_bdy_20220101-000000.000"

#for cmem in `seq -f %04g $mem` ; do 
for cmem in `seq -f %04g 11 50` ; do 
  mkdir -p $cmem
  cp -r mean_orig/* $cmem/
  python3 ./init_perturb.py $cmem/$filebase  
done
