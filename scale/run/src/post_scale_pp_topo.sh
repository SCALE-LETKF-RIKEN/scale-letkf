#!/bin/bash
#===============================================================================
#
#  Script to post-process the SCALE model outputs.
#  November 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.all

if (($# < 5)); then
  cat >&2 << EOF

[post_scale_pp_topo.sh] Post-process the SCALE model outputs.

Usage: $0 MYRANK STIME TMPDIR

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  TMPDIR   Temporary directory to run the model

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
TMPDIR="$1"

#===============================================================================

if ((LOG_OPT <= 2)); then
  mkdir -p $TMPOUT/${STIME}/log/scale_topo
  if [ -f "$TMPDIR/pp_LOG${SCALE_LOG_SFX}" ]; then
    mv -f $TMPDIR/pp_LOG${SCALE_LOG_SFX} $TMPOUT/${STIME}/log/scale_topo/pp_LOG${SCALE_LOG_SFX}
  fi
fi

#===============================================================================

exit 0
