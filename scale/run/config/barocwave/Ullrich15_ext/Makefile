################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(shell cd ../../../.. ; pwd)
TESTDIR     = $(TOPDIR)/scale-rm/test


# user-defined source files
CODE_DIR    = .
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
TPROC       = 20

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
