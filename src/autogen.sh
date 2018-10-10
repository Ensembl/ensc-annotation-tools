#!/bin/bash
###############################################################################
# Simple script to bootstrap and create the configure script, should be run
# when any control files have been changed (e.g. new source files added which
# led to changes to Makefile.am files) including:
#
#    configure.ac
#    Makefile.am
#    Makefile.in
#
############################################################

SCRIPT_NAME=$(basename $0)
INITIAL_DIR=$(pwd)
SCRIPT_DIR=$(dirname $0)
if ! echo $SCRIPT_DIR | egrep -q "(^)/" ; then
   BASE_DIR=$INITIAL_DIR/$SCRIPT_DIR
else
   BASE_DIR=$SCRIPT_DIR
fi


# set up version number, this is arcane and horrible and all
# autoconf's fault. See http://sources.redhat.com/automake/automake.html#Rebuilding
# and the stuff about AC_INIT
# NOTE that ensc_version.m4 is referenced from configure.ac
#
version_macro_file='ensc_annotation_tools_version.m4'
rm -f $version_macro_file

ENSC_VERSION=`git describe --always --abbrev=1`

echo 'dnl ensc_annotation_tools_version.m4 generated by autogen.sh script. Do not hand edit.'  > $version_macro_file
echo "m4_define([VERSION_NUMBER], [$ENSC_VERSION])" >> $version_macro_file
echo "EnsC-annotation-tools version is: $ENSC_VERSION"


# then for now we just run autoconf
autoreconf -fi -v || echo "autoreconf failed...."



