#!/bin/bash
LIBPATH=/Users/sambaran/updated-nbody7/GPU2/
NBPATH=/Users/sambaran/updated-nbody7/GPU2/run_versions/
##########
cp -p $LIBPATH/gpunb.metallib $LIBPATH/gpupot.metallib .
$NBPATH/nbody7b.metal < input 1> run.out 2> err.out
#$NBPATH/nbody7b.metal < input.rst 1> run.out 2> err.out
