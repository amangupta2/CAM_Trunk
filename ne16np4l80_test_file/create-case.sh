#!/bin/bash

# Wall clock limit is max 12:00:00 = 12hrs for Cheyenne.
# Make sure to mention the project name.
export CIMEROOT=~/models/code/cslam_trunk/cime

export TEST_CASE_NAME=ne16np4l80k2h3_cslam_test

$CIMEROOT/scripts/create_newcase \
    --case $TEST_CASE_NAME \
    --compset FHS94 \
    --res ne16pg3_ne16pg3_mg37 \
    --machine cheyenne \
    --compiler intel \
    --mpilib openmpi \
    --walltime 12:00:00 \
    --input-dir /glade/scratch/agupta/model_output/CESM/inputdata \
    --run-unsupported  \
    --project UNYU0002

