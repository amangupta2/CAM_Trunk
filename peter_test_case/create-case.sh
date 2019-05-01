#!/bin/bash

# Wall clock limit is max 12:00:00 = 12hrs for Cheyenne.
# Make sure to mention the project name.
export CIMEROOT=~/models/code/git_trunk/cime

export TEST_CASE_NAME=ne16np4l40k2h3_aoa_relax_tss_analytic_ueq

$CIMEROOT/scripts/create_newcase \
    --case $TEST_CASE_NAME \
    --compset FHS94 \
    --res ne16_ne16_mg17 \
    --machine cheyenne \
    --compiler intel \
    --mpilib openmpi \
    --walltime 12:00:00 \
    --input-dir /glade/scratch/agupta/model_output/CESM/inputdata \
    --run-unsupported  \
    --project UNYU0002 

# --q verylong (for spectral element)
# --res ne16_ne16_mg17 \ OR ne30_ne30_mg17
# OR ne120_ne120_mg16
# OR --res f19_f19 \ f09_f09
# OR --res T42_T42 \
# To mention queue :  --queue cpu48 \
