#!/bin/bash

cd /afs/cern.ch/user/f/folguera/workdir/EXO/CMSSW_8_0_11/src/MassResolution
eval `scramv1 runtime -sh`

makeResiduals_trkpt.py -i $1 -o $2


