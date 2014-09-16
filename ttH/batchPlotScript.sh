#!/bin/bash

hostname
cat /etc/redhat-release
uname -i

#export CURRENTDIR=`pwd` #get current directory on lxbatch machine

#setup atlas software
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
export DQ2_LOCAL_SITE_ID=DESY-HH_SCRATCHDISK 
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

#localSetupROOT
#asetup 17.3.12,slc5 --testarea /afs/desy.de/user/e/eckardt/athena/17.3.12/

cd /afs/desy.de/user/e/eckardt/Plotting/ttH/

asetup 19.1.1,64,slc6

python truth.py

