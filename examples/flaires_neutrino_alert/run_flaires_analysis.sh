#!/bin/bash

MYWORKDIR=$PWD
SCRIPTDIR=$0
SCRIPTPATH=${MYWORKDIR}${SCRIPTDIR:1}
EXECUTABLE=${SCRIPTPATH%.*}".py"
VENV="/data/user/gsommani/alertstack-venv"

CVFMS_PREFIX="/cvmfs/icecube.opensciencegrid.org"
CVMFS_DISTRIBUTION="py3-v4.4.0"
CVMFS_BASE=${CVFMS_PREFIX}/${CVMFS_DISTRIBUTION}
NU_SKYMAP_DIR="/data/ana/realtime/alert_catalog_v3/fits/"

eval $(${CVMFS_BASE}/setup.sh)
export NU_SKYMAP_DIR=${NU_SKYMAP_DIR}

source ${VENV}"/bin/activate"
python ${EXECUTABLE} $@

echo "Job complete!"