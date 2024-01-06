#!/bin/bash

module load miniconda/3-22.11 > /dev/null 2>&1
module load singularity/3.8.3
module load fsl/6.0.3

fsLicense="/appl/freesurfer-7.1.1/license.txt"

scriptPath=$(readlink -f "$0")
scriptDir=$(dirname "${scriptPath}")
# Repo base dir under which we find bin/, containers/ and scripts/
repoDir=${scriptDir%/bin}

function usage() {
  echo "Usage:
  $0 -v synb0_version -L logdir [-n num_cores=1] -- [options to bidsSynB0.py]

  $0 -h for help
  "
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

function help() {
cat << HELP
  `usage`

  This is a wrapper script to submit processing to bsub.

  bidsSynBOLD options are below. The following args are set automatically by this wrapper:

     -c / --container (set by -v option to this wrapper)
     -n / --num-threads (set by -n option to this wrapper)
     --fs-license-file (hard-coded to ${fsLicense})

  `conda run -p /project/ftdc_pipeline/ftdc-picsl/miniconda/envs/ftdc-picsl-cp311 ${repoDir}/scripts/bidsSynBOLD.py -h`

HELP

}

numThreads=1

while getopts "L:n:v:h" opt; do
  case $opt in
    L) logDir=$OPTARG;;
    n) numThreads=$OPTARG;;
    v) synBOLDVersion=$OPTARG;;
    h) help; exit 1;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;
  esac
done

shift $((OPTIND-1))

date=`date +%Y%m%d`

# Makes python output unbuffered, so we can tail the log file and see progress
# and errors in order
export PYTHONUNBUFFERED=1

if [[ -z ${logDir} ]]; then
    echo "Please specify a log directory with -L"
    exit 1
fi

if [[ ! -d ${logDir} ]]; then
    mkdir -p ${logDir}
fi

bsub -n $numThreads -cwd . -o "${logDir}/synBOLD_${date}_%J.txt" \
    conda run -p /project/ftdc_pipeline/ftdc-picsl/miniconda/envs/ftdc-picsl-cp311 ${repoDir}/scripts/bidsSynBOLD.py \
      --container ${repoDir}/containers/synbold-disco-${synBOLDVersion}.sif \
      --num-threads $numThreads \
      --fs-license-file ${fsLicense} \
      $*