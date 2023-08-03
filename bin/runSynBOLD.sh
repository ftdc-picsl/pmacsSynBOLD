#!/bin/bash -e

module load singularity/3.8.3

cleanup=1
fsLicense="/appl/freesurfer-7.1.1/license.txt"

function usage() {
  echo "Usage:
  $0 [-h] [-B src:dest,...,src:dest] [-c 1/0] \\
    -v synBOLDVersion -i /path/to/input -o /path/to/output -- [extra args]
"
}

function help() {
    usage
  echo "Wrapper for SyNBOLD distortion correction.

Use absolute paths for input and output, as these have to be mounted in the container.

Using the options below, specify paths on the local file system. These will be bound automatically
to locations inside the container. If needed, you can add extra mount points with '-B'.

Any args after the '--' should reference paths within the container. Currently, supported extra
arg are:
  --notopup prevents topup running automatically after the distortion correction.
  --stripped input T1w is skull-stripped.

Required args:

  -i /path/to/input
     Input directory on the local file system. Will be bound to /INPUTS inside the container. This must contain
     T1.nii.gz, b0.nii.gz, acqparams.txt. It may also contain a brain mask called T1_mask.nii.gz, but this is
     not required.

  -o /path/to/outputDir
     Output directory on the local files system. Will be bound to /OUTPUTS inside the container.

  -v version
     version. The script will look for containers/synBOLD-[version].sif.


Options:

  -B src:dest[,src:dest,...,src:dest]
     Use this to add mount points to bind inside the container, that aren't handled by other options.
     'src' is an absolute path on the local file system and 'dest' is an absolute path inside the container.
     Several bind points are always defined inside the container including \$HOME, \$PWD (where script is
     executed from), and /tmp (more on this below).

  -c 1/0
     Cleanup the working dir after running the prep (default = $cleanup).

  -h
     Prints this help message.


*** Hard-coded configuration ***

The FreeSurfer license file is sourced from ${fsLicense} .

The singularity module sets the singularity temp dir to be on /scratch. To avoid conflicts with other jobs,
the script makes a temp dir specifically for this job under /scratch. By default it is removed after
the code finishes, but this can be disabled with '-c 0'.

The total number of threads for OMP and ITK is set to \$LSB_DJOB_NUMPROC, which is the number of slots requested
with the -n argument to bsub.


*** Additional args ***

From the SynBOLD README:

--no_topup

Skip the application of FSL's topup susceptibility correction. As a default, we run topup for you.

--motion_corrected

Lets the pipeline know that supplied distorted bold image has already been motion corrected. As a default, we motion correct the distorted image.

--skull_stripped

Lets the container know the supplied T1 has already been skull-stripped. As a default, we assume it is not skull stripped. Please note this feature requires a well-stripped T1 as stripping artifacts can affect performance.

--no_smoothing

By default, Gaussian kernel with standard deviation 1.15 mm were applied to BOLD_d to match the smoothness of the synthetic image. This flag alters the default behavior and skip the smoothing step.

--custom_cnf

User can configurate their own cnf file for topup.

--total_readout_time

User can specify readout time by adding a desired number after (with a space) this flag. Any argument after this flag will be saved as the new total_readout_time. By default, the readout time of distortion BOLD image is 1.

*** Citation ***

Yu T, Cai LY, Morgan VL, Goodale SE, Englot DJ, Chang CE, Landman BA, Schilling KG. SynBOLD-DisCo: Synthetic BOLD images for distortion correction of fMRI without additional calibration scans. Proc SPIE Int Soc Opt Eng. 2023 Feb;12464:1246417. doi: 10.1117/12.2653647. Epub 2023 Apr 3. PMID: 37465092; PMCID: PMC10353777.

"
}

if [[ $# -eq 0 ]]; then
  usage
  exit 1
fi

scriptPath=$(readlink -f "$0")
scriptDir=$(dirname "${scriptPath}")
# Repo base dir under which we find bin/ and containers/
repoDir=${scriptDir%/bin}

userBindPoints=""
containerVersion=""

while getopts "B:c:f:i:m:o:t:v:h" opt; do
  case $opt in
    B) userBindPoints=$OPTARG;;
    c) cleanup=$OPTARG;;
    h) help; exit 1;;
    i) inputDir=$OPTARG;;
    o) outputDir=$OPTARG;;
    v) containerVersion=$OPTARG;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;
  esac
done

shift $((OPTIND-1))

image="${repoDir}/containers/synbold-disco-${containerVersion}.sif"

if [[ ! -f $image ]]; then
  echo "Cannot find requested container $image"
  exit 1
fi

if [[ -z "${LSB_JOBID}" ]]; then
  echo "This script must be run within a (batch or interactive) LSF job"
  exit 1
fi

sngl=$( which singularity ) ||
    ( echo "Cannot find singularity executable. Try module load singularity"; exit 1 )

if [[ ! -d "$inputDir" ]]; then
  echo "Cannot find input directory $inputDir"
  exit 1
fi

requiredInputs=("T1.nii.gz" "BOLD_d.nii.gz")

for input in ${requiredInputs[@]}; do
  if [[ ! -f "${inputDir}/${input}" ]]; then
    echo "Cannot find required input ${inputDir}/${input}"
    exit 1
  fi
done

if [[ ! -d "${outputDir}" ]]; then
  mkdir -p "$outputDir"
fi

if [[ ! -d "${outputDir}" ]]; then
  echo "Could not find or create output directory ${outputDir}"
  exit 1
fi

# Set a job-specific temp dir
jobTmpDir=$( mktemp -d -p ${SINGULARITY_TMPDIR} synb0.${LSB_JOBID}.XXXXXXXX.tmpdir ) ||
    ( echo "Could not create job temp dir ${jobTmpDir}"; exit 1 )

# Not all software uses TMPDIR
# module DEV/singularity sets SINGULARITYENV_TMPDIR=/scratch
# We will make a temp dir there and bind to /tmp in the container
export SINGULARITYENV_TMPDIR="/tmp"

# Control threads
export SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${LSB_DJOB_NUMPROC}
export SINGULARITYENV_OMP_NUMTHREADS=${LSB_DJOB_NUMPROC}

# singularity args
singularityArgs="--cleanenv --no-home \
  -B ${jobTmpDir}:/tmp \
  -B ${fsLicense}:/extra/freesurfer/license.txt \
  -B ${inputDir}:/INPUTS \
  -B ${outputDir}:/OUTPUTS"

if [[ -n "$userBindPoints" ]]; then
  singularityArgs="$singularityArgs \
  -B $userBindPoints"
fi

userArgs="$*"

echo "
--- args passed through to container ---
$*
---
"

echo "
--- Script options ---
synb0 image            : $image
Input directory        : $inputDir
Output directory       : $outputDir
Cleanup temp           : $cleanup
User bind points       : $userBindPoints
---
"

echo "
--- Container details ---"
singularity inspect $image
echo "---
"

cmd="singularity run \
  $singularityArgs \
  $image \
  $userArgs"

echo "
--- full command ---
$cmd
---
"

singExit=0
( $cmd ) || singExit=$?

if [[ $singExit -eq 0 ]]; then
  echo "Container exited with status 0"
fi

if [[ $cleanup -eq 1 ]]; then
  echo "Removing temp dir ${jobTmpDir}"
  rm -rf ${jobTmpDir}
else
  echo "Leaving temp dir ${jobTmpDir}"
fi

exit $singExit
