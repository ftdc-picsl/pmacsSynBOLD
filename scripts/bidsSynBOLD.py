#!/usr/bin/env python

import argparse
import bids
import json
import os
import shutil
import subprocess
import sys
import tempfile

def _filter_pybids_none_any(dct):
    import bids
    return {
        k: bids.layout.Query.NONE
        if v is None
        else (bids.layout.Query.ANY if v == "*" else v)
        for k, v in dct.items()
    }

def get_bids_filters(value):
    from bids.layout import Query
    from json import loads
    from pathlib import Path

    if value and Path(value).exists():
        try:
            filters = loads(Path(value).read_text(), object_hook=_filter_pybids_none_any)
        except Exception as e:
            raise Exception("Unable to parse BIDS filter file. Check that it is "
                            "valid JSON.")
    else:
        raise Exception("Unable to load BIDS filter file " + value)

    # unserialize pybids Query enum values
    for acq, _filters in filters.items():
        filters[acq] = {
            k: getattr(Query, v[7:-4])
            if not isinstance(v, Query) and "Query" in v
            else v
            for k, v in _filters.items()
        }
    return filters

def get_bold_images(layout, subject_label, session_label, filter_criteria):
    bold_images = layout.get(subject=subject_label, session=session_label, suffix='bold', extension=['nii.gz'],
                             **filter_criteria)

    # Images grouped by task and run, this is similar to HCP style
    # Run synbold once per group, make an fmap with intendedfor all bolds in that group
    grouped_images = {}

    for file in bold_images:

        group_label = f"task-{file.get_entities()['task']}"

        # check for optional run label
        if 'run' in file.get_entities():
            group_label = f"{group_label}_run-{file.get_entities()['run']}"

        if group_label not in grouped_images:
            grouped_images[group_label] = []
        grouped_images[group_label].append(file)

    return grouped_images


def get_t1w_skull_stripped(dataset, participant_label, session_label, t1w_filename):
    # Look in the mask dataset for a T1w mask matching the T1w image
    mask_dir = os.path.join(dataset, f"sub-{participant_label}", f"ses-{session_label}", 'anat')
    # Get all json files in the mask directory
    mask_sidecars = [f for f in os.listdir(mask_dir) if f.endswith('_mask.json')]

    t1w_mask = None
    t1w_skull_stripped_path = None

    found_mask = False

    for mask_sidecar in mask_sidecars:
        # Load the sidecar
        with open(os.path.join(mask_dir, mask_sidecar)) as json_file:
            mask_json = json.load(json_file)
            # Check if the T1w image matches the T1w image in the mask sidecar
            for source in mask_json['Sources']:
                if source.endswith(t1w_filename):
                    # Found a match
                    t1w_mask = os.path.join(mask_dir, mask_sidecar.replace('.json', '.nii.gz'))
                    found_mask = True
                    break
            if found_mask:
                break

    if t1w_mask is not None:
        print("Using T1w mask " + t1w_mask)
        t1w_skull_stripped_path = os.path.join(working_dir, "t1w_skull_stripped.nii.gz")
        t1w_is_skull_stripped = True
        subprocess.run(["fslmaths", t1w_path, "-mas", t1w_mask, t1w_skull_stripped_path])

    return t1w_skull_stripped_path

# Helps with CLI help formatting
class RawDefaultsHelpFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass

# parse arguments with argparse
# if no args, print usage
parser = argparse.ArgumentParser(formatter_class=RawDefaultsHelpFormatter,
                                 prog="bidsSynBOLD", add_help = False, description='''
Runs synBOLD on BIDS data.

A BIDs filter can be used to select a subset of the data.

Limitations:

 * Only one T1w image per subject/session is supported. If more than one is found, the
script will exit. Use a filter or the --t1w-image-suffix option to select a specific T1w image.

 * Data must be organized as sub-<participant_label>/ses-<session_label>.

 * Axial scans are assumed, phase encoding axis must be AP (j-) or PA (j).

 * By default, images are grouped by task. So task-rest and task-somethingelse will have their own field map.
   But there is no provision to group by other entities like run, acq, etc. This can be done with a BIDS filter.


Requires: FSL, singularity
                                 ''')

required = parser.add_argument_group('Required arguments')
required.add_argument("-c", "--container", help="synBOLD container to use", type=str, required=True)
required.add_argument("--bids-dataset", help="Input BIDS dataset", type=str, required=True)
required.add_argument("--fs-license-file", help="Freesurfer license file", type=str, required=True)
required.add_argument("--participant-label", help="participant to process", type=str, required=True)
required.add_argument("--session-label", help="session to process", type=str, required=True)

optional = parser.add_argument_group('Optional arguments')
optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
optional.add_argument("-f", "--bids-filter", help="BIDS filter file", type=str, default=None)
optional.add_argument("-n", "--num-threads", help="Number of computational threads", type=int, default=1)
optional.add_argument("-t", "--total-readout-time", help="Total readout time for the BOLD. Some older DICOM files "
                      "do not provide this information, so it can be specified manually. Ignored if the BIDS "
                      "sidecar contains total readout time", type=str, default=None)
optional.add_argument("--combine-all-bold", help="Combine all bold images into one group. By default, images will be "
                      "grouped by task and run, so each run of each task gets its own field map", action='store_true')
optional.add_argument("--bold-input-volumes", help="Number of volumes to use from the bold input. The first N volumes "
                      "will be taken from the first input in each group, and averaged. Set to 0 to use an entire "
                      "series.", type=int, default=1)
optional.add_argument("--bold-input-4d", help="Use the entire 4D time series as input to synBOLD. Use this to have "
                      "synBOLD motion-correct the data before defining an average volume as a reference image",
                      action='store_true')
optional.add_argument("--no-smooth", help="Do not smooth the distorted reference image before running topup",
                      action='store_true')
optional.add_argument("--direct-field-mapping", help="Use direct field mapping instead of producing EPI image pairs",
                      action='store_true')
optional.add_argument("--t1w-image-suffix", help="Use a specific T1w head image suffix. Eg, 'acq-mprage_T1w.nii.gz' selects "
                      "sub-participant/ses-session/sub-participant_ses-session_acq-mprage_T1w.nii.gz'. "
                      "Using this overrides BIDS filters for the T1w", type=str, default=None)
optional.add_argument("--t1w-mask-dataset", help="BIDS dataset to use for brain masking the T1w. If not specified, "
                      "the current dataset will be searched for a brain mask. The last resort is to pass the unstripped "
                      "T1w to synBOLD. It is highly recommended to pass a high quality brain mask.", type=str, default=None)
optional.add_argument("-w", "--work-dir", help="Temp dir for intermediate output, defaults to system "
                      "TMPDIR if defined, otherwise '/scratch'", type=str, default=os.environ.get('TMPDIR', '/scratch'))

args = parser.parse_args()

if shutil.which('singularity') is None:
    raise RuntimeError('singularity executable not found')

# open the BIDS dataset
layout = bids.BIDSLayout(args.bids_dataset, validate=False)

# Get filter if provided
bids_filters = {}
if args.bids_filter is not None:
    bids_filters = get_bids_filters(args.bids_filter)

# Get all dMRI data files for the given subject and session
bold_groups = get_bold_images(layout, args.participant_label, args.session_label, bids_filters)

if (args.combine_all_bold):
    # Combine all bold into one group
    bold_groups = {'combined': [item for sublist in bold_groups.values() for item in sublist]}

if args.t1w_image_suffix is not None:
    t1w_files = [os.path.join(args.bids_dataset, f"sub-{args.participant_label}", f"ses-{args.session_label}", 'anat',
                f"sub-{args.participant_label}_ses-{args.session_label}_{args.t1w_image_suffix}")]
else:
    # return type files gets actual files not BIDSFile objects
    t1w_files = layout.get(subject=args.participant_label, session=args.session_label, suffix='T1w',
                           return_type='file', extension=['nii.gz'], **bids_filters)

if (len(t1w_files) > 1):
    print("More than one T1w image found for subject " + args.participant_label + " session " + args.session_label +
          ". Need a more specific filter")
    sys.exit(1)

# Check t1w exists
if len(t1w_files) == 0:
    print("No T1w image found for subject " + args.participant_label + " session " + args.session_label)
    sys.exit(1)
if not os.path.exists(t1w_files[0]):
    print("T1w image " + t1w_files[0] + " not found")
    sys.exit(1)

t1w_path = t1w_files[0]

# Get filename from the full path
t1w_filename = os.path.basename(t1w_path)

# If this is None, use synBOLD's internal brain masking
t1w_skull_stripped_path = None
t1w_is_skull_stripped = False

# will be cleaned up after the script finishes
working_dir_tmpdir = tempfile.TemporaryDirectory(prefix=f"bids-synBOLD.", dir=args.work_dir, ignore_cleanup_errors=True)
working_dir = working_dir_tmpdir.name

# Print the files we're processing
print("Processing subject " + args.participant_label + " session " + args.session_label)
print("T1w: " + t1w_path)
print("BOLD groups: " + str(bold_groups))

if args.t1w_mask_dataset is not None:
    # Look in the mask dataset for a T1w mask matching the T1w image
    t1w_skull_stripped_path = get_t1w_skull_stripped(args.t1w_mask_dataset, args.participant_label, args.session_label,
                                                     t1w_filename)
else:
    # Search current dataset for a brain mask
    t1w_skull_stripped_path = get_t1w_skull_stripped(args.bids_dataset, args.participant_label, args.session_label,
                                                     t1w_filename)

if t1w_skull_stripped_path is not None:
    t1w_is_skull_stripped = True
else:
    print("No T1w brain mask found for subject " + args.participant_label + " session " + args.session_label)

# Run synBOLD on each group of BOLD
for group in bold_groups:
    print("Running synBOLD on group " + group)
    # Check all images in the group have the same phase encoding direction
    pe_direction = None
    pe_direction_consistent = True
    group_bold_images = bold_groups[group]
    for bold_image in group_bold_images:
        print("Checking phase encoding direction for " + bold_image.filename)
        # Get the phase encoding direction from the image metadata
        # Metadata is stored in the JSON sidecar file, same file name but with .json instead of
        # .nii.gz
        sidecar_file = bold_image.path.replace('.nii.gz', '.json')
        with open(sidecar_file) as sidecar_fh:
            sidecar = json.load(sidecar_fh)
            if pe_direction is not None:
                if pe_direction != sidecar['PhaseEncodingDirection']:
                    print("Phase encoding direction for " +  bold_image.filename + " does not match previous images in group")
                    pe_direction_consistent = False
            else:
                pe_direction = sidecar['PhaseEncodingDirection']
                if pe_direction not in ['j', 'j-']:
                    print(f"Phase encoding direction {pe_direction} in {group} not supported by synBOLD. Skipping group")
                    continue

    if not pe_direction_consistent:
        print("Phase encoding direction not consistent for group " + group + ". Skipping group")
        continue

    print("Phase encoding direction for group is " + pe_direction)

    # Run synBOLD on the group using the first series

    bold_ref = group_bold_images[0]

    # Get total readout time from the first image or use command line alternative (needed for older DICOM files)
    bold_ref_total_readout_time = args.total_readout_time
    bold_ref_effective_echo_spacing = 0.0
    bold_ref_sidecar_file = bold_ref.path.replace('.nii.gz', '.json')

    with open(bold_ref_sidecar_file) as sidecar_fh:
        sidecar = json.load(sidecar_fh)
        try:
            bold_ref_total_readout_time = sidecar['TotalReadoutTime']
        except KeyError:
            print("No total readout time in sidecar for " + bold_ref.path)
            if bold_ref_total_readout_time is None:
                print("No total readout time in sidecar and no readout time specified on command line")
                sys.exit(1)
        try:
            phase_encode_dir = sidecar['PhaseEncodingDirection']
            if phase_encode_dir not in ['j', 'j-']:
                print("Phase encoding direction " + phase_encode_dir + " not supported by synBOLD")
                sys.exit(1)
        except KeyError:
            print("No phase encoding direction in sidecar for " + bold_ref.path)
            sys.exit(1)
        try:
            bold_ref_effective_echo_spacing = sidecar['EffectiveEchoSpacing']
        except KeyError:
            print("WARNING: No effective echo spacing in sidecar for " + bold_ref.path)

    synBOLD_env = os.environ.copy()
    synBOLD_env['SINGULARITYENV_TMPDIR'] = '/tmp'
    synBOLD_env['SINGULARITYENV_OMP_NUM_THREADS'] = str(args.num_threads)
    synBOLD_env['SINGULARITYENV_ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(args.num_threads)

    # These are inputs and output for this group under the top working directory
    tmp_input_dir = os.path.join(working_dir, f"{group}_synBOLD_input")
    tmp_output_dir = os.path.join(working_dir, f"{group}_synBOLD_output")
    # Mount this to /tmp for the container
    tmp_singularity_dir = os.path.join(working_dir, f"{group}_synBOLD_tmpdir")

    os.makedirs(tmp_input_dir)
    os.makedirs(tmp_output_dir)
    os.makedirs(tmp_singularity_dir)

    # Get the first image from the bold time series - this is the reference image for synbold
    # input to synBOLD can be 3D or 4D. If 4D, it will do moco unless told otherwise
    bold_input = os.path.join(tmp_input_dir, 'BOLD_d.nii.gz')

    if args.bold_input_4d:
        shutil.copy(bold_ref.path, bold_input)
    else:
        num_ref_volumes = args.bold_input_volumes

        if num_ref_volumes == 0:
            # fslroi uses -1 to mean all volumes
            num_ref_volumes = -1

        subprocess.run(['fslroi', bold_ref.path, bold_input, '0', str(num_ref_volumes)], env=synBOLD_env)

    # fslroi exits 0 if the file doesn't exist, so check manually
    if not os.path.exists(bold_input):
        print("Could not extract reference BOLD image from " + bold_ref.path)
        sys.exit(1)

    t1_input = os.path.join(tmp_input_dir, 'T1.nii.gz')

    if t1w_is_skull_stripped:
        shutil.copy(t1w_skull_stripped_path, t1_input)
    else:
        shutil.copy(t1w_path, t1_input)

    # Get synBOLD output and copy to fmap/
    synBOLD_cmd_list = ['singularity', 'run', '--cleanenv', '--no-home', '-B', f"{os.path.realpath(tmp_input_dir)}:/INPUTS",
                        '-B', f"{os.path.realpath(tmp_output_dir)}:/OUTPUTS",
                        '-B', f"{os.path.realpath(tmp_singularity_dir)}:/tmp",
                        '-B', f"{os.path.realpath(args.fs_license_file)}:/opt/freesurfer/license.txt",
                        args.container, '--total_readout_time', str(bold_ref_total_readout_time)]

    if t1w_is_skull_stripped:
        synBOLD_cmd_list.append('--skull_stripped')

    if args.no_smooth:
        synBOLD_cmd_list.append('--no_smoothing')

    if not args.direct_field_mapping:
        synBOLD_cmd_list.append('--no_topup')

    print("---synBOLD call---\n" + " ".join(synBOLD_cmd_list) + "\n---")

    subprocess.run(synBOLD_cmd_list, env=synBOLD_env)

    # output_dir is fmap/ under the session directory in the dataset
    output_dir = os.path.join(args.bids_dataset, f"sub-{args.participant_label}", f"ses-{args.session_label}", 'fmap')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.direct_field_mapping:

        fmap_file_name = None
        magnitude_file_name = None

        if group == 'combined':
            fmap_file_name = f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_fieldmap.nii.gz"
        else:
            fmap_file_name = f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_{group}_fieldmap.nii.gz"

        magnitude_file_name = fmap_file_name.replace('fieldmap', 'magnitude')

        shutil.copy(os.path.join(tmp_output_dir, 'BOLD_u_3D.nii.gz'), os.path.join(output_dir, magnitude_file_name))
        shutil.copy(os.path.join(tmp_output_dir, 'topup_results_field.nii.gz'), os.path.join(output_dir, fmap_file_name))

        magnitude_sidecar_file = magnitude_file_name.replace('.nii.gz', '.json')

        fmap_sidecar_file = fmap_file_name.replace('.nii.gz', '.json')

        # Set magnitude sidecar

        magnitude_sidecar = {
            "EffectiveEchoSpacing": 0.0,
            "Modality": "MR",
            "SeriesDescription": "synBOLD magnitude image",
            "TotalReadoutTime": 0.0000001
        }

        with open(os.path.join(output_dir, magnitude_sidecar_file), 'w') as magnitude_sidecar_fh:
            json.dump(magnitude_sidecar, magnitude_sidecar_fh, indent=2, sort_keys=True)

        # Set fmap sidecar
        fmap_sidecar = {
            "Modality": "MR",
            "SeriesDescription": "synBOLD fieldmap",
            "Units": "Hz"
        }

        fmap_sidecar['IntendedFor'] = [os.path.join(f"ses-{args.session_label}", 'func', file.filename)
                                    for file in group_bold_images]

        with open(os.path.join(output_dir, fmap_sidecar_file), 'w') as fmap_sidecar_fh:
            json.dump(fmap_sidecar, fmap_sidecar_fh, indent=2, sort_keys=True)

    else:
        # Copy the distorted and undistorted 3D images to the fmap directory
        if pe_direction == 'j':
            bold_pe_label = 'PA'
            rpe_label = 'AP'
            rpe_direction = 'j-'
        if pe_direction == 'j-':
            bold_pe_label = 'AP'
            rpe_label = 'PA'
            rpe_direction = 'j'
        if group == 'combined':
            epi_bold_file_name = \
                f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_dir-{bold_pe_label}_epi.nii.gz"
            rpe_bold_file_name = \
                f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_dir-{rpe_label}_epi.nii.gz"
        else:
            epi_bold_file_name = \
                f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_{group}_dir-{bold_pe_label}_epi.nii.gz"
            rpe_bold_file_name = \
                f"sub-{args.participant_label}_ses-{args.session_label}_acq-synBOLD_{group}_dir-{rpe_label}_epi.nii.gz"

        shutil.copy(os.path.join(tmp_output_dir, 'BOLD_d_3D_smoothed.nii.gz'), os.path.join(output_dir, epi_bold_file_name))
        shutil.copy(os.path.join(tmp_output_dir, 'BOLD_s_3D.nii.gz'), os.path.join(output_dir, rpe_bold_file_name))

        epi_bold_sidecar_file = epi_bold_file_name.replace('.nii.gz', '.json')

        epi_bold_sidecar = {
            "EffectiveEchoSpacing": bold_ref_effective_echo_spacing,
            "Modality": "MR",
            "PhaseEncodingDirection": pe_direction,
            "SeriesDescription": "synBOLD distorted image",
            "TotalReadoutTime": bold_ref_total_readout_time
        }

        with open(os.path.join(output_dir, epi_bold_sidecar_file), 'w') as epi_sidecar_fh:
            json.dump(epi_bold_sidecar, epi_sidecar_fh, indent=2, sort_keys=True)

        rpe_bold_sidecar_file = rpe_bold_file_name.replace('.nii.gz', '.json')

        rpe_bold_sidecar = {
            "EffectiveEchoSpacing": 0.0,
            "Modality": "MR",
            "PhaseEncodingDirection": rpe_direction,
            "SeriesDescription": "synBOLD undistorted image",
            "TotalReadoutTime": 0.0000001
        }

        with open(os.path.join(output_dir, rpe_bold_sidecar_file), 'w') as rpe_sidecar_fh:
            json.dump(rpe_bold_sidecar, rpe_sidecar_fh, indent=2, sort_keys=True)



