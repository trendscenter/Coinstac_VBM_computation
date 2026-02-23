#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer includes the interface adapter(IA) for parsing json args to read structural T1w scans (formats:BIDS, nifti files, dicoms)
This layer sends the output to vbm_use_cases_layer with the appropriate inputs to run the pipeine using nipype interface

Sample run examples:
python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/3D_T1"]]}}
3D_T1 contains T1w dicoms

python3 run_fmri.py python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/sub1_t1w.nii","/input/local0/simulatorRun/sub1_t1w.nii.gz"]]}}

python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/BIDS_DIR"]]}}
BIDS_DIR contains bids data

success=True means program finished execution , despite the success or failure of the code
This is to indicate to coinstac that program finished execution
"""
import contextlib
import os
import sys
import traceback
import warnings


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

    oldstdchannel = None
    dest_file = None
    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


import ujson as json
import numpy as np, scipy.io, spm_matrix as s

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
# Load Nipype spm interface #
from nipype.interfaces import spm
import vbm_use_cases_layer,vbm_standalone_use_cases_layer

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

# Store all paths/options/output names in a single config dictionary.
from vbm_config import template_dict
from vbm_args import args_parser
"""
More info. on keys in template_dict

spm_path is path to spm software inside docker .
SPM is Statistical Parametric Mapping toolbox for matlab
Info. from http://www.fil.ion.ucl.ac.uk/spm/
"Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to test hypotheses about functional imaging data.
These ideas have been instantiated in software that is called SPM.
The SPM software package has been designed for the analysis of brain imaging data sequences.
The sequences can be a series of images from different cohorts, or time-series from the same subject.
The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG."

tpm_path is the path where the SPM structural template nifti file is stored
This file is used to :
1) Perform segmentation in the VBM pipeline
2) Compute correlation value to smoothed, warped grey matter from output of pipeline, which is stored in the vbm_qc_filename

transf_mat_path is the path to the transformation matrix used in running the reorient step of the pipeline
scan_type is the type of structural scans on which is accepted by this pipeline
FWHM_SMOOTH is the full width half maximum smoothing kernel value in mm in x,y,z directions
vbm_output_dirname is the name of the output directory to which the outputs from this pipeline are written to
vbm_qc_filename is the name of the VBM quality control text file , which is placed in vbm_output_dirname

For nifti files , it is assumed that they are T1w (T1 weighted) type of scans
FWHM_SMOOTH is an optional parameter that can be passed as json in args['input']['opts']

json output description
                    "message"-This string is used by coinstac to display output message to the user on the UI after computation is finished
                    "download_outputs"-Zipped directory where outputs are stored
                    "display"-base64 encoded string of segmented gray matter normalized output nifti
"""
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")



def software_check():
    """This function returns the spm standalone version installed inside the docker
    """
    spm.SPMCommand.set_mlab_paths(
        matlab_cmd=template_dict['matlab_cmd'], use_mcr=True)
    return (spm.SPMCommand().version)

def convert_reorientparams_save_to_mat_script():
    try:
        pi = 22 / 7
        scipy.io.savemat('/computation/transform.mat',
                         mdict={'M': np.around(s.spm_matrix([template_dict['options_reorient_params_x_mm'],
                                                             template_dict['options_reorient_params_y_mm'],
                                                             template_dict['options_reorient_params_z_mm'],
                                                             template_dict['options_reorient_params_pitch'] * (
                                                                         pi / 180),
                                                             template_dict['options_reorient_params_roll'] * (pi / 180),
                                                             template_dict['options_reorient_params_yaw'] * (pi / 180),
                                                             template_dict['options_reorient_params_x_scaling'],
                                                             template_dict['options_reorient_params_y_scaling'],
                                                             template_dict['options_reorient_params_z_scaling'],
                                                             template_dict['options_reorient_params_x_affine'],
                                                             template_dict['options_reorient_params_y_affine'],
                                                             template_dict['options_reorient_params_z_affine']], 1),
                                               decimals=4)[0]})
    except Exception as e:
        sys.stderr.write('Unable to convert reorientation params to transform.mat Error_log:'+str(e)+str(traceback.format_exc()))


def data_parser(args):
    """ This function parses the type of data i.e BIDS, nifti files or Dicoms
    and passes them to vbm_use_cases_layer.py
    """
    # if template_dict['standalone']:
    #     data = [args['state']['baseDirectory'] + '/' + file_names for file_names in args['input']['site_data']]
    # else:
    #     data = [args['state']['baseDirectory'] + '/' + subject[0] for subject in args['input']['covariates'][0][0][1:]]
    # We are only concerned with pre-processing
    data = [args['state']['baseDirectory'] + '/' + file_names for file_names in args['input']['covariates']]

    WriteDir = args['state']['outputDirectory']

    # Check if data has nifti files
    for x in data:
        if not os.path.isfile(x):
            raise Exception("File does not exist or can't be read: " + str(x));
    if not os.access(WriteDir, os.W_OK):
        raise Exception("Output write permissions denied: " + str(WriteDir))
    nifti_paths = data
    # if template_dict['standalone']:
    #     computation_output = vbm_standalone_use_cases_layer.setup_pipeline(
    #         data=nifti_paths,
    #         write_dir=WriteDir,
    #         data_type='nifti',
    #         **template_dict)
    #     sys.stdout.write(computation_output)
    # else:
    #     computation_output = vbm_use_cases_layer.setup_pipeline(
    #         data=nifti_paths,
    #         write_dir=WriteDir,
    #         data_type='nifti',
    #         **template_dict)
    #     sys.stdout.write(computation_output)
    # We are only concerned with pre-processing
    covariates = args['input']['covariates']

    computation_output = vbm_standalone_use_cases_layer.setup_pipeline(
        data=nifti_paths,
        write_dir=WriteDir,
        covars=covariates,
        data_type='nifti',
        **template_dict)
    return computation_output

def start(args):
    # Check if spm is running
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm_check = software_check()
    if spm_check != template_dict['spm_version']:
        raise EnvironmentError("spm unable to start in vbm docker")

    # Read json args
    # args = json.loads(sys.stdin.read())

    # Parse args
    args_output = args_parser(args, template_dict)
    if args_output is not None:
        return args_output

    #Convert reorient params to mat file if they exist
    convert_reorientparams_save_to_mat_script()

    # Parse input data and run the code
    return data_parser(args)
