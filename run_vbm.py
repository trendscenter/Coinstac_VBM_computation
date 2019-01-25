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


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

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
import warnings, os, glob, sys
import nibabel as nib

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
# Load Nipype spm interface #
from nipype.interfaces import spm
import vbm_use_cases_layer

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

#Create a dictionary to store all paths to softwares,templates & store parameters, names of output files

template_dict = {
    'spm_version':
    '12.7169',
    'matlab_cmd':
    '/opt/spm12/run_spm12.sh /opt/mcr/v92 script',
    'spm_path':
    '/opt/spm12/fsroot',
    'tpm_path':
    '/opt/spm12/fsroot/spm/spm12/tpm/TPM.nii',
    'transf_mat_path':
    os.path.join('/computation', 'transform.mat'),
    'scan_type':
    'T1w',
    'FWHM_SMOOTH': [10, 10, 10],
    'bounding_box':
    '',
    'reorient_params':
    '',
    'BIAS_REGULARISATION':
    0.0001,
    'FWHM_GAUSSIAN_SMOOTH_BIAS':
    60,
    'vbm_output_dirname':
    'vbm_spm12',
    'output_zip_dir':
    'vbm_outputs',
    'qa_flagged_filename':
    'QA_flagged_subjects.txt',
    'display_image_name':
    'wc1Re.png',
    'display_pngimage_name':
    'Gray matter (Normalized)',
    'cut_coords': (0, 0, 0),
    'display_nifti':
    'wc1Re.nii',
    'qc_nifti':
    'swc1*nii',
    'vbm_qc_filename':
    'vbm_corr_value.txt',
    'outputs_manual_name':
    'outputs_description.txt',
    'coinstac_display_info':
    'Please read outputs_description.txt for description of pre-processed output files and quality_control_readme.txt for quality control measurement.'
    'These files are placed under the pre-processed data.',
    'flag_warning':
    ' QC warning: Only half of the data is pre-processed or passed the QA, please check the data!',
    'bids_outputs_manual_content':
    "Prefixes descriptions for segmented images:c1-Gray matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'nifti_outputs_manual_content':
    "sub-1,sub-2,sub-* denotes each nifti file with respect to the order of the nifti paths given"
    "\nPrefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'dicoms_outputs_manual_content':
    "sub-1,sub-2,sub-* denotes each session's dicom directory with respect to the order of the dicom directory paths given"
    "\nPrefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
    "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
    "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'qc_readme_name':
    'quality_control_readme.txt',
    'qc_readme_content':
    "In each subject's anat/vbm_spm12 directory,vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12 adult brain tissue probability maps (TPM.nii) file from SPM12 toolbox"
    "\nIf your subjects are kids/adolescent, take this into considersation as correlation value may not pass threshold. For an adult study, its suggested that scans with correlation value <=0.91 "
    "\nshould be manually looked into for possible reorientation to the ac-pc line. Subjects that do not pass this QA metric will be saved in vbm_outputs/QA_flagged_subjects.txt"
}
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


def args_parser(args):
    """ This function extracts options from arguments
    """
    if 'options' in args['input']:
        template_dict['FWHM_SMOOTH'] = [float(args['input']['options'])]*3

    if 'registration_template' in args['input']:
        if os.path.isfile(args['input']['registration_template']) and (str(
            ((nib.load(template_dict['tpm_path'])).shape)) == str(
                ((nib.load(args['input']['registration_template'])).shape))):
            template_dict['tpm_path'] = args['input']['registration_template']
        else:
            sys.stdout.write(
                json.dumps({
                    "output": {
                        "message": "Non-standard Registration template "
                    },
                    "cache": {},
                    "success": True
                }))
            sys.exit()


def data_parser(args):
    """ This function parses the type of data i.e BIDS, nifti files or Dicoms
    and passes them to vbm_use_cases_layer.py
    """
    data = args['input']['data']
    WriteDir = args['state']['outputDirectory']

    # Check if data is BIDS
    if os.path.isfile(os.path.join(data[0],
                                   'dataset_description.json')) and os.access(
                                       WriteDir, os.W_OK):
        cmd = "bids-validator {0}".format(data[0])
        bids_process = os.popen(cmd).read()
        bids_dir = data[0]
        if bids_process and template_dict['scan_type'] in bids_process:
            computation_output = vbm_use_cases_layer.setup_pipeline(
                data=bids_dir,
                write_dir=WriteDir,
                data_type='bids',
                **template_dict)
    # Check if data has nifti files
    elif [x
          for x in data if os.path.isfile(x)] and os.access(WriteDir, os.W_OK):
        nifti_paths = data
        computation_output = vbm_use_cases_layer.setup_pipeline(
            data=nifti_paths,
            write_dir=WriteDir,
            data_type='nifti',
            **template_dict)
        sys.stdout.write(computation_output)
    # Check if inputs are dicoms
    elif [x
          for x in data if os.path.isdir(x)] and os.access(WriteDir, os.W_OK):
        dicom_dirs = list()
        for dcm in data:
            if os.path.isdir(dcm) and os.listdir(dcm):
                dicom_file = glob.glob(dcm + '/*')[0]
                dicom_header_info = os.popen('strings' + ' ' + dicom_file +
                                             '|grep DICM').read()
                if 'DICM' in dicom_header_info: dicom_dirs.append(dcm)
        computation_output = vbm_use_cases_layer.setup_pipeline(
            data=dicom_dirs,
            write_dir=WriteDir,
            data_type='dicoms',
            **template_dict)
        sys.stdout.write(computation_output)
    else:
        sys.stdout.write(
            json.dumps({
                "output": {
                    "message":
                    "Input data not found/Can not write to target directory"
                },
                "cache": {},
                "success": True
            }))


if __name__ == '__main__':
    # Check if spm is running
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm_check = software_check()
    if spm_check != template_dict['spm_version']:
        raise EnvironmentError("spm unable to start in vbm docker")

    #Read json args
    args = json.loads(sys.stdin.read())

    #Parse args
    args_parser(args)

    #Parse input data
    data_parser(args)
