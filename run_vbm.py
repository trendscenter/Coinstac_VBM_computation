#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer includes the interface adapter(IA) for parsing json args to read ni pre-processing structural T1w scans (accepts BIDS format)
This layer sends the output to vbm_use_cases_layer with the appropriate inputs to run the pipeine using nipype interface

Sample run for bids input data:
python3 run_vbm.py '{"input":{"opts":{"fwhm": 7}, "BidsDir":"/computation/test_dir/bids_input_data","WriteDir":"/computation/test_dir/bids_output"}}'

Sample run for input data of nifti paths in text file:
python3 run_vbm.py '{"input":{"opts":{"fwhm": 7}, "NiftiPathsTxt":"/computation/test_dir/nifti_paths.txt","WriteDir":"/computation/test_dir/nifti_outputs"}}'
"""

import ujson as json
import warnings, os, sys, time
import nibabel as nib

## Load Nipype spm interface ##
from nipype.interfaces import spm

#Import use cases layer and entities layer code
import vbm_use_cases_layer

#Create a dict to store all paths to softwares,templates & store parameters, names of output files

template_dict={ 'spm_version':'12.7169',
                'matlab_cmd':'/opt/spm12/run_spm12.sh /opt/mcr/v92 script',
                'spm_path':'/opt/spm12',
                'tpm_path':'/opt/spm12/spm12_mcr/spm/spm12/tpm/TPM.nii',
                'transf_mat_path':'/computation/transform.mat',
                'scan_type':'T1w',
                'FWHM_SMOOTH':[6, 6, 6],
                'BIAS_REGULARISATION':0.0001,
                'FWHM_GUASSIAN_SMOOTH_BIAS':60,
                'vbm_output_dirname':'vbm_spm12',
                'vbm_qc_filename':'vbm_corr_value.txt',
                'outputs_manual_name':'outputs_description.txt',
                'outputs_manual_content':"Prefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
                "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
                "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
                'qc_readme_name':'quality_control_readme.txt',
                'qc_readme_content':"vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox"
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
"""
with warnings.catch_warnings():
    warnings.simplefilter("ignore")


def process_bids(args):
    """Runs the pre-processing pipeline on structural T1w scans in BIDS data
        Args:
            args (dictionary): {"input":{
                                        "BidsDir": {
                                                        "type": "directory",
                                                        "label": "Input Bids Directory",
                                                        },
                                        "WriteDir": {
                                                        "type": "directory",
                                                        "label": "Input Bids Directory",
                                                        }
                                        }
                                }
        Returns:
            computation_output (json): {"output": {
                                                  "success": {
                                                    "type": "boolean"
                                                  },
                                                   "vbmdirs": {
                                                    "type": "array",
                                                    "contains": ["string"]
                                                  },
                                                  "wc1files": {
                                                    "type": "array",
                                                    "contains": ["string"]
                                                  }
                                                  }
                                        }
        Comments:
            After verifying the BIDS format , the bids_dir and write_dir along with pre-processing specific pipeline options
            are sent to vbm_use_cases_layer for running the pipeline

        """
    BidsDir = args['input']['BidsDir']
    WriteDir = args['input']['WriteDir']

    if 'opts' in args['input']:opts=args['input']['opts']
    else:opts=None

    # Check if input_bids_dir is in BIDS format using bids-validator tool
    cmd = "bids-validator {0}".format(BidsDir)
    bids_process = os.popen(cmd).read()

    # Check if it has T1w data and write permissions to Write dir
    if bids_process and template_dict['scan_type'] in bids_process and os.access(WriteDir, os.W_OK):
        return vbm_use_cases_layer.run_bids_pipeline(bids_dir=BidsDir, write_dir=WriteDir,pipeline_opts=opts,**template_dict)

    else:
        sys.stderr.write('Incompatible Bids directory format or write directory does not have permissions')



def process_niftis(args):
    """Runs the pre-processing pipeline on structural T1w nifti scans from paths in the text file
            Args:
                args (dictionary): {"input":{
                                            "NiftiFile": {
                                                            "type": "directory",
                                                            "label": "text file with complete T1w nifti paths",
                                                            },
                                            "WriteDir": {
                                                            "type": "directory",
                                                            "label": "Input Bids Directory",
                                                            }
                                            }
                                    }
            Returns:
                computation_output (json): {"output": {
                                                      "success": {
                                                        "type": "boolean"
                                                      },
                                                       "vbmdirs": {
                                                        "type": "array",
                                                        "contains": ["string"]
                                                      },
                                                      "wc1files": {
                                                        "type": "array",
                                                        "contains": ["string"]
                                                      }
                                                      }
                                            }
            Comments:
                After verifying the nifti paths , the paths to nifti files and write_dir along with pre-processing specific pipeline options
                are sent to vbm_use_cases_layer for running the pipeline

            """

    paths_file = args['input']['NiftiPathsTxt']
    WriteDir = args['input']['WriteDir']

    if 'opts' in args['input']:opts=args['input']['opts']
    else:opts=None

    # Reach each line in nifti_paths.txt into niftis variable
    niftis=[]
    with open(paths_file, "r") as f:
        for line in f:
            niftis.append(line.rstrip(('\n')))

    if os.access(WriteDir, os.W_OK):
        return vbm_use_cases_layer.run_nifti_pipeline(nii_files=niftis, write_dir=WriteDir, pipeline_opts=opts, **template_dict)
    else:
        sys.stderr.write('write directory does not have permissions')

def software_check():
    """This function returns the spm standalone version installed inside the docker
    """
    spm.SPMCommand.set_mlab_paths(matlab_cmd=template_dict['matlab_cmd'], use_mcr=True)
    return(spm.SPMCommand().version)

if __name__=='__main__':

    # Check if spm is running and give 15 secs to connect spm
    spm_check = software_check()
    time.sleep(15)
    if spm_check is None: 
        raise ValueError("spm unable to start in vbm docker")

    # The following block of code assigns the appropriate pre-processing function for input data format, based on Bids or nifti file paths in text file #
    #args = json.loads(sys.stdin.read()) This is not working for reading in json strings

    args = json.loads(sys.argv[1])

    if args and spm_check==template_dict['spm_version']:
        if 'BidsDir' in args['input'] and 'WriteDir' in args['input']:
            computation_output=process_bids(args)
            sys.stdout.write(computation_output)
        elif 'NiftiPathsTxt' in args['input'] and 'WriteDir' in args['input']:
            computation_output=process_niftis(args)
            sys.stdout.write(computation_output)
        else:sys.stderr.write('Cannot read input data from json')

    else:
        raise ValueError("Incompatible json or spm unable to start in vbm docker")