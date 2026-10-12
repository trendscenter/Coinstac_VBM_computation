#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""VBM standalone preprocessing pipeline orchestration."""

import base64
import glob
import math
import os
import shutil
import sys
import warnings

import nibabel as nib
import nipype.pipeline.engine as pe
import numpy as np
from nilearn import plotting

from . import vbm_entities_layer, vbm_spm12_file_output
from .vbm_io_utils import remove_tmp_files, write_readme_files
from .pipeline_common import create_workflow_input, stdchannel_redirected
from .pipeline_nodes import create_pipeline_nodes, smooth_images

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

# Stop printing nipype.workflow info to stdout.
from nipype import logging

logging.getLogger('nipype.workflow').setLevel('CRITICAL')


def setup_pipeline(data='', write_dir='', covars='', data_type=None, **template_dict):
    """setup the pre-processing pipeline on T1W scans
        Args:
            data (array) : Input data
            write_dir (string): Directory to write outputs
            data_type (string): BIDS, niftis, dicoms
            covariates (object): Input Covariates if they exist
            template_dict ( dictionary) : Dictionary that stores all the paths, file names, software locations
        Returns:
            computation_output (json): {"output": {
                                                  "success": {
                                                    "type": "boolean"
                                                  },
                                                   "message": {
                                                    "type": "string",
                                                  },
                                                   "download_outputs": {
                                                    "type": "string",
                                                  },
                                                   "display": {
                                                    "type": "string",
                                                  }
                                                  }
                                        }
        Comments:
            After setting up the pipeline here , the pipeline is run with run_pipeline function

        """
    [reorient, datasink, vbm_preprocess] = create_pipeline_nodes(**template_dict)

    if data_type == 'nifti':
        smri_data = data
        return run_pipeline(
            write_dir,
            smri_data,
            reorient,
            datasink,
            vbm_preprocess,
            covars,
            data_type='nifti',
            **template_dict)
    elif data_type == 'dicoms':
        smri_data = data
        return run_pipeline(
            write_dir,
            smri_data,
            reorient,
            datasink,
            vbm_preprocess,
            covars,
            data_type='dicoms',
            **template_dict)


def nii_to_image_converter(write_dir, label, **template_dict):
    """This function converts nifti to png image for displaying on coinstac web gui
    in this case : wc1*.nii
    """
    file = os.path.join(write_dir, template_dict['display_nifti'])

    mask = nib.load(file)
    new_data = mask.get_fdata()
    clipped_img = nib.Nifti1Image(new_data, mask.affine, mask.header)

    plotting.plot_anat(
        clipped_img,
        cut_coords=(0, 0, 0),
        annotate=False,
        draw_cross=False,
        output_file=os.path.join(write_dir,
                                 template_dict['display_image_name']),
        display_mode='ortho',
        title=label + ' ' + template_dict['display_pngimage_name'],
        colorbar=False)


def get_corr(segmented_file, write_dir, sub_id, **template_dict):
    """This function computes correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox """

    def extract_data(file):
        a_data = nib.load(file)
        t_data = a_data.get_fdata()
        if len(t_data.shape) == 4:
            tx, ty, tz, im = t_data.shape
            st_data = t_data[:, :, :, 0]
        else:
            tx, ty, tz = t_data.shape
            st_data = t_data
        stn_data = st_data.reshape((tx * ty * tz, 1), order='F')
        stn_data = np.nan_to_num(stn_data)
        return stn_data

    cstn_data = extract_data(template_dict['tpm_path'])
    cre_data = extract_data(segmented_file)
    indices = np.logical_and(cstn_data != 0, cre_data != 0)
    fcstn_data, fcre_data = cstn_data[indices], cre_data[indices]

    a = fcstn_data - np.mean(fcstn_data)
    b = fcre_data - np.mean(fcre_data)
    covalue = (a * b).sum() / math.sqrt((a * a).sum() * (b * b).sum())
    write_path = os.path.dirname(segmented_file)

    with open(os.path.join(write_path, template_dict['vbm_qc_filename']),
              'w') as fp:
        fp.write("%3.2f\n" % (covalue))
        fp.close()

    # Flag subjects with <0.90 correlation value.
    if round(covalue, 2) < template_dict['correlation_value']:
        with open(
                os.path.join(write_dir, template_dict['qa_flagged_filename']),
                'w') as fp:
            fp.write("%s\n" % (sub_id))
            fp.close()


def run_pipeline(write_dir,
                 smri_data,
                 reorient,
                 datasink,
                 vbm_preprocess,
                 covars,
                 data_type=None,
                 **template_dict):
    """This function runs pipeline"""

    count_success = 0
    write_dir = write_dir + '/' + template_dict['output_zip_dir']
    error_log = dict()

    for each_sub in smri_data:
        sub_id = (each_sub.split('/')[-1]).split('.')[0]

        try:
            if data_type == 'nifti':
                session = ''
                n1_img = nib.load(each_sub)

            if data_type == 'dicoms':
                session = ''
                vbm_out = os.path.join(write_dir, sub_id, session, 'anat')
                os.makedirs(vbm_out, exist_ok=True)

                # Run DICOM to NIfTI conversion.
                from nipype.interfaces.spm.utils import DicomImport
                import nipype.pipeline.engine as pe
                dcm_nii_convert = pe.Node(
                    interface=DicomImport(), name='converter')
                dcm_nii_convert.inputs.in_files = glob.glob(
                    os.path.join(each_sub, '*'))
                dcm_nii_convert.inputs.output_dir = vbm_out
                with stdchannel_redirected(sys.stderr, os.devnull):
                    dcm_nii_convert.run()
                n1_img = nib.load(glob.glob(os.path.join(vbm_out, '*.nii'))[0])

            vbm_out = os.path.join(write_dir, sub_id, session, 'anat')
            os.makedirs(vbm_out, exist_ok=True)

            if n1_img:
                # DICOM branch already writes NIfTI output.
                if data_type != 'dicoms':
                    nib.save(n1_img, os.path.join(vbm_out, sub_id))

                # Create vbm_spm12 dir under the specific sub-id/anat
                os.makedirs(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname']),
                    exist_ok=True)

                nifti_file = glob.glob(os.path.join(vbm_out, '*.nii'))[0]

                # Edit reorient node inputs
                reorient.node.inputs.in_file = nifti_file
                reorient.node.inputs.out_file = vbm_out + "/" + template_dict[
                    'vbm_output_dirname'] + "/Re.nii"

                # Edit datasink node inputs
                datasink.node.inputs.base_directory = vbm_out

                # Run the nipype pipeline
                with stdchannel_redirected(sys.stderr, os.devnull):
                    vbm_preprocess.run()

                # Smooth modulated images from segmentation node spm.Smooth()
                smooth_images(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname']), **template_dict)

                # Calculate correlation coefficient of swc1*nii to SPM12 TPM.nii
                segmented_file = glob.glob(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname'],
                                 template_dict['qc_nifti']))
                get_corr(segmented_file[0], write_dir, sub_id, **template_dict)

                # Write readme files
                write_readme_files(write_dir, data_type, **template_dict)

                # Convert wc1*.nii to wc1*.png
                label = sub_id + session
                nii_to_image_converter(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname']),
                    label, **template_dict)

        except Exception as e:
            error_log.update({sub_id: str(e)})
            continue

        else:
            count_success = count_success + 1

            if count_success == 1:
                shutil.copy(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname'],
                                 template_dict['display_image_name']),
                    os.path.dirname(write_dir))

        finally:
            remove_tmp_files()

    vbm_spm12_file_output.make_file_output(write_dir, template_dict, covars)

    if os.path.isfile(
            os.path.join(
                os.path.dirname(write_dir),
                template_dict['display_image_name'])):
        # Zip output files.
        shutil.make_archive(
            os.path.join(
                os.path.dirname(write_dir), template_dict['output_zip_dir']),
            'zip', write_dir)

        download_outputs_path = write_dir + '.zip'

        output_message = "VBM preprocessing completed. " + str(
            count_success) + "/" + str(
                len(smri_data)
            ) + " subjects completed successfully." + template_dict[
                'coinstac_display_info']

        preprocessed_percentage = (count_success / len(smri_data)) * 100

        if os.path.isfile(
                os.path.join(write_dir, template_dict['qa_flagged_filename'])):
            qa_percentage = (len(
                open(
                    os.path.join(write_dir,
                                 template_dict['qa_flagged_filename'])).
                readlines()) / len(smri_data)) * 100
            if (qa_percentage <= template_dict['qc_threshold']) or (preprocessed_percentage <= template_dict['qc_threshold']):
                output_message = output_message + template_dict['flag_warning']
        else:
            if (preprocessed_percentage <= template_dict['qc_threshold']):
                output_message = output_message + template_dict['flag_warning']

        if bool(error_log):
            output_message = output_message + " Error log:" + str(error_log)

        with open(
                os.path.join(
                    os.path.dirname(write_dir),
                    template_dict['display_image_name']), "rb") as imageFile:
            encoded_image_str = str(base64.b64encode(imageFile.read()))

        return {
            "output": {
                "message": output_message,
                "download_outputs": download_outputs_path,
                "display": encoded_image_str,
                "covariates": covars
            },
            "cache": {},
            "success": True
        }
    else:
        return {
            "output": {
                "message": " Error log:" + str(error_log)
            },
            "cache": {},
            "success": True
        }
