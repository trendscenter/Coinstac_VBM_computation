#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer runs the pre-processing VBM (Voxel Based Morphometry) pipeline based on the inputs from interface adapter layer
This layer uses entities layer to modify nodes of the pipeline as needed
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


import sys, os, glob, shutil, math, base64, warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
import ujson as json
from bids.grabbids import BIDSLayout
import nibabel as nib
import nipype.pipeline.engine as pe
import numpy as np
from nilearn import plotting

import vbm_entities_layer

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')


def execute_pipeline(bids_dir='',
                     write_dir='',
                     nii_files='',
                     data_type=None,
                     pipeline_opts=None,
                     **template_dict):
    """Runs the pre-processing pipeline on structural T1w scans in BIDS data
        Args:
            bids_dir (string) : Input BIds directory
            write_dir (string): Directory to write outputs
            pipeline_opts ( integer) : Options to change pipeline
            template_dict ( dictionary) : Dictionary that stores all the paths, file names, software locations
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
            After setting up the pipeline here , the pipeline is run

        """

    [reorient, datasink, vbm_preprocess] = create_pipeline_nodes(
        pipeline_opts, **template_dict)

    if data_type == 'bids':
        # Runs the pipeline on each subject, this algorithm runs serially
        layout = BIDSLayout(bids_dir)
        smri_data = layout.get(
            type=template_dict['scan_type'], extensions='.nii.gz')

        return run_pipeline(
            write_dir,
            smri_data,
            reorient,
            datasink,
            vbm_preprocess,
            data_type='bids',
            **template_dict)
    elif data_type == 'nifti':
        # Runs the pipeline on each nifti file, this algorithm runs serially
        smri_data = nii_files
        return run_pipeline(
            write_dir,
            smri_data,
            reorient,
            datasink,
            vbm_preprocess,
            data_type='nifti',
            **template_dict)


def remove_tmp_files():
    """this function removes any tmp files in the docker"""

    for a in glob.glob('/var/tmp/*'):
        os.remove(a)

    for b in glob.glob(os.getcwd() + '/crash*'):
        os.remove(b)

    for c in glob.glob(os.getcwd() + '/tmp*'):
        shutil.rmtree(c, ignore_errors=True)

    for d in glob.glob(os.getcwd() + '/__pycache__'):
        shutil.rmtree(d, ignore_errors=True)

    shutil.rmtree(os.getcwd() + '/vbm_preprocess', ignore_errors=True)

    if os.path.exists(os.getcwd() + '/pyscript.m'):
        os.remove(os.getcwd() + '/pyscript.m')


def write_readme_files(write_dir='', data_type=None, **template_dict):
    """This function writes readme files"""

    # Write a text file with info. on each of the output nifti files
    if data_type == 'bids':
        with open(
                os.path.join(write_dir,
                             template_dict['bids_outputs_manual_name']),
                'w') as fp:
            fp.write(template_dict['bids_outputs_manual_content'])
            fp.close()
    elif data_type == 'nifti':
        with open(
                os.path.join(write_dir,
                             template_dict['nifti_outputs_manual_name']),
                'w') as fp:
            fp.write(template_dict['nifti_outputs_manual_content'])
            fp.close()

    # Write a text file with info. on quality control correlation coefficent
    with open(os.path.join(write_dir, template_dict['qc_readme_name']),
              'w') as fp:
        fp.write(template_dict['qc_readme_content'])
        fp.close()


def nii_to_string_converter(write_dir, label, **template_dict):
    """This function converts nifti to base64 string"""
    import nibabel as nib
    from nilearn import plotting
    import os, base64

    file = os.path.join(write_dir, template_dict['display_nifti'])

    mask = nib.load(file)
    new_data = mask.get_data()
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


#Compute corrcoef
def get_corr(segmented_file, **template_dict):
    """This function computes correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox """

    def extract_data(file):
        a_data = nib.load(file)
        t_data = a_data.get_data()
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


def create_pipeline_nodes(pipeline_opts, **template_dict):
    """This function creates and modifies nodes of the pipeline from entities layer with nipype

        smooth.node.inputs.fwhm: (a list of from 3 to 3 items which are a float or a float)
        3-list of fwhm for each dimension
        This is the size of the Gaussian (in mm) for smoothing the preprocessed data by. This is typically between about 4mm and 12mm.

        segment.node.inputs.channel_info: (a tuple of the form: (a float, a float, a tuple of the
        form: (a boolean, a boolean)))
        A tuple with the following fields:
         - bias regularisation (0-10)
         - FWHM of Gaussian smoothness of bias
         - which maps to save (Corrected, Field) - a tuple of two boolean
        values
    """

    # 1 Reorientation node and settings #
    reorient = vbm_entities_layer.Reorient(**template_dict)

    # 2 Segementation Node and settings #
    segment = vbm_entities_layer.Segment(**template_dict)

    def create_tissue(tpm_path,
                      tissue_id,
                      write_native_maps,
                      write_dartel_maps,
                      write_unmodulated_maps,
                      write_modulated_maps,
                      num_gaussians=None):
        """
        tissue_id is tissue probability image for this class
        1-grey matter
        2-white matter
        3-CSF
        4-bone
        5-soft tissues
        6-air (background)

        write_native_maps,write_dartel_maps- which maps to save [Native, DARTEL] - a tuple of two boolean
        values
        write_unmodulated_maps,z-which maps to save [Unmodulated, Modulated] - a tuple of two
        boolean values

        num_gaussians is the number of Gaussians used to represent the intensity distribution for each type of tissue and can be greater than one.
        In other words, a tissue probability map may be shared by several clusters. The assumption of
        a single Gaussian distribution for each class does not hold for a number of reasons.
        Typical numbers of Gaussians could be 1 for grey matter, 1 for white matter, two for CSF, three for bone, four for other soft tissues and two for air (background).
        """

        NUM_gaussianS = 1

        if num_gaussians is not None:
            num_gaussians = num_gaussians
        else:
            num_gaussians = NUM_gaussianS
        tissue_type = (tpm_path, tissue_id)
        return (tissue_type, num_gaussians, (write_native_maps,
                                             write_dartel_maps),
                (write_unmodulated_maps, write_modulated_maps))

    Tis1 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=1,
        num_gaussians=1,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )

    Tis2 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=2,
        num_gaussians=1,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )

    Tis3 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=3,
        num_gaussians=2,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )

    Tis4 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=4,
        num_gaussians=3,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )

    Tis5 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=5,
        num_gaussians=4,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )

    Tis6 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=6,
        num_gaussians=2,
        write_native_maps=True,
        write_dartel_maps=False,
        write_unmodulated_maps=True,
        write_modulated_maps=True,
    )
    segment.node.inputs.tissues = [Tis1, Tis2, Tis3, Tis4, Tis5, Tis6]

    # 3 Smoothing Node & Settings #
    smooth = vbm_entities_layer.Smooth(**template_dict)

    # 4 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir #
    datasink = vbm_entities_layer.Datasink()

    # 5 Modify Pipeline based on opts to update smoothing fwhm ( full width half maximum ) in mm in x,y,z directions
    if pipeline_opts is not None:
        fwhm = float(pipeline_opts)
        smooth.node.inputs.fwhm = [fwhm] * 3

    ## 6 Create the pipeline/workflow and connect the nodes created above ##
    vbm_preprocess = pe.Workflow(name="vbm_preprocess")

    list_norm_images = vbm_entities_layer.List_Normalized_Images()

    vbm_preprocess.connect([
        create_workflow_input(
            source=reorient.node,
            target=segment.node,
            source_output='out_file',
            target_input='channel_files'),
        create_workflow_input(
            source=segment.node,
            target=list_norm_images.node,
            source_output='normalized_class_images',
            target_input='normalized_class_images'),
        create_workflow_input(
            source=list_norm_images.node,
            target=smooth.node,
            source_output='list_norm_images',
            target_input='in_files'),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='modulated_class_images',
            target_input=template_dict['vbm_output_dirname']),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='native_class_images',
            target_input=template_dict['vbm_output_dirname'] + '.@1'),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='normalized_class_images',
            target_input=template_dict['vbm_output_dirname'] + '.@2'),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='transformation_mat',
            target_input=template_dict['vbm_output_dirname'] + '.@3'),
        create_workflow_input(
            source=smooth.node,
            target=datasink.node,
            source_output='smoothed_files',
            target_input=template_dict['vbm_output_dirname'] + '.@4')
    ])
    return [reorient, datasink, vbm_preprocess]


def create_workflow_input(source, target, source_output, target_input):
    return (source, target, [(source_output, target_input)])


def run_pipeline(write_dir,
                 smri_data,
                 reorient,
                 datasink,
                 vbm_preprocess,
                 data_type=None,
                 **template_dict):
    """This function runs pipeline on the current case"""

    id = 0  # id for assigning sub-id incase of nifti files in txt format
    count_success = 0  # variable for counting how many subjects were successfully run
    write_dir = write_dir + '/' + template_dict['output_zip_dir']  # Store outputs in this directory for zipping the directory
    error_log = dict()  # Dictionary for storing error log for each subject

    for each_sub in smri_data:

        try:

            # Extract subject id and name of nifti file
            if data_type == 'bids':
                sub_id = 'sub-' + each_sub.subject
                session_id = getattr(each_sub, 'session', None)

                if session_id is not None:
                    session = 'ses-' + getattr(each_sub, 'session', None)
                else:
                    session = ''

                nii_output = ((
                    each_sub.filename).split('/')[-1]).split('.gz')[0]
                n1_img = nib.load(each_sub.filename)

            if data_type == 'nifti':
                id = id + 1
                sub_id = 'sub-' + str(id)
                nii_output = ((each_sub).split('/')[-1]).split('.gz')[0]
                n1_img = nib.load(each_sub)

            # Directory in which vbm outputs will be written
            vbm_out = os.path.join(write_dir, sub_id, session, 'anat')

            # Create output dir for sub_id
            os.makedirs(vbm_out, exist_ok=True)

            nifti_file = os.path.join(vbm_out, nii_output)

            if n1_img:
                nib.save(n1_img, os.path.join(vbm_out, nii_output))

                # Create vbm_spm12 dir under the specific sub-id/anat
                os.makedirs(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname']),
                    exist_ok=True)

                # Edit reorient node inputs
                reorient.node.inputs.in_file = nifti_file
                reorient.node.inputs.out_file = vbm_out + "/" + template_dict['vbm_output_dirname'] + "/Re.nii"

                # Edit datasink node inputs
                datasink.node.inputs.base_directory = vbm_out

                # Run the nipype pipeline
                with stdchannel_redirected(sys.stderr, os.devnull):
                    vbm_preprocess.run()

                # Calculate correlation coefficient of swc1*nii to SPM12 TPM.nii
                segmented_file = glob.glob(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname'],
                                 template_dict['qc_nifti']))
                get_corr(segmented_file[0], **template_dict)

                # Write readme files
                write_readme_files(write_dir, data_type, **template_dict)

                label = sub_id + session
                nii_to_string_converter(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname']),
                    label, **template_dict)

        except Exception as e:
            # If fails raise the exception,print exception error
            #sys.stderr.write(str(e))
            error_log.update({sub_id: str(e)})
            continue

        else:
            # If the try block succeeds, increase the count
            count_success = count_success + 1

            if count_success == 1:
                shutil.copy(
                    os.path.join(vbm_out, template_dict['vbm_output_dirname'],
                                 template_dict['display_image_name']),
                    os.path.dirname(write_dir))

        finally:
            remove_tmp_files()

    #Zip output files
    shutil.make_archive(
        os.path.join(
            os.path.dirname(write_dir), template_dict['output_zip_dir']),
        'zip', write_dir)

    #shutil.rmtree(write_dir, ignore_errors=True)

    download_outputs_path = write_dir + '.zip'

    output_message = "VBM preprocessing completed. " + str(
        count_success) + "/" + str(
            len(smri_data)) + " subjects" + " completed successfully."

    if bool(error_log):
        output_message = output_message + " Error log:" + str(error_log)

    with open(
            os.path.join(
                os.path.dirname(write_dir),
                template_dict['display_image_name']), "rb") as imageFile:
        encoded_image_str = base64.b64encode(imageFile.read())

    return json.dumps({
        "output": {
            "message": output_message,
            "download_outputs": download_outputs_path,
            "display": encoded_image_str
        },
        "cache": {},
        "success": True
    })
