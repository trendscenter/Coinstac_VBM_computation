#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer runs the pre-processing VBM (Voxel Based Morphometry) pipeline based on the inputs from interface adapter layer
This layer uses entities layer to modify nodes of the pipeline as needed
"""

import sys,os,glob,shutil
import ujson as json
from bids.grabbids import BIDSLayout
import nibabel as nib
import nipype.pipeline.engine as pe


#Import Enities layer
import vbm_entities_layer

def run_bids_pipeline( bids_dir='',write_dir='',pipeline_opts=None, **template_dict):
    """This function runs the bids use case1"""

    [reorient, datasink, vbm_preprocess]=create_pipeline_nodes(pipeline_opts, **template_dict)

    # Runs the pipeline on each subject, this algorithm runs serially
    layout = BIDSLayout(bids_dir)
    smri_data = layout.get(extensions=template_dict['scan_type']+'.nii.gz')
    return run_pipeline(write_dir, smri_data, reorient, datasink, vbm_preprocess,flag='bids', **template_dict)

def run_nifti_pipeline( nii_files='',write_dir='',pipeline_opts=None, **template_dict):
    """This function runs the nifti use case2"""

    [reorient, datasink, vbm_preprocess]=create_pipeline_nodes(pipeline_opts, **template_dict)

    # Runs the pipeline on each nifti file, this algorithm runs serially
    smri_data = nii_files
    return run_pipeline(write_dir, smri_data, reorient, datasink, vbm_preprocess,flag='nifti', **template_dict)


def remove_tmp_files():
    """this function removes any tmp files in the docker"""

    if (os.path.exists('/var/tmp')):
        shutil.rmtree('/var/tmp/*', ignore_errors=True)

    for c in glob.glob(os.getcwd() + '/crash*'):
        os.remove(c)

    for f in glob.glob(os.getcwd() + '/tmp*'):
        shutil.rmtree(f, ignore_errors=True)

    for f in glob.glob(os.getcwd() + '/__pycache__'):
        shutil.rmtree(f, ignore_errors=True)

    if os.path.exists(os.getcwd() + '/vbm_preprocess'):
        shutil.rmtree(os.getcwd() + '/vbm_preprocess', ignore_errors=True)

    if os.path.exists(os.getcwd() + '/pyscript.m'):
        os.remove(os.getcwd() + '/pyscript.m')

def write_readme_files(vbm_out,**template_dict):
    """This function writes readme files"""

    # Write a text file with info. on each of the output nifti files
    with open(os.path.join(vbm_out,template_dict['vbm_output_dirname'],template_dict['outputs_manual_name']), 'w') as fp:
        fp.write(template_dict['outputs_manual_content'])
        fp.close()

    # Write a text file with info. on quality control correlation coefficent
    with open(os.path.join(vbm_out,template_dict['vbm_output_dirname'],template_dict['qc_readme_name']), 'w') as fp:
        fp.write(template_dict['qc_readme_content'])
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
    reorient = vbm_entities_layer.Reorient()
    reorient.node.inputs.mat = template_dict['transf_mat_path']
    reorient.node.inputs.paths = template_dict['spm_path']

    # 2 Segementation Node and settings #
    segment = vbm_entities_layer.Segment()
    segment.node.inputs.paths = template_dict['spm_path']

    # transform the foll. into readable function and assign names to values
    segment.node.inputs.channel_info = (template_dict['BIAS_REGULARISATION'], template_dict['FWHM_GUASSIAN_SMOOTH_BIAS'], (False, False))



    def create_tissue(tpm_path, tissue_id, w, x, y, z, num_guassians=None):
        """
        tissue_id is tissue probability image for this class
        1-grey matter
        2-white matter
        3-CSF
        4-bone
        5-soft tissues
        6-air (background)

        w,x- which maps to save [Native, DARTEL] - a tuple of two boolean
        values
         y,z-which maps to save [Unmodulated, Modulated] - a tuple of two
        boolean values

        num_guassians is the number of Gaussians used to represent the intensity distribution for each type of tissue and can be greater than one.
        In other words, a tissue probability map may be shared by several clusters. The assumption of
        a single Gaussian distribution for each class does not hold for a number of reasons.
        Typical numbers of Gaussians could be 1 for grey matter, 1 for white matter, two for CSF, three for bone, four for other soft tissues and two for air (background).
        """

        NUM_GUASSIANS = 1

        if num_guassians is not None:
            num_guassians = num_guassians
        else:
            num_guassians = NUM_GUASSIANS
        tissue_type = (tpm_path, tissue_id)
        return (tissue_type, num_guassians, (w, x), (y, z))

    Tis1 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=1,
        num_guassians=1,
        w=True,
        x=False,
        y=True,
        z=True,
    )

    Tis2 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=2,
        num_guassians=1,
        w=True,
        x=False,
        y=True,
        z=True,
    )

    Tis3 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=3,
        num_guassians=2,
        w=True,
        x=False,
        y=True,
        z=True,
    )

    Tis4 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=4,
        num_guassians=3,
        w=True,
        x=False,
        y=True,
        z=True,
    )

    Tis5 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=5,
        num_guassians=4,
        w=True,
        x=False,
        y=True,
        z=True,
    )

    Tis6 = create_tissue(
        tpm_path=template_dict['tpm_path'],
        tissue_id=6,
        num_guassians=2,
        w=True,
        x=False,
        y=True,
        z=True,
    )
    segment.node.inputs.tissues = [Tis1, Tis2, Tis3, Tis4, Tis5, Tis6]

    # 3 Smoothing Node & Settings #
    smooth = vbm_entities_layer.Smooth()
    smooth.node.inputs.paths = template_dict['spm_path']
    smooth.node.inputs.fwhm = template_dict['FWHM_SMOOTH']

    # 4 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir #
    datasink = vbm_entities_layer.Datasink()

    # 5 Modify Pipeline based on opts to update smoothing fwhm ( full width half maximum ) in mm in x,y,z directions
    if pipeline_opts is not None:
        fwhm = float(pipeline_opts['fwhm'])
        smooth.node.inputs.fwhm = [fwhm] * 3


    ## 6 Create the pipeline/workflow and connect the nodes created above ##
    vbm_preprocess = pe.Workflow(name="vbm_preprocess")

    list_norm_images=vbm_entities_layer.List_normalized_images()

    vbm_preprocess.connect([
        create_workflow_input(
            source=reorient.node,
            target=segment.node,
            source_output='out_file',
            target_input='channel_files'
        ),
        create_workflow_input(
            source=segment.node,
            target=list_norm_images.node,
            source_output='normalized_class_images',
            target_input='normalized_class_images'
        ),
        create_workflow_input(
            source=list_norm_images.node,
            target=smooth.node,
            source_output='list_norm_images',
            target_input='in_files'
        ),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='modulated_class_images',
            target_input=template_dict['vbm_output_dirname']
        ),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='native_class_images',
            target_input=template_dict['vbm_output_dirname'] + '.@1'
        ),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='normalized_class_images',
            target_input=template_dict['vbm_output_dirname'] + '.@2'
        ),
        create_workflow_input(
            source=segment.node,
            target=datasink.node,
            source_output='transformation_mat',
            target_input=template_dict['vbm_output_dirname'] + '.@3'
        ),
        create_workflow_input(
            source=smooth.node,
            target=datasink.node,
            source_output='smoothed_files',
            target_input=template_dict['vbm_output_dirname'] + '.@4'
        )
    ])
    return [reorient, datasink, vbm_preprocess]


def create_workflow_input(source, target, source_output, target_input):
    return (source, target, [(source_output, target_input)])



def run_pipeline(write_dir, smri_data, reorient, datasink, vbm_preprocess, flag=None,**template_dict):
    """This function runs pipeline on the current case"""

    id = 0  # id for assigning sub-id incase of nifti files in txt format
    count_success = 0  # variable for counting how many subjects were successfully run

    # create dirs array to store output directories where vbm spm12 data is written to
    dirs = []

    # create wc1files array to store paths to wc1 files for each subject
    wc1files = []

    for each_sub in smri_data:


        try:

            # Extract subject id and name of nifti file
            if flag == 'bids':
                sub_id = 'sub-' + each_sub.subject
                nii_output = ((each_sub.filename).split('/')[-1]).split('.gz')[0]
                n1_img = nib.load(each_sub.filename)

            if flag == 'nifti':
                id = id + 1
                sub_id = 'sub-' + str(id)
                nii_output = ((each_sub).split('/')[-1]).split('.gz')[0]
                n1_img = nib.load(each_sub)

            # Directory in which vbm outputs will be written to
            vbm_out = os.path.join(write_dir, sub_id, 'anat')

            # Create output dir for sub_id
            os.makedirs(vbm_out, exist_ok=True)

            nifti_file = os.path.join(vbm_out, nii_output)

            if n1_img:
                nib.save(n1_img, os.path.join(vbm_out, nii_output))

                # Create vbm_spm12 dir under the specific sub-id/anat
                os.makedirs(os.path.join(vbm_out,template_dict['vbm_output_dirname']), exist_ok=True)

                # Edit reorient node inputs
                reorient.node.inputs.in_file = nifti_file
                reorient.node.inputs.out_file = vbm_out + "/"+template_dict['vbm_output_dirname']+"/Re.nii"

                # Edit datasink node inputs
                datasink.node.inputs.base_directory = vbm_out

                # Run the nipype pipeline
                vbm_preprocess.run()

                # Append vbm output dirs and wc*nii files for json output
                dirs.append(os.path.join(vbm_out ,template_dict['vbm_output_dirname']))
                wc1files.append(glob.glob(os.path.join(vbm_out,template_dict['vbm_output_dirname'],"wc1*nii"))[0])

                # Calculate correlation coefficient of swc1*nii to SPM12 TPM.nii
                segmented_file = glob.glob(os.path.join(vbm_out,template_dict['vbm_output_dirname'],"swc1*nii"))
                vbm_entities_layer.get_corr(template_dict['tpm_path'], segmented_file[0],template_dict['vbm_qc_filename'])

                # Write readme files
                write_readme_files(vbm_out,**template_dict)

        except Exception as e:
            # If fails raise the exception,print exception err string to stderr
            sys.stderr.write(str(e))
            continue

        else:
            # If the try block succeeds, increase the count
            count_success = count_success + 1

        finally:
            remove_tmp_files()

    '''
    Calculate how many nifti's successfully got run through the pipeline, this may help in colloborative projects
    where some of the projects may have low quality data
    '''
    completed_percent = (count_success/len(smri_data))*100;


    '''
    Status=True means program finished execution , despite the success or failure of the code
    This is to indicate to coinstac that program finished execution
    '''

    status = True
    return json.dumps({"output": {"success": status, "vbmdirs": dirs, "wc1files": wc1files, "complete%":completed_percent}})



