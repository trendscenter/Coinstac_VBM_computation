"""Shared SPM workflow construction for both VBM use cases."""

import glob
import os
import sys

import nipype.pipeline.engine as pe

from . import vbm_entities_layer
from .pipeline_common import create_workflow_input, stdchannel_redirected


def create_pipeline_nodes(**template_dict):
    """This function creates and modifies nodes of the pipeline from entities layer with nipype
    """

    #  Reorientation node and settings #
    reorient = vbm_entities_layer.Reorient(**template_dict)

    #  Segementation Node and settings #
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

        NUM_gaussians = 1

        if num_gaussians is not None:
            num_gaussians = num_gaussians
        else:
            num_gaussians = NUM_gaussians
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

    # Lists normalized images
    list_norm_images = vbm_entities_layer.List_Normalized_Images()

    #  Smoothing Node & Settings #
    smooth = vbm_entities_layer.Smooth(**template_dict)

    #  Datsink Node that collects segmented, smoothed files and writes to temp_write_dir #
    datasink = vbm_entities_layer.Datasink()

    # Create the pipeline/workflow and connect the nodes created above #
    vbm_preprocess = pe.Workflow(name="vbm_preprocess")

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


def smooth_images(write_dir, **template_dict):
    """This function runs smoothing on input images. Ex: modulated images"""
    from nipype.interfaces import spm
    from nipype.interfaces.io import DataSink
    smooth = pe.Node(interface=spm.Smooth(), name='smooth')
    smooth.inputs.paths = template_dict['spm_path']
    smooth.inputs.implicit_masking = template_dict['implicit_masking']
    smooth.inputs.in_files = glob.glob(os.path.join(write_dir, 'mwc*.nii'))
    smooth.inputs.fwhm = template_dict['FWHM_SMOOTH']
    vbm_smooth_modulated_images = pe.Workflow(
        name="vbm_smooth_modulated_images")
    datasink = pe.Node(interface=DataSink(), name='datasink')
    datasink.inputs.base_directory = write_dir
    vbm_smooth_modulated_images.connect([(smooth, datasink, [('smoothed_files',
                                                              write_dir)])])
    with stdchannel_redirected(sys.stderr, os.devnull):
        vbm_smooth_modulated_images.run()
