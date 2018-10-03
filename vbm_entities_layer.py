#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer defines the nodes of pre-processing pipeline and correlation computation function
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.spm as spm
spm.terminal_output = 'file'
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')


## 1 Reorientation node & settings ##
class Reorient:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.ApplyTransform(), name='reorient')
        self.node.inputs.mat = template_dict['transf_mat_path']
        self.node.inputs.paths = template_dict['spm_path']


## 2 Segementation Node and settings ##
class Segment:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.NewSegment(), name='segmentation')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.channel_info = (
            template_dict['BIAS_REGULARISATION'],
            template_dict['FWHM_GAUSSIAN_SMOOTH_BIAS'], (False, False))


def transform_list(normalized_class_images):
    return [each[0] for each in normalized_class_images]


class List_Normalized_Images:
    def __init__(self):
        self.node = pe.Node(
            interface=Function(
                input_names='normalized_class_images',
                output_names='list_norm_images',
                function=transform_list),
            name='List_normalized_images')


## 4 Smoothing Node & Settings ##
class Smooth:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.Smooth(), name='smoothing')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.fwhm = template_dict['FWHM_SMOOTH']


## 5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir ##
class Datasink:
    def __init__(self):
        self.node = pe.Node(interface=DataSink(), name='sinker')
