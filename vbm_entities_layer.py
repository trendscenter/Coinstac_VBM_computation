#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer defines the nodes of pre-processing pipeline and correlation computation function
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.spm as spm
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function

import nibabel as nb
import numpy as np
import math,os

## 1 Reorientation node & settings ##
class Reorient:
	def __init__(self):
		self.node = pe.Node(interface=spm.ApplyTransform(), name='reorient')

## 2 Segementation Node and settings ##
class Segment:
	def __init__(self):
		self.node = pe.Node(interface=spm.NewSegment(), name='segmentation')


def transform_list(normalized_class_images):
		return [each[0] for each in normalized_class_images]

class List_normalized_images:
	def __init__(self):
		self.node=pe.Node(interface=Function(
			input_names='normalized_class_images',
			output_names='list_norm_images',
			function=transform_list
		), name='List_normalized_images')

## 4 Smoothing Node & Settings ##
class Smooth:
	def __init__(self):
		self.node = pe.Node(interface=spm.Smooth(), name='smoothing')

## 5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir ##
class Datasink:
	def __init__(self):
		self.node = pe.Node(interface=DataSink(), name='sinker')

#Explore corrcoef
def get_corr(template,segmented_file,corr_file_name):
	"""This function computes correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox """

	def get_data(file):
		a_data = nb.load(file)
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

	cstn_data = get_data(template)
	cre_data = get_data(segmented_file)
	indices = np.logical_and(cstn_data!=0,cre_data!=0)
	fcstn_data, fcre_data = cstn_data[indices], cre_data[indices]

	a=fcstn_data-np.mean(fcstn_data)
	b=fcre_data-np.mean(fcre_data)
	covalue = (a * b).sum() / math.sqrt((a * a).sum() * (b * b).sum())
	write_path= os.path.dirname(segmented_file)

	with open(os.path.join(write_path,corr_file_name),'w') as fp:
		fp.write("%3.2f\n"%(covalue))
		fp.close()

