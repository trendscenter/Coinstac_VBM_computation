#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""COINSTAC entrypoint for VBM preprocessing orchestration."""
import os
import sys
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

# Stop printing nipype.workflow info to stdout.
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

# Store all paths/options/output names in a single config dictionary.
from vbm_config import template_dict
from vbm_args import args_parser
from vbm_data import data_parser
from vbm_runtime import (
    convert_reorientparams_save_to_mat_script,
    software_check,
    stdchannel_redirected,
)


def start(args):
    # Validate SPM runtime availability before any pipeline work.
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm_check = software_check(template_dict)
    if spm_check != template_dict['spm_version']:
        raise EnvironmentError("spm unable to start in vbm docker")

    # Parse args and update runtime template settings.
    args_output = args_parser(args, template_dict)
    if args_output is not None:
        return args_output

    # Write reorientation transform matrix consumed by SPM.
    convert_reorientparams_save_to_mat_script(template_dict)

    # Run standalone VBM preprocessing.
    return data_parser(args, template_dict)
