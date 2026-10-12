#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""COINSTAC entrypoint for VBM preprocessing orchestration."""
import copy
import os
import sys
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

# Stop printing nipype.workflow info to stdout.
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

# Store all paths/options/output names in a single config dictionary.
from vbm.vbm_config import template_dict
from vbm.vbm_args import args_parser
from vbm.vbm_data import data_parser
from vbm.vbm_runtime import (
    convert_reorientparams_save_to_mat_script,
    software_check,
    stdchannel_redirected,
)


def start(args):
    # Use per-run mutable config to avoid leaking state across invocations.
    runtime_template = copy.deepcopy(template_dict)

    # Validate SPM runtime availability before any pipeline work.
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm_check = software_check(runtime_template)
    if spm_check != runtime_template['spm_version']:
        raise EnvironmentError("spm unable to start in vbm docker")

    # Parse args and update runtime template settings.
    args_output = args_parser(args, runtime_template)
    if args_output is not None:
        return args_output

    # Write reorientation transform matrix consumed by SPM.
    convert_reorientparams_save_to_mat_script(runtime_template)

    # Run standalone VBM preprocessing.
    return data_parser(args, runtime_template)
