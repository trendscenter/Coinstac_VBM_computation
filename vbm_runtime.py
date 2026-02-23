import contextlib
import os
import sys
import traceback

import numpy as np
import scipy.io
import spm_matrix as s
from nipype.interfaces import spm


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """Temporarily redirect stdout/stderr to a destination file."""
    oldstdchannel = None
    dest_file = None
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


def software_check(template_dict):
    """Return the SPM standalone version configured in the runtime."""
    spm.SPMCommand.set_mlab_paths(
        matlab_cmd=template_dict['matlab_cmd'], use_mcr=True)
    return spm.SPMCommand().version


def convert_reorientparams_save_to_mat_script(template_dict):
    """Write reorientation parameters to transform.mat for SPM."""
    try:
        pi = 22 / 7
        scipy.io.savemat(
            '/computation/transform.mat',
            mdict={
                'M': np.around(
                    s.spm_matrix([
                        template_dict['options_reorient_params_x_mm'],
                        template_dict['options_reorient_params_y_mm'],
                        template_dict['options_reorient_params_z_mm'],
                        template_dict['options_reorient_params_pitch'] * (pi / 180),
                        template_dict['options_reorient_params_roll'] * (pi / 180),
                        template_dict['options_reorient_params_yaw'] * (pi / 180),
                        template_dict['options_reorient_params_x_scaling'],
                        template_dict['options_reorient_params_y_scaling'],
                        template_dict['options_reorient_params_z_scaling'],
                        template_dict['options_reorient_params_x_affine'],
                        template_dict['options_reorient_params_y_affine'],
                        template_dict['options_reorient_params_z_affine'],
                    ], 1),
                    decimals=4)[0]
            })
    except Exception as e:
        sys.stderr.write(
            'Unable to convert reorientation params to transform.mat Error_log:'
            + str(e)
            + str(traceback.format_exc()))
