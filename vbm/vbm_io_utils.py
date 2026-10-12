"""Small IO helpers shared by VBM pipeline modules."""

import glob
import os
import shutil


def remove_tmp_files():
    """Remove temporary workflow folders in the working directory."""
    for path in glob.glob(os.path.join(os.getcwd(), "tmp*")):
        shutil.rmtree(path, ignore_errors=True)


def write_readme_files(write_dir='', data_type=None, **template_dict):
    """Write output description and QC readme files."""
    if data_type == 'nifti':
        with open(os.path.join(write_dir, template_dict['outputs_manual_name']), 'w') as fp:
            fp.write(template_dict['nifti_outputs_manual_content'])
    elif data_type == 'dicoms':
        with open(os.path.join(write_dir, template_dict['outputs_manual_name']), 'w') as fp:
            fp.write(template_dict['dicoms_outputs_manual_content'])

    with open(os.path.join(write_dir, template_dict['qc_readme_name']), 'w') as fp:
        fp.write(template_dict['qc_readme_content'])
