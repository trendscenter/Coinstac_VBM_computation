import os

import vbm_standalone_use_cases_layer


def data_parser(args, template_dict):
    """Parse input data and run preprocessing pipeline."""
    # We are only concerned with pre-processing.
    base_dir = args['state']['baseDirectory']
    covariates = args['input']['covariates']
    data = [base_dir + '/' + file_name for file_name in covariates]
    write_dir = args['state']['outputDirectory']

    for path in data:
        if not os.path.isfile(path):
            raise Exception("File does not exist or can't be read: " + str(path))
    if not os.access(write_dir, os.W_OK):
        raise Exception("Output write permissions denied: " + str(write_dir))

    return vbm_standalone_use_cases_layer.setup_pipeline(
        data=data,
        write_dir=write_dir,
        covars=covariates,
        data_type='nifti',
        **template_dict)
