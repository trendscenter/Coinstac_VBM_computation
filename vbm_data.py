import os

import vbm_standalone_use_cases_layer


def data_parser(args, template_dict):
    """Parse input data and run preprocessing pipeline."""
    # We are only concerned with pre-processing.
    data = [args['state']['baseDirectory'] + '/' + file_names for file_names in args['input']['covariates']]
    write_dir = args['state']['outputDirectory']

    for path in data:
        if not os.path.isfile(path):
            raise Exception("File does not exist or can't be read: " + str(path))
    if not os.access(write_dir, os.W_OK):
        raise Exception("Output write permissions denied: " + str(write_dir))

    covariates = args['input']['covariates']
    return vbm_standalone_use_cases_layer.setup_pipeline(
        data=data,
        write_dir=write_dir,
        covars=covariates,
        data_type='nifti',
        **template_dict)
