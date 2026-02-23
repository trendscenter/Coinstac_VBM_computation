import os

import nibabel as nib


def args_parser(args, template_dict):
    """Extract options from input arguments into template_dict."""
    input_args = args.get('input', {})

    if 'options_smoothing_x_mm' in input_args:
        template_dict['FWHM_SMOOTH'][0] = float(input_args['options_smoothing_x_mm'])
    if 'options_smoothing_y_mm' in input_args:
        template_dict['FWHM_SMOOTH'][1] = float(input_args['options_smoothing_y_mm'])
    if 'options_smoothing_z_mm' in input_args:
        template_dict['FWHM_SMOOTH'][2] = float(input_args['options_smoothing_z_mm'])
    if 'options_smoothing_implicit_masking' in input_args:
        template_dict['implicit_masking'] = input_args['options_smoothing_implicit_masking']

    if input_args.get('standalone'):
        template_dict['standalone'] = input_args['standalone']
    else:
        template_dict['covariates'] = input_args.get('covariates', [])
        template_dict['regression_data'] = input_args.get('data', [])

    if 'options_reorient_params_x_mm' in input_args:
        template_dict['options_reorient_params_x_mm'] = float(input_args['options_reorient_params_x_mm'])
    if 'options_reorient_params_y_mm' in input_args:
        template_dict['options_reorient_params_y_mm'] = float(input_args['options_reorient_params_y_mm'])
    if 'options_reorient_params_z_mm' in input_args:
        template_dict['options_reorient_params_z_mm'] = float(input_args['options_reorient_params_z_mm'])
    if 'options_reorient_params_pitch' in input_args:
        template_dict['options_reorient_params_pitch'] = float((input_args['options_reorient_params_pitch']))
    if 'options_reorient_params_roll' in input_args:
        template_dict['options_reorient_params_roll'] = float((input_args['options_reorient_params_roll']))
    if 'options_reorient_params_yaw' in input_args:
        template_dict['options_reorient_params_yaw'] = float((input_args['options_reorient_params_yaw']))
    if 'options_reorient_params_x_scaling' in input_args:
        template_dict['options_reorient_params_x_scaling'] = float(input_args['options_reorient_params_x_scaling'])
    if 'options_reorient_params_y_scaling' in input_args:
        template_dict['options_reorient_params_y_scaling'] = float(input_args['options_reorient_params_y_scaling'])
    if 'options_reorient_params_z_scaling' in input_args:
        template_dict['options_reorient_params_z_scaling'] = float(input_args['options_reorient_params_z_scaling'])
    if 'options_reorient_params_x_affine' in input_args:
        template_dict['options_reorient_params_x_affine'] = float(input_args['options_reorient_params_x_affine'])
    if 'options_reorient_params_y_affine' in input_args:
        template_dict['options_reorient_params_y_affine'] = float(input_args['options_reorient_params_y_affine'])
    if 'options_reorient_params_z_affine' in input_args:
        template_dict['options_reorient_params_z_affine'] = float(input_args['options_reorient_params_z_affine'])

    if 'options_BIAS_REGULARISATION' in input_args:
        template_dict['BIAS_REGULARISATION'] = float(input_args['options_BIAS_REGULARISATION'])

    if 'options_FWHM_GAUSSIAN_SMOOTH_BIAS' in input_args:
        template_dict['FWHM_GAUSSIAN_SMOOTH_BIAS'] = input_args['options_FWHM_GAUSSIAN_SMOOTH_BIAS']

    if 'options_affine_regularization' in input_args:
        template_dict['affine_regularization'] = input_args['options_affine_regularization']

    if 'options_warping_regularization' in input_args:
        if len(input_args['options_warping_regularization']) == 5:
            template_dict['warping_regularization'] = input_args['options_warping_regularization']

    if 'options_sampling_distance' in input_args:
        template_dict['sampling_distance'] = float(input_args['options_sampling_distance'])

    if 'options_mrf_weighting' in input_args:
        template_dict['mrf_weighting'] = float(input_args['options_mrf_weighting'])

    if 'options_cleanup' in input_args:
        template_dict['cleanup'] = int(input_args['options_cleanup'])

    if 'regression_resample_voxel_size' in input_args:
        voxel_size = float(input_args['regression_resample_voxel_size'])
        template_dict['regression_resample_voxel_size'] = (voxel_size, voxel_size, voxel_size)

    if 'registration_template' in input_args:
        registration_template = input_args['registration_template']
        if os.path.isfile(registration_template) and (
            str(nib.load(template_dict['tpm_path']).shape)
            == str(nib.load(registration_template).shape)
        ):
            template_dict['tpm_path'] = registration_template
        else:
            return {
                "output": {
                    "message": "Non-standard Registration template "
                },
                "cache": {},
                "success": True
            }
