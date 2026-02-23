import os


# Central runtime/config dictionary used across parsing and execution modules.
template_dict = {
    # Runtime/software paths.
    'spm_version': '12.7771',
    'matlab_cmd': '/opt/spm12/run_spm12.sh /opt/mcr/R2022b script',
    'spm_path': '/opt/spm12/fsroot',
    'tpm_path': '/opt/spm12/fsroot/spm/spm12/tpm/TPM.nii',
    'transf_mat_path': os.path.join('/computation', 'transform.mat'),

    # Input/data mode.
    'scan_type': 'T1w',
    'standalone': False,
    'covariates': list(),
    'regression_data': list(),
    'regression_file_input_type': 'swc1',
    'regression_dir_name': 'regression_input_files',
    'regression_file': 'swc1Re.nii',
    'regression_resample_voxel_size': None,
    'regression_resample_method': 'Li',

    # Reorientation + smoothing.
    'FWHM_SMOOTH': [10.0, 10.0, 10.0],
    'bounding_box': '',
    'options_reorient_params_x_mm': 0,
    'options_reorient_params_y_mm': 0,
    'options_reorient_params_z_mm': 0,
    'options_reorient_params_pitch': 0,
    'options_reorient_params_roll': 0,
    'options_reorient_params_yaw': 0,
    'options_reorient_params_x_scaling': 1,
    'options_reorient_params_y_scaling': 1,
    'options_reorient_params_z_scaling': 1,
    'options_reorient_params_x_affine': 0,
    'options_reorient_params_y_affine': 0,
    'options_reorient_params_z_affine': 0,

    # Segmentation options.
    'BIAS_REGULARISATION': 0.0001,
    'FWHM_GAUSSIAN_SMOOTH_BIAS': 60,
    'affine_regularization': 'mni',
    'warping_regularization': [0, 1e-3, 0.5, 0.05, 0.2],
    'sampling_distance': 3.0,
    'mrf_weighting': 1.0,
    'cleanup': 1,
    'implicit_masking': False,

    # Output names and quality-control settings.
    'correlation_value': 0.90,
    'vbm_output_dirname': 'vbm_spm12',
    'output_zip_dir': 'vbm_outputs',
    'log_filename': 'vbm_log.txt',
    'qa_flagged_filename': 'QA_flagged_subjects.txt',
    'display_image_name': 'wc1Re.png',
    'display_pngimage_name': 'Gray matter (Normalized)',
    'cut_coords': (0, 0, 0),
    'display_nifti': 'wc1Re.nii',
    'qc_nifti': 'swc1*nii',
    'vbm_append_string': 'vbm_prepoc',
    'qc_threshold': 70,
    'vbm_qc_filename': 'vbm_corr_value.txt',
    'outputs_manual_name': 'outputs_description.txt',
    'coinstac_display_info': (
        'Please read outputs_description.txt for description of pre-processed output files and '
        'quality_control_readme.txt for quality control measurement.'
        'These files are placed under the pre-processed data.'
    ),
    'flag_warning': (
        ' QC warning: Atleast 30% of input data did not pass QA or could not be pre-processed, '
        'please check the data, vbm_log.txt and QA_flagged_subjects.txt'
    ),

    # Manual/readme text content.
    'bids_outputs_manual_content': (
        "Prefixes descriptions for segmented images:c1-Gray matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
        "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) Default:[10 10 10]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf"
    ),
    'nifti_outputs_manual_content': (
        "subID-1,subID-2,subID-* denotes each nifti file with respect to the order of the nifti paths given"
        "\nPrefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
        "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf"
    ),
    'dicoms_outputs_manual_content': (
        "subID-1,subID-2,subID-* denotes each session's dicom directory with respect to the order of the dicom directory paths given"
        "\nPrefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air(background)"
        "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf"
    ),
    'qc_readme_name': 'quality_control_readme.txt',
    'qc_readme_content': (
        "In each subject's anat/vbm_spm12 directory,vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12 adult brain tissue probability maps (TPM.nii) file from SPM12 toolbox"
        "\nIf your subjects are kids/adolescent, take this into considersation as correlation value may not pass threshold. For an adult study, its suggested that scans with correlation value <=0.91 "
        "\nshould be manually looked into for possible reorientation to the ac-pc line. Subjects that do not pass this QA metric will be saved in vbm_outputs/QA_flagged_subjects.txt"
    ),
}
