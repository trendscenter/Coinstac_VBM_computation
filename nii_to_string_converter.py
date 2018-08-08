#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def nii_to_string_converter(write_dir, **template_dict):

    import nibabel as nib
    from nilearn import plotting
    import os, base64

    file = os.path.join(write_dir, template_dict['display_nifti'])

    mask = nib.load(file)
    new_data = mask.get_data()
    clipped_img = nib.Nifti1Image(new_data, mask.affine, mask.header)

    plotting.plot_anat(
        clipped_img,
        cut_coords=(0, 0, 0),
        annotate=False,
        draw_cross=False,
        output_file=os.path.join(write_dir,
                                 template_dict['display_pngimage_name']),
        display_mode='ortho',
        title=template_dict['display_image_name'],
        colorbar=False)

    byte_string = (base64.b64encode(
        open(
            os.path.join(write_dir, template_dict['display_pngimage_name']),
            "rb").read()))[1:]

    return byte_string
