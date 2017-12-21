"""************************** corr.py ****************"""
""" This script is called by run_vbm_bids.py """
""" This script computes the correlation value between wc1*nii from the vbm output to the gray matter template from spm12"""

import nibabel as nb
import numpy as np
import math
import os

VBM_CORR_VALUE_FILE_NAME='vbm_corr_value.txt'

def get_data(path):
    a_data = nb.load(path)
    t_data = a_data.get_data()
    if len(t_data.shape) == 4:
        tx,ty,tz,im = t_data.shape
        st_data = t_data[:,:,:,0]
    else:
        tx,ty,tz = t_data.shape
        st_data = t_data
    stn_data = st_data.reshape((tx*ty*tz,1), order = 'F')
    stn_data = np.nan_to_num(stn_data)
    return stn_data


def get_corr(template,segmented_file):
    cstn_data = get_data(template)
    cre_data = get_data(segmented_file)
    indices = np.logical_and(cstn_data!=0,cre_data!=0)
    fcstn_data, fcre_data = cstn_data[indices], cre_data[indices]   
    covalue = (np.corrcoef(fcstn_data,fcre_data)).item()
    write_path= os.path.dirname(segmented_file)
    with open(os.path.join(write_path, VBM_CORR_VALUE_FILE_NAME),'w') as fp:
	    fp.write("%3.2f\n"%(covalue))
    fp.close()
    return "%3.2f" % (covalue)

