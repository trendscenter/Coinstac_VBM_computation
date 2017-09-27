#!/usr/bin/env python

import nibabel as nb
import numpy as np
import math

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

def mean2(x):
        y = np.sum(x) / np.size(x);
        return y

def corr2(a,b):
        a = a - mean2(a)
        b = b - mean2(b)
        r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum())
        return r

def get_corr(template,segmented_file):
    cstn_data = get_data(template)
    cre_data = get_data(segmented_file)
    indices = np.logical_and(cstn_data!=0,cre_data!=0)
    fcstn_data, fcre_data = cstn_data[indices], cre_data[indices]   
    coval = corr2(fcstn_data,fcre_data)
    write_path=segmented_file.rsplit('/',1)[0]
    with open(write_path+'/vbm_corr_value.txt','w') as fp:
	    fp.write("%3.2f\n"%(coval.item()))
    fp.close()
    return "%3.2f" % (coval.item())

