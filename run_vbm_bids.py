# !/usr/bin/env python

#Example run of the code- python3 run_vbm_bids.py '{"input_bids_dir":"/root/data/bids_input","temp_write_dir":"/root/coinstac/tmp/"}'

import glob, os,sys,json,argparse,shutil
import nibabel as nib

# Load Nipype interfaces
from nipype.interfaces import spm
from nipype.interfaces.spm.utils import DicomImport, ApplyTransform
from nipype.interfaces.spm import NewSegment, Smooth
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe

#Import corr.py (correlation) script
import corr


#Read and extract json args
args= json.loads(sys.argv[1])
input_bids_dir = args['input_bids_dir']
temp_write_dir = args['temp_write_dir']



#Check if input_bids_dir is in BIDS format using bids-validator tool and check if it has T1w data and write permissions to tmp write dir
cmd = "bids-validator {0}".format(input_bids_dir)
a=os.popen(cmd).read()

if (a) and ('T1w' in a) and (os.access(temp_write_dir, os.W_OK)):
        # Get the paths to the T1w files to run the algorithm
        smri_data = glob.glob(input_bids_dir + '/sub*/*/*T1w*.nii.gz')

        # Loop through each of the T1w*.nii.gz file to run the algorithm, this algorithm runs serially
        i=0
        while i<len(smri_data):
            gz=smri_data[i]
            i=i+1

            #Extract subject directory name from the T1w*.nii.gz files
            sub_id = gz.split('/')[-3]

            vbm_out = temp_write_dir + '/' + sub_id + '/anat'
            nii_output = (gz.split('/')[-1]).split('.gz')[0]

            # Create output dir for sub_id
            if not os.path.exists(vbm_out):
                    os.makedirs(vbm_out)

            nifti_file = vbm_out + '/' + nii_output

            # Connect spm12 standalone to nipype
            matlab_cmd = '/opt/spm12/run_spm12.sh /opt/mcr/v92 script'
            spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

            # Set the paths to the SPM12 Template and transform.mat
            transf_mat_path = '/root/vbm_scripts/transform.mat'
            tpm_path = '/root/vbm_scripts/TPM.nii'

            # Create vbm_spm12 dir under the specific sub-id/anat
            if not os.path.exists(vbm_out + "/vbm_spm12"):
                os.makedirs(vbm_out + "/vbm_spm12")

            # Create the VBM pipeline using Nipype
            #1 Reorientation node and settings
            reorient = pe.Node(interface=ApplyTransform(), name='reorient')
            reorient.inputs.mat = transf_mat_path
            if os.path.exists(nifti_file): reorient.inputs.in_file = nifti_file
            reorient.inputs.out_file = vbm_out + "/vbm_spm12/Re.nii"

            #2 Segementation Node and settings
            segmentation = pe.Node(interface=NewSegment(), name='segmentation')
            segmentation.inputs.channel_info = (0.0001, 60, (False, False))
            Tis1 = ((tpm_path, 1), 1, (True, False), (True, True))
            Tis2 = ((tpm_path, 2), 1, (True, False), (True, True))
            Tis3 = ((tpm_path, 3), 2, (True, False), (True, True))
            Tis4 = ((tpm_path, 4), 3, (True, False), (True, True))
            Tis5 = ((tpm_path, 5), 4, (True, False), (True, True))
            Tis6 = ((tpm_path, 6), 2, (True, False), (True, True))
            segmentation.inputs.tissues = [Tis1, Tis2, Tis3, Tis4, Tis5, Tis6]


            #3 Function & Node to transform the list of normalized class images to a compatible version for smoothing
            def transform_list(normalized_class_images):
                    return [each[0] for each in normalized_class_images]

            list_normalized_images = pe.Node(interface=Function(input_names='normalized_class_images', output_names='list_norm_images', function=transform_list),name='list_normalized_images')

            #4 Smoothing Node & Settings
            smoothing = pe.Node(interface=Smooth(), name='smoothing')
            smoothing.inputs.fwhm = [10, 10, 10]

            #5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir
            datasink = pe.Node(interface=DataSink(), name='sinker')
            datasink.inputs.base_directory = vbm_out

            #6 Create the pipeline/workflow and connect the nodes created above
            vbm_preprocess = pe.Workflow(name="vbm_preprocess")
            vbm_preprocess.connect([(reorient, segmentation, [('out_file', 'channel_files')]), \
                                                            (segmentation, list_normalized_images,
                                                             [('normalized_class_images', 'normalized_class_images')]), \
                                                            (list_normalized_images, smoothing, [('list_norm_images', 'in_files')]), \
                                                            (segmentation, datasink,
                                                             [('modulated_class_images', 'vbm_spm12'), ('native_class_images', 'vbm_spm12.@1'), \
                                                              ('normalized_class_images', 'vbm_spm12.@2'), ('transformation_mat', 'vbm_spm12.@3')]), \
                                                            (smoothing, datasink, [('smoothed_files', 'vbm_spm12.@4')])])

            try:
                    # Remove any tmp files in the docker
                    if (os.path.exists('/var/tmp')): shutil.rmtree('/var/tmp', ignore_errors=True)
                    for c in glob.glob(os.getcwd() + '/crash*'): os.remove(c)
                    for f in glob.glob(os.getcwd() + '/tmp*'): shutil.rmtree(f, ignore_errors=True)
                    for f in glob.glob(os.getcwd() + '/__pycache__'): shutil.rmtree(f, ignore_errors=True)
                    shutil.rmtree(os.getcwd() +'/vbm_preprocess', ignore_errors=True)
                    os.remove(os.getcwd() + '/pyscript.m')

                    # Pipeline execution starts here..

                    # gunzip *T1w*.gz files
                    n1_img = nib.load(gz)
                    nib.save(n1_img, vbm_out + '/' + nii_output)

                    # Run the VBM pipeline
                    if os.path.exists(nifti_file): res = vbm_preprocess.run()

                    #Calculate correlation coefficent of swc1*nii to SPM12 TPM.nii
                    segmented_file = glob.glob(vbm_out + "/vbm_spm12/swc1*nii")
                    corr_value = corr.get_corr(tpm_path, segmented_file[0])

                    #Write a text file that desribes what the output files are
                    with open(vbm_out + '/vbm_spm12/outputs_description.txt', 'w') as fp:
                        fp.write("Prefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air"
                                 "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
                                 "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf")
                    fp.close()

                    #Write a text file with info. on quality control correlation coefficent
                    with open(vbm_out + '/vbm_spm12/quality_control_readme.txt', 'w') as fp:
                        fp.write("vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox")
                    fp.close()

            except Exception as e:
                    # If fails raise the exception,set status False,write json output and print exception err string to std.err
                    status = False
                    sys.stderr.write(str(e))
                    continue

            else:
                    # If succeeds, set status True
                    status = True

            finally:
                    # Remove any tmp files in the docker
                    if (os.path.exists('/var/tmp')): shutil.rmtree('/var/tmp', ignore_errors=True)
                    for c in glob.glob(os.getcwd() + '/crash*'): os.remove(c)
                    for f in glob.glob(os.getcwd() + '/tmp*'): shutil.rmtree(f, ignore_errors=True)
                    for f in glob.glob(os.getcwd() + '/__pycache__'): shutil.rmtree(f, ignore_errors=True)
                    shutil.rmtree(os.getcwd() +'/vbm_preprocess', ignore_errors=True)
                    os.remove(os.getcwd() +'/pyscript.m')

                    # On the last subject in the BIDS directory , write the success status output to json object
                    if gz==smri_data[-1]:
                        sys.stdout.write(json.dumps({'success': status}, sort_keys=True, indent=4, separators=(',', ': ')))


# If input_bids_dir is not in BIDS format and does not have T1w data and no write permissions to tmp write dir then
# Set the Status to False, write the error message to stderr and output the json object on stdout
else:
    status = False
    sys.stderr.write("Make sure data is in BIDS format,T1w images exist and space is available on the system to write outputs")
    sys.stdout.write(json.dumps({'success': status}, sort_keys=True, indent=4, separators=(',', ': ')))
