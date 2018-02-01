"""************************** run_vbm_bids.py ****************"""
""" This script runs vbm pipeline on BIDS anatomical data using spm12 standalone and Matlab common runtime"""
"""Example run of the code- python3 /computation/run_vbm_bids.py --run '{"inputBidsDir":"/computation/test_bids_input_data","SmoothingValue":"[6, 6, 6]","tempWriteDir":"/computation/bids_vbm_output"}'"""
""" Input args: --run json (this json structure may involve different field for different run) """
""" output: json """

## script name: run_vbm_bids.py ##
## import dependent libraries ##
import glob, os, sys, json, argparse, shutil, ast
import nibabel as nib

## Load Nipype interfaces ##
from nipype.interfaces import spm
from nipype.interfaces.spm.utils import DicomImport, ApplyTransform
from nipype.interfaces.spm import NewSegment, Smooth
from nipype.interfaces.io import DataSink
from nipype.interfaces.utility import Function
import nipype.pipeline.engine as pe

## Import corr.py (correlation) script ##
import corr

SCAN_TYPE = 'T1w'
smooth_mm_value = [10, 10, 10]
# Set the paths to the SPM12 Template and transform.mat
transf_mat_path = '/computation/transform.mat'
spm_path=["/opt/spm12/fsroot"]
tpm_path = '/opt/spm12/fsroot/spm/spm12/tpm/TPM.nii'

if __name__=='__main__':
    # Read and extract json args
    parser = argparse.ArgumentParser(description='read in coinstac args for local computation')
    parser.add_argument('--run', type=str, help='grab coinstac args')
    args = parser.parse_args()
    args.run = json.loads(args.run)
    input_bids_dir = args.run['inputBidsDir']
    temp_write_dir = args.run['tempWriteDir']
    smooth_mm_value = ast.literal_eval(args.run['SmoothingValue'])

    # Check if input_bids_dir is in BIDS format using bids-validator tool and check if it has T1w data and write permissions to tmp write dir
    cmd = "bids-validator {0}".format(input_bids_dir)
    bids_process = os.popen(cmd).read()

    if bids_process and (SCAN_TYPE in bids_process) and os.access(temp_write_dir, os.W_OK):
        # Get the paths to the T1w files to run the algorithm
        glob_str = os.path.join(input_bids_dir, 'sub*', '*', '*T1w*.nii.gz')
        smri_data = glob.glob(glob_str)

        ## Loop through each of the T1w*.nii.gz file to run the algorithm, this algorithm runs serially ##
        i = 0  # variable for looping
        count_success = 0  # variable for counting how many subjects were successfully run
        
        #create dirs array to store output directories where vbm spm12 data is written to
        dirs=[]
        
        #create wc1files array to store paths to wc1 files for each subject
        wc1files=[]
        
        while i < len(smri_data):
            gzip_file_path = smri_data[i]
            i = i + 1

            # Extract subject directory name from the T1w*.nii.gz files
            sub_path = (os.path.dirname(os.path.dirname(gzip_file_path))).split('/')
            sub_id='/'.join(sub_path[2:len(sub_path)])

            vbm_out = temp_write_dir + '/' + sub_id + '/anat'
            nii_output = (gzip_file_path.split('/')[-1]).split('.gz')[0]

            # Create output dir for sub_id
            os.makedirs(vbm_out, exist_ok=True)

            nifti_file = os.path.join(vbm_out, nii_output)
            n1_img = nib.load(gzip_file_path)
            nib.save(n1_img, os.path.join(vbm_out, nii_output))
            
            

            # Connect spm12 standalone to nipype
            matlab_cmd = '/opt/spm12/run_spm12.sh /opt/mcr/v92 script'
            spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)

            # Create vbm_spm12 dir under the specific sub-id/anat
            os.makedirs(os.path.join(vbm_out + "/vbm_spm12"), exist_ok=True)

            # Create the VBM pipeline using Nipype
            ## 1 Reorientation node and settings ##
            reorient = pe.Node(interface=ApplyTransform(), name='reorient')
            reorient.inputs.mat = transf_mat_path
            if os.path.exists(nifti_file): reorient.inputs.in_file = nifti_file
            reorient.inputs.out_file = vbm_out + "/vbm_spm12/Re.nii"
            reorient.inputs.paths=spm_path

            ## 2 Segementation Node and settings ##
            segmentation = pe.Node(interface=NewSegment(), name='segmentation')
            segmentation.inputs.channel_info = (0.0001, 60, (False, False))
            Tis1 = ((tpm_path, 1), 1, (True, False), (True, True))
            Tis2 = ((tpm_path, 2), 1, (True, False), (True, True))
            Tis3 = ((tpm_path, 3), 2, (True, False), (True, True))
            Tis4 = ((tpm_path, 4), 3, (True, False), (True, True))
            Tis5 = ((tpm_path, 5), 4, (True, False), (True, True))
            Tis6 = ((tpm_path, 6), 2, (True, False), (True, True))
            segmentation.inputs.tissues = [Tis1, Tis2, Tis3, Tis4, Tis5, Tis6]
            segmentation.inputs.paths=spm_path


            ## 3 Function & Node to transform the list of normalized class images to a compatible version for smoothing ##
            def transform_list(normalized_class_images):
                return [each[0] for each in normalized_class_images]


            interface = Function(
                input_names='normalized_class_images',
                output_names='list_norm_images',
                function=transform_list
            )
            list_normalized_images = pe.Node(interface=interface, name='list_normalized_images')

            ## 4 Smoothing Node & Settings ##
            smoothing = pe.Node(interface=Smooth(), name='smoothing')
            smoothing.inputs.fwhm = smooth_mm_value
            smoothing.inputs.paths=spm_path
            

            ## 5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir ##
            datasink = pe.Node(interface=DataSink(), name='sinker')
            datasink.inputs.base_directory = vbm_out

            ## 6 Create the pipeline/workflow and connect the nodes created above ##
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
                ## Remove any tmp files in the docker ##
                if (os.path.exists('/var/tmp')):
                    shutil.rmtree('/var/tmp', ignore_errors=True)

                for c in glob.glob(os.getcwd() + '/crash*'):
                    os.remove(c)

                for f in glob.glob(os.getcwd() + '/tmp*'):
                    shutil.rmtree(f, ignore_errors=True)

                for f in glob.glob(os.getcwd() + '/__pycache__'):
                    shutil.rmtree(f, ignore_errors=True)

                if os.path.exists(os.getcwd() + '/vbm_preprocess'):
                    shutil.rmtree(os.getcwd() + '/vbm_preprocess', ignore_errors=True)

                if os.path.exists(os.getcwd() + '/pyscript.m'):
                    os.remove(os.getcwd() + '/pyscript.m')

                ## Pipeline execution starts here.. ##

                # Run the VBM pipeline
                if os.path.exists(nifti_file):
                    res = vbm_preprocess.run()
                    dirs.append(os.path.join(vbm_out + "/vbm_spm12"))
                    wc1files.append(glob.glob(vbm_out + "/vbm_spm12/wc1*nii")[0])

                # Calculate correlation coefficent of swc1*nii to SPM12 TPM.nii
                segmented_file = glob.glob(vbm_out + "/vbm_spm12/swc1*nii")
                corr_value = corr.get_corr(tpm_path, segmented_file[0])

                ## Write a text file that desribes what the output files are ##
                with open(vbm_out + '/vbm_spm12/outputs_description.txt', 'w') as fp:
                    fp.write(
                        "Prefixes descriptions for segmented images:c1-Grey matter,c2-White matter,c3-Cerebro spinal fluid,c4-Bone,c5-Soft tissue,c6-Air"
                        "\nw-Normalized\nm-Modulated\ns-Smoothed with fwhm(mm) [10 10 10]\nFor more info. please refer to spm12 manual here: "
                        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf")
                fp.close()

                # Write a text file with info. on quality control correlation coefficent
                with open(vbm_out + '/vbm_spm12/quality_control_readme.txt', 'w') as fp:
                    fp.write(
                        "vbm_corr_value.txt gives the correlation value of the swc1*nii file with spm12/tpm/TPM.nii file from SPM12 toolbox")
                fp.close()

            except Exception as e:
                # If fails raise the exception,set status False,write json output and print exception err string to std.err
                status = True
                sys.stderr.write(str(e))
                continue

            else:
                # If the try block succeeds, increase the count
                count_success = count_success + 1

            finally:
                ## Remove any tmp files in the docker ##
                if (os.path.exists('/var/tmp')): shutil.rmtree('/var/tmp', ignore_errors=True)
                for c in glob.glob(os.getcwd() + '/crash*'): os.remove(c)
                for f in glob.glob(os.getcwd() + '/tmp*'): shutil.rmtree(f, ignore_errors=True)
                for f in glob.glob(os.getcwd() + '/__pycache__'): shutil.rmtree(f, ignore_errors=True)
                if os.path.exists(os.getcwd() + '/vbm_preprocess'): shutil.rmtree(os.getcwd() + '/vbm_preprocess',
                                                                                  ignore_errors=True)
                if os.path.exists(os.getcwd() + '/pyscript.m'): os.remove(os.getcwd() + '/pyscript.m')

                # On the last subject in the BIDS directory , write the success status output to json object
                if gzip_file_path == smri_data[-1]:
                    if count_success > 0: status = True  # If atleast 1 scan in the BIDS directory finishes successfully
                    sys.stdout.write(json.dumps({"output": {"success": status,"vbmdirs":dirs,"wc1files":wc1files}}, indent=4, separators=(',', ': ')))


    # If input_bids_dir is not in BIDS format and does not have T1w data and no write permissions to tmp write dir then
    # Set the Status to False, write the error message to stderr and output the json object on stdout
    else:
        status = True
        sys.stderr.write(
            "Make sure data is in BIDS format,T1w images exist and space is available on the system to write outputs")
        sys.stdout.write(json.dumps({"output": {"success": status}}, indent=4, separators=(',', ': ')))
