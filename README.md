# Coinstac_VBM_computation
This computation runs Voxel Based Morphometry on structural T1 weighted MRI scans(BIDS format) using SPMv12 standalone and MATLAB Runtimev713. 

The computation runs serially and each scan takes approximately 10 mins to run on a system with 2.3 GHz,i5 equivalent processor and 8GB RAM. Each scan output directory takes about 150MB space after running this algorithm. Please make sure to have the space and resources to run this algorithm accordingly. Please read the contents of outputs_description.txt from the output directory to learn about the descriptions of output filed produced.

Test bids data is placed in /computation/test_bids_input_data

Example run below:

python3 /computation/run_vbm_bids.py --run '{"inputBidsDir":"/computation/test_bids_input_data","tempWriteDir":"/computation"}'

At the end of the computation , the output is written to stdout in json with the property "success" set to true/false indicating whether the algorithm ran successfully or failed
