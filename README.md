The foll. should work for testing instead of building docker from scratch
docker pull spanta28/coinstac_vbm_docker 

# Coinstac_VBM_computation
This computation runs Voxel Based Morphometry on structural T1 weighted MRI scans(BIDS format) using SPMv12 standalone and MATLAB Runtimev713. 

The computation runs serially and each scan takes approximately 10 mins to run on a system with 2.3 GHz,i5 equivalent processor and 8GB RAM. Each scan output directory takes about 150MB space after running this algorithm. Please make sure to have the space and resources to run this algorithm accordingly. Please read the contents of outputs_description.txt from the output directory to learn about the descriptions of output filed produced.

Test bids data is placed in /computation/test_dir/bids_input_data

At the end of the computation , the output is written to stdout in json with the property "success" set to true/false indicating whether the algorithm ran successfully or failed

Dockerfile_build_from_scratch builds everything from scratch, downloads packages from online
Dockerfile sources from trendscenter/coinstac_nipype_spm_mcr_afni_centos7 and just uploads the current scripts into the docker.

For vbm+regression, use test/inputspec_vbm_regression.json as inputspec.json

if you use sample data in test directory, example output below:
{
  "data": [
    [
      "regression_input_files/subID-1_swc1_4mm.nii",
      "regression_input_files/subID-2_swc1_4mm.nii"
    ],
    [
      "niftifile"
    ]
  ],
  "covariates": [
    [
      [
        [
          "niftifile",
          "isControl",
          "age",
          "sex"
        ],
        [
          "regression_input_files/subID-1_swc1_4mm.nii",
          true,
          28,
          "M"
        ],
        [
          "regression_input_files/subID-2_swc1_4mm.nii",
          true,
          39,
          "M"
        ]
      ]
    ],
    [
      "isControl",
      "age",
      "sex"
    ],
    [
      "boolean",
      "number",
      "string"
    ]
  ]
}

This is how outputDirectory looks:
cd test/output/
ls */*/*
local0/simulatorRun/vbm_outputs.zip	local0/simulatorRun/wc1Re.png

local0/simulatorRun/regression_input_files:
subID-1_swc1_4mm.nii	subID-2_swc1_4mm.nii

local0/simulatorRun/vbm_outputs:
outputs_description.txt		quality_control_readme.txt	subID-1				subID-2				vbm_log.txt

vbm+regression optional params : 
options_smoothing,reorientation params, regression_resampling_voxel_size etc. regression_resampling_voxel_size can be varied 
to check the performance of coinstac in terms of how much data flow it can handle. Default value for this is 4mm.
For more info. check compspec.json
Use inputspec_vbm_regression_all_options.json as inputspec.json


For standalone vbm, use test/inputspec_vbm_standalone.json as inputspec.json

example outputDirectory below:
cd test/output/
ls */*/*
local0/simulatorRun/vbm_outputs.zip	local0/simulatorRun/wc1Re.png

local0/simulatorRun/vbm_outputs:
outputs_description.txt		quality_control_readme.txt	subID-1				subID-2				vbm_log.txt

For standalone vbm options, use test/inputspec_vbm_standalone_all_options.json as inputspec.json

The optional params info is detailed in vbm optional params description.xlsx
Type of info. available is : Name	Label	Description	Json data type	Algorithm dev type for inner processing	Default value	Range (if applicable)	Step increments	pre=processing step
