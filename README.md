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

## Output comparison helper
Use `compare_outputs.py` to compare two output archives during refactors.
It ignores macOS metadata files and compares SHA-256 hashes of decompressed
contents. Every file, including `.mat`, is strict by default. Any `--tolerate`
patterns are explicit exceptions and should not be used to prove compatibility.

Example:
`python3 compare_outputs.py test/output.zip test/output/local0/simulatorRun/vbm_outputs.zip --strip-prefix output/local0/simulatorRun/vbm_outputs --ignore "output/local0/simulatorRun/vbm_outputs.zip" --ignore "output/local0/simulatorRun/wc1Re.png" --tolerate "*.mat"`

## Full checkpoint regression
Use `scripts/full_regression.sh` to run the full workflow in one command:
clean output directory, build image, run simulator, and compare artifacts.

Example:
`scripts/full_regression.sh`

## One-scan byte compatibility regression

Build the MATLAB/SPM runtime once, then compare the committed implementation
against the current working tree without the simulator or another image build:

```sh
docker build --platform linux/amd64 -t coinstacteam/vbm_pre .
python3 scripts/single_file_regression.py --baseline-ref 208591d
```

The default scan is `test/input/local0/simulatorRun/A00063430.nii.gz`. Use
`--input /path/to/one-scan.nii.gz` to select another scan. Both implementations
run in fresh containers with the same pinned image, options, input, thread
settings, and internal paths. The test calls `run_vbm.start` and verifies actual
completion of exactly one subject, required artifacts, the downloadable archive,
the display PNG, and the returned JSON. Existing simulator outputs are preserved.

The first check preprocesses one file twice (baseline and candidate). Later
checks can reuse the baseline and preprocess just once:

```sh
python3 scripts/single_file_regression.py --baseline-ref 208591d \
  --baseline-run test/regression/<first-run-directory>
```

Every run retains source snapshots, logs, outputs, input/runtime provenance, and
`comparison.json` in an ignored `test/regression/` directory. Baseline reuse
requires the same input bytes, image ID, code revision, and options.

Comparison reads actual file bytes, rather than ZIP CRCs or floating-point
tolerances. ZIP timestamps/compression are packaging metadata; every decompressed
entry is compared. The only permitted file-byte normalization is MATLAB v5's
creation timestamp in its 116-byte descriptive header, reported in the comparison
JSON. The platform, remaining header, and entire MAT payload must match exactly.
Stable temporary paths in the test runner prevent random paths embedded in MAT
data from changing. Production preprocessing code is not patched for this.
`--raw-mat-header` also requires the creation timestamp to match (normally fails
across separate runs). Missing/added files and image, QC, covariate, PNG, response,
or MAT payload changes fail the test. Exit codes: 0 match, 1 mismatch, 2 run/setup
failure. A one-scan fixture covers this NIfTI path and default settings, not all
optional configurations or DICOM inputs.

Fast tests need only host Python:

```sh
python3 -m unittest discover -s test -p 'test_*.py'
```

## Smaller runtime image

`Dockerfile.slim` assembles a slim Python image from the installed runtime in the
original image. It copies the same MATLAB R2022b Update 9, SPM r7771, Python
packages, and plotting fonts, and leaves compilers, development headers, and
download tools behind. It explicitly copies production code instead of the
whole working directory. MATLAB libraries and SPM's CTF are retained until their
removal can be validated; the MATLAB Runtime remains a substantial size floor.

```sh
docker build --platform linux/amd64 -f Dockerfile.slim -t coinstacteam/vbm_pre:slim .
python3 scripts/single_file_regression.py --baseline-ref 208591d \
  --baseline-run test/regression/<first-run-directory> \
  --candidate-image coinstacteam/vbm_pre:slim
```

The original image is needed only at build time; the final slim image includes
its runtime. This reduces the distributable image, though retaining both images
and build caches locally still consumes disk space.

For the CPU-only variant, `--build-arg STRIP_GPU_LIBS=1` additionally removes an
explicit set of NVIDIA cuDNN/cuBLAS/cuFFT/cuSOLVER/cuSPARSE/cuRAND/NVRTC libraries
in the source stage. The reference runtime and default slim image keep these
libraries. Run the same regression against this variant before using it:

```sh
docker build --platform linux/amd64 -f Dockerfile.slim \
  --build-arg STRIP_GPU_LIBS=1 -t coinstacteam/vbm_pre:cpu .
python3 scripts/single_file_regression.py --baseline-ref 208591d \
  --baseline-run test/regression/<first-run-directory> \
  --candidate-image coinstacteam/vbm_pre:cpu
```

This variant retains CUDA's core runtime, MATLAB's own libraries, SPM's CTF,
and the same software versions. It is intended for this CPU SPM computation;
GPU-dependent MATLAB applications need the full runtime.

For the smallest verified variant, select the runtime products declared by SPM's
existing `requiredMCRProducts.txt`:

```sh
docker build --platform linux/amd64 -f Dockerfile.slim \
  --build-arg SELECT_SPM_PRODUCTS=1 -t coinstacteam/vbm_pre:bespoke .
python3 scripts/single_file_regression.py --baseline-ref 208591d \
  --baseline-run test/regression/<first-run-directory> \
  --candidate-image coinstacteam/vbm_pre:bespoke
```

`scripts/select_runtime_products.py` uses the installed MathWorks product,
component, and file manifests to retain SPM's Core, Graphics, Non Interactive
MATLAB, and Numerics products plus the runtime foundation. The selected runtime
files total 2.15 GiB. Binary files remain unchanged, licenses and notices are
retained, and installation indexes describe the selected products. The build
also checks that standalone SPM starts. The selection report is stored at
`/opt/mcr/runtime-selection.json` inside the image.

This assembles a subset of the existing R2022b Update 9 installation; it does
not recompile SPM or generate a MathWorks custom installer. No MATLAB Compiler
installation is needed. Other MATLAB applications may require additional
products. The original image and Docker build caches still consume local space.

Verified locally using the single A00063430 scan (bespoke check: 2026-10-12):

| Image | Docker image size | One-scan result |
| --- | --- | --- |
| `coinstacteam/vbm_pre:latest` | 11 GB | Refactor matches committed baseline |
| `coinstacteam/vbm_pre:slim` | 9.94 GB | Built; not separately scanned |
| `coinstacteam/vbm_pre:cpu` | 7.38 GB | Matches full-runtime baseline |
| `coinstacteam/vbm_pre:bespoke` | 3.50 GB | Matches full-runtime baseline |

All three verified comparisons matched all 100 archive entries after only MAT
creation-header timestamp normalization, plus the display PNG and returned JSON
exactly. The bespoke image reduces the distributable size by about 68% compared
with the original, and 53% compared with the CPU variant. Its comparison report
is `test/regression/20261012T055340Z-5d76zw15/comparison.json`.
The 30 fast tests also pass.
These results cover one NIfTI scan with default preprocessing options.
