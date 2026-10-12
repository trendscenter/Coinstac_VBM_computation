"""Container-side one-scan runner; never imported by production code."""

import json
from pathlib import Path
import shutil
import sys
import sysconfig


def validate_result(result, output_dir):
    """COINSTAC success alone does not mean the scan actually completed."""
    message = result.get("output", {}).get("message", "")
    if not result.get("success") or "1/1 subjects completed successfully." not in message or "Error log:" in message:
        raise RuntimeError("The single scan did not complete: " + message)
    output = output_dir / "vbm_outputs"
    subjects = [p for p in output.iterdir() if p.is_dir() and p.name != "covariates"]
    if len(subjects) != 1:
        raise RuntimeError("Expected exactly one preprocessed subject")
    anat = subjects[0] / "anat/vbm_spm12"
    required = ["Re.nii", "Re_seg8.mat", "wc1Re.png", "vbm_corr_value.txt"]
    required += ["%s%dRe.nii" % (prefix, tissue) for prefix in ("c", "wc", "mwc", "swc", "smwc") for tissue in range(1, 7)]
    for name in required:
        if not (anat / name).is_file() or (anat / name).stat().st_size == 0:
            raise RuntimeError("Missing/empty preprocessing artifact: " + name)
    if not (output_dir / "vbm_outputs.zip").is_file():
        raise RuntimeError("Missing downloadable output archive")


def main():
    shutil.copytree("/source", "/computation", dirs_exist_ok=True)
    # Use the exact interface override belonging to each source snapshot.
    shutil.copyfile(
        "/computation/preprocess.py",
        Path(sysconfig.get_paths()["purelib"]) / "nipype/interfaces/spm/preprocess.py",
    )
    # Fresh containers have no cached nodes. Set only the workflow work root:
    # SPM embeds it in Re_seg8.mat. Other temporary-file APIs remain untouched.
    import nipype.pipeline.engine as pe
    original_init = pe.Workflow.__init__

    def stable_workflow(self, name, base_dir=None, *args, **kwargs):
        original_init(self, name, base_dir or "/tmp/vbm-regression", *args, **kwargs)

    pe.Workflow.__init__ = stable_workflow
    sys.path.insert(0, "/computation")
    import run_vbm

    args = json.loads(Path("/regression/args.json").read_text())
    result = run_vbm.start(args)
    Path("/output/result.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    validate_result(result, Path("/output"))
    print(result["output"]["message"])


if __name__ == "__main__":
    main()
