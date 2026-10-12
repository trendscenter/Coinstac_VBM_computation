#!/usr/bin/env python3
"""Run one scan per implementation in an identical Docker runtime and compare bytes.

No simulator, rebuild, host dependency installs, or existing-output deletion.
Runs and provenance are retained under test/regression/ for inspection/reuse.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
CHUNK = 1024 * 1024


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def snapshot(destination, ref=None):
    if ref:
        names = git("ls-tree", "-r", "--name-only", "-z", ref).decode().split("\0")
    else:
        names = git("ls-files", "-c", "-o", "--exclude-standard", "-z").decode().split("\0")
    for name in sorted(set(names)):
        # The runtime needs code, not the full multi-subject data set.
        if not name or name.startswith(("test/", "test_dir/", "scripts/")):
            continue
        if Path(name).suffix not in (".py", ".m", ".mat"):
            continue
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if ref:
            target.write_bytes(git("show", ref + ":" + name))
        else:
            shutil.copyfile(ROOT / name, target)


def bytes_equal(left, right):
    while True:
        a, b = left.read(CHUNK), right.read(CHUNK)
        if a != b:
            return False
        if not a:
            return True


def mat_header(data):
    # MATLAB v5 has a 116-byte descriptive field; the remaining header includes
    # format/version/endian data and must remain exact. Never normalize payloads.
    if len(data) < 128 or not data.startswith(b"MATLAB 5.0 MAT-file") or b"Created on:" not in data[:116]:
        raise ValueError("Unrecognized MATLAB v5 header")
    return data[:116]


def compare_archives(baseline, candidate, normalize_mat_header=True):
    differences, normalized = [], []
    with ZipFile(baseline) as a, ZipFile(candidate) as b:
        def entries(archive):
            files = {}
            for info in archive.infolist():
                if info.is_dir():
                    continue
                if info.filename in files:
                    raise ValueError("Duplicate archive entry: " + info.filename)
                files[info.filename] = info
            if not files:
                raise ValueError("Empty preprocessing archive")
            return files

        aa, bb = entries(a), entries(b)
        for name in sorted(aa.keys() | bb.keys()):
            if name not in aa or name not in bb:
                differences.append(name + " (missing/added)")
                continue
            with a.open(aa[name]) as x, b.open(bb[name]) as y:
                if bytes_equal(x, y):
                    continue
            if normalize_mat_header and name.endswith(".mat"):
                x, y = a.read(aa[name]), b.read(bb[name])
                hx, hy = mat_header(x), mat_header(y)
                # Preserve platform and all header text preceding the timestamp.
                if hx.split(b"Created on:")[0] == hy.split(b"Created on:")[0] and x[116:] == y[116:]:
                    normalized.append({"path": name, "baseline_header": hx.decode("ascii"),
                                       "candidate_header": hy.decode("ascii")})
                    continue
            differences.append(name + " (bytes changed)")
    return {"files": len(aa), "differences": differences, "normalized_mat_headers": normalized}


def run(image, source, output, fixture, harness, label):
    output.mkdir()
    command = ["docker", "run", "--rm", "--platform", "linux/amd64", "--network", "none",
               "--mount", "type=tmpfs,target=/computation",
               "--workdir", "/computation", "-e", "PYTHONHASHSEED=0", "-e", "MPLBACKEND=Agg",
               "-e", "OMP_NUM_THREADS=1", "-e", "OPENBLAS_NUM_THREADS=1", "-e", "MKL_NUM_THREADS=1"]
    for host, container, readonly in (
        (source, "/source", True), (output, "/output", False),
        (fixture, "/input", True), (harness, "/regression", True),
    ):
        command += ["--mount", "type=bind,source=%s,target=%s%s" % (
            host, container, ",readonly" if readonly else "")]
    command += ["--entrypoint", "python", image, "/regression/regression_container.py"]
    print("Preprocessing one scan (%s); log: %s" % (label, output / "run.log"), flush=True)
    with (output / "run.log").open("w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("%s failed; inspect %s" % (label, output / "run.log"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", default="coinstacteam/vbm_pre")
    parser.add_argument("--candidate-image", help="Compare another runtime image as well (e.g. the slim image)")
    parser.add_argument("--baseline-ref", default="HEAD", help="Known-good code revision (default HEAD)")
    parser.add_argument("--baseline-run", type=Path, help="Reuse a prior run's baseline; preprocess only the candidate")
    parser.add_argument("--input", type=Path, default=ROOT / "test/input/local0/simulatorRun/A00063430.nii.gz")
    parser.add_argument("--raw-mat-header", action="store_true", help="Require timestamp bytes to match as well")
    args = parser.parse_args()
    fixture = args.input.resolve(strict=True)
    if not fixture.name.endswith((".nii", ".nii.gz")):
        parser.error("--input must be a single .nii or .nii.gz file")
    # Pin both runs to the same immutable runtime, even if a tag changes midway.
    image = subprocess.check_output(["docker", "image", "inspect", args.image,
                                     "--format", "{{.Id}}"], text=True).strip()
    candidate_image = subprocess.check_output(["docker", "image", "inspect", args.candidate_image,
                                               "--format", "{{.Id}}"], text=True).strip() if args.candidate_image else image
    ref = git("rev-parse", args.baseline_ref + "^{commit}").decode().strip()
    spec = {
        "state": {"baseDirectory": "/input", "outputDirectory": "/output"},
        "input": {"options_smoothing_x_mm": 10.0, "options_smoothing_y_mm": 10.0,
                  "options_smoothing_z_mm": 10.0, "standalone": False,
                  "covariates": {fixture.name: {"isControl": "True", "age": "28", "sex": "M", "site": "MGM"}}},
    }
    provenance = {"image": image, "baseline_ref": ref,
                  "input_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(), "args": spec,
                  "runner_sha256": hashlib.sha256(Path(__file__).with_name("regression_container.py").read_bytes()).hexdigest()}
    previous = args.baseline_run.resolve() if args.baseline_run else None
    if previous:
        old = json.loads((previous / "provenance.json").read_text())
        if old != provenance:
            raise ValueError("Baseline runtime, revision, input, or options differ; generate a fresh baseline")
    runs = ROOT / "test/regression"
    runs.mkdir(parents=True, exist_ok=True)
    # mkdtemp also avoids output reuse, races, and stale successful artifacts.
    import tempfile
    run_dir = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ-"), dir=runs))
    print("Regression artifacts: " + str(run_dir), flush=True)
    (run_dir / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    harness = run_dir / "harness"
    harness.mkdir()
    shutil.copyfile(Path(__file__).with_name("regression_container.py"), harness / "regression_container.py")
    (harness / "args.json").write_text(json.dumps(spec))
    inputs = run_dir / "input"
    inputs.mkdir()
    shutil.copyfile(fixture, inputs / fixture.name)
    baseline = previous / "baseline" if previous else run_dir / "baseline"
    if not previous:
        snapshot(run_dir / "baseline-source", ref)
        run(image, run_dir / "baseline-source", baseline, inputs, harness, "baseline")
    snapshot(run_dir / "candidate-source")
    candidate = run_dir / "candidate"
    run(candidate_image, run_dir / "candidate-source", candidate, inputs, harness, "candidate")
    report = compare_archives(baseline / "vbm_outputs.zip", candidate / "vbm_outputs.zip", not args.raw_mat_header)
    report["candidate_image"] = candidate_image
    for name in ("wc1Re.png", "result.json"):
        with (baseline / name).open("rb") as x, (candidate / name).open("rb") as y:
            if not bytes_equal(x, y):
                report["differences"].append(name + " (bytes changed)")
    (run_dir / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print("FAIL" if report["differences"] else "PASS")
    return bool(report["differences"])


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(str(error), file=sys.stderr)
        raise SystemExit(2)
