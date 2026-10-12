"""Remove an explicit set of NVIDIA compute libraries from a build-stage copy.

This is optional for the CPU SPM image, never applied to the reference runtime.
Keep CUDA's small core runtime and MATLAB's own libraries. In particular,
libicudata is an ICU Unicode library and must not be mistaken for CUDA.
"""

from pathlib import Path
import sys

GPU_LIBRARY_PREFIXES = (
    "libcudnn", "libcublas", "libcufft", "libcusolver", "libcusparse",
    "libcurand", "libnvrtc",
)


def prune(runtime_dir):
    directory = runtime_dir / "bin/glnxa64"
    if not directory.is_dir():
        raise ValueError("Not a MATLAB Linux runtime: " + str(runtime_dir))
    removed = 0
    for path in sorted(directory.iterdir()):
        if path.name.startswith(GPU_LIBRARY_PREFIXES) and (path.is_file() or path.is_symlink()):
            size = 0 if path.is_symlink() else path.stat().st_size
            path.unlink()
            removed += size
            print("Removed " + path.name)
    print("Removed %.2f GiB of GPU compute libraries" % (removed / 1024 ** 3))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: prune_runtime_gpu.py /path/to/MATLAB/runtime")
    prune(Path(sys.argv[1]))
