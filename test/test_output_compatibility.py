import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from zipfile import ZipFile, ZIP_DEFLATED, ZIP_STORED

import compare_outputs
from scripts.prune_runtime_gpu import prune

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("single_file_regression", ROOT / "scripts/single_file_regression.py")
regression = importlib.util.module_from_spec(spec)
spec.loader.exec_module(regression)
runner_spec = importlib.util.spec_from_file_location("regression_container", ROOT / "scripts/regression_container.py")
runner = importlib.util.module_from_spec(runner_spec)
runner_spec.loader.exec_module(runner)


class OutputCompatibilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def archive(self, name, files, compression=ZIP_STORED):
        path = self.root / name
        with ZipFile(path, "w", compression=compression) as archive:
            for filename, content in files.items():
                archive.writestr(filename, content)
        return path

    def compare(self, left, right, **kwargs):
        return regression.compare_archives(self.archive("a.zip", left), self.archive("b.zip", right), **kwargs)

    @staticmethod
    def mat(timestamp, payload=b"scientific-data"):
        return ("MATLAB 5.0 MAT-file, Platform: GLNXA64, Created on: " + timestamp).encode().ljust(116, b" ") + b"\0" * 8 + b"\0\1IM" + payload

    def test_zip_container_metadata_and_compression_do_not_affect_content(self):
        a = self.archive("a.zip", {"scan.nii": b"same bytes"}, ZIP_STORED)
        b = self.archive("b.zip", {"scan.nii": b"same bytes"}, ZIP_DEFLATED)
        self.assertNotEqual(a.read_bytes(), b.read_bytes())
        self.assertEqual(regression.compare_archives(a, b)["differences"], [])

    def test_one_changed_image_byte_fails(self):
        report = self.compare({"scan.nii": b"1234"}, {"scan.nii": b"1235"})
        self.assertEqual(report["differences"], ["scan.nii (bytes changed)"])

    def test_generic_archive_comparator_does_not_tolerate_mat_by_default(self):
        with mock.patch("sys.argv", ["compare_outputs.py", "a.zip", "b.zip"]):
            self.assertEqual(compare_outputs.parse_args().tolerate, [])
        a = self.archive("a.zip", {"seg.mat": b"original"})
        b = self.archive("b.zip", {"seg.mat": b"modified"})
        result = compare_outputs.compare_manifests(
            compare_outputs.build_manifest(a, [], ()),
            compare_outputs.build_manifest(b, [], ()), (),
        )
        self.assertEqual(result[2], ["seg.mat"])
        self.assertEqual(result[3], [])

    def test_added_and_missing_files_fail(self):
        report = self.compare({"old.nii": b"x"}, {"new.nii": b"x"})
        self.assertEqual(len(report["differences"]), 2)

    def test_empty_archives_cannot_pass(self):
        with self.assertRaisesRegex(ValueError, "Empty"):
            self.compare({}, {})

    def test_timestamp_only_is_reported_and_normalized(self):
        report = self.compare({"seg.mat": self.mat("Monday")}, {"seg.mat": self.mat("Tuesday")})
        self.assertEqual(report["differences"], [])
        self.assertEqual(report["normalized_mat_headers"][0]["path"], "seg.mat")

    def test_raw_timestamp_mode_fails(self):
        report = self.compare({"seg.mat": self.mat("Monday")}, {"seg.mat": self.mat("Tuesday")}, normalize_mat_header=False)
        self.assertEqual(len(report["differences"]), 1)

    def test_mat_data_changes_cannot_be_hidden_by_normalization(self):
        report = self.compare({"seg.mat": self.mat("Monday", b"a")}, {"seg.mat": self.mat("Tuesday", b"b")})
        self.assertEqual(len(report["differences"]), 1)
        self.assertEqual(report["normalized_mat_headers"], [])

    def test_embedded_temporary_path_changes_fail(self):
        report = self.compare({"seg.mat": self.mat("Monday", b"/tmp/one")}, {"seg.mat": self.mat("Tuesday", b"/tmp/two")})
        self.assertEqual(len(report["differences"]), 1)

    def test_mat_version_and_platform_changes_fail(self):
        a = self.mat("Monday")
        for b in (a[:126] + b"XX" + a[128:], a.replace(b"GLNXA64", b"WINDOWS")):
            report = self.compare({"seg.mat": a}, {"seg.mat": b})
            self.assertEqual(len(report["differences"]), 1)

    def test_unknown_mat_formats_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unrecognized"):
            self.compare({"seg.mat": b"a"}, {"seg.mat": b"b"})

    def test_duplicate_zip_entries_are_rejected(self):
        a = self.archive("a.zip", {"scan.nii": b"a"})
        with ZipFile(a, "a") as archive:
            with self.assertWarns(UserWarning):
                archive.writestr("scan.nii", b"b")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            regression.compare_archives(a, a)

    def test_pipeline_success_true_cannot_hide_a_failed_scan(self):
        with self.assertRaisesRegex(RuntimeError, "did not complete"):
            runner.validate_result({"success": True, "output": {"message": " Error log: segmentation failed"}}, self.root)

    def test_partial_completion_and_error_messages_are_rejected(self):
        for message in ("0/1 subjects completed successfully.", "1/1 subjects completed successfully. Error log: failed"):
            with self.assertRaisesRegex(RuntimeError, "did not complete"):
                runner.validate_result({"success": True, "output": {"message": message}}, self.root)

    def test_missing_preprocessed_image_is_rejected(self):
        (self.root / "vbm_outputs/subject/anat/vbm_spm12").mkdir(parents=True)
        with self.assertRaisesRegex(RuntimeError, "Missing/empty"):
            runner.validate_result({"success": True, "output": {"message": "1/1 subjects completed successfully."}}, self.root)

    def test_multiple_subjects_are_rejected(self):
        (self.root / "vbm_outputs/subject1").mkdir(parents=True)
        (self.root / "vbm_outputs/subject2").mkdir()
        with self.assertRaisesRegex(RuntimeError, "exactly one"):
            runner.validate_result({"success": True, "output": {"message": "1/1 subjects completed successfully."}}, self.root)

    def test_cpu_image_pruning_preserves_unicode_and_core_runtime_libraries(self):
        directory = self.root / "bin/glnxa64"
        directory.mkdir(parents=True)
        retained = ("libicudata.so.69", "libcudart.so.11", "libmwgpuruntime.so")
        for name in (*retained, "libcublas.so.11.3.1.68"):
            (directory / name).write_bytes(b"library")
        (directory / "libcublas.so.11").symlink_to("libcublas.so.11.3.1.68")
        with mock.patch("builtins.print"):
            prune(self.root)
        self.assertEqual({p.name for p in directory.iterdir()}, set(retained))


if __name__ == "__main__":
    unittest.main()
