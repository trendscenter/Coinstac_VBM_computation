import copy
import sys
import types
import unittest
from unittest import mock

# Keep fast-loop tests runnable outside the docker image.
if "nibabel" not in sys.modules:
    sys.modules["nibabel"] = types.SimpleNamespace(load=lambda _: None)
if "vbm.vbm_standalone_use_cases_layer" not in sys.modules:
    sys.modules["vbm.vbm_standalone_use_cases_layer"] = types.SimpleNamespace(
        setup_pipeline=lambda **_: {"success": True}
    )

from vbm import vbm_args, vbm_data
from vbm.vbm_config import template_dict as base_template_dict


class _FakeNifti:
    def __init__(self, shape):
        self.shape = shape


class RefactorSmokeTests(unittest.TestCase):
    def setUp(self):
        self.template = copy.deepcopy(base_template_dict)

    def test_args_parser_updates_selected_values(self):
        args = {
            "input": {
                "options_smoothing_x_mm": 8,
                "options_smoothing_y_mm": 9,
                "options_smoothing_z_mm": 10,
                "regression_resample_voxel_size": 4,
                "standalone": True,
            }
        }

        result = vbm_args.args_parser(args, self.template)

        self.assertIsNone(result)
        self.assertEqual(self.template["FWHM_SMOOTH"], [8.0, 9.0, 10.0])
        self.assertEqual(self.template["regression_resample_voxel_size"], (4.0, 4.0, 4.0))

    def test_args_parser_rejects_nonstandard_registration_template(self):
        args = {"input": {"registration_template": "/tmp/custom_tpm.nii"}}

        with mock.patch("vbm.vbm_args.os.path.isfile", return_value=True), mock.patch(
            "vbm.vbm_args.nib.load",
            side_effect=[_FakeNifti((121, 145, 121)), _FakeNifti((10, 10, 10))],
        ):
            result = vbm_args.args_parser(args, self.template)

        self.assertIsInstance(result, dict)
        self.assertEqual(result["success"], True)
        self.assertIn("Non-standard Registration template", result["output"]["message"])

    def test_data_parser_missing_file_raises(self):
        args = {
            "state": {"baseDirectory": "/input/local0/simulatorRun", "outputDirectory": "/tmp/out"},
            "input": {"covariates": ["missing_subject.nii.gz"]},
        }

        with mock.patch("vbm.vbm_data.os.path.isfile", return_value=False):
            with self.assertRaises(Exception) as context:
                vbm_data.data_parser(args, self.template)

        self.assertIn("File does not exist", str(context.exception))

    def test_data_parser_output_not_writable_raises(self):
        args = {
            "state": {"baseDirectory": "/input/local0/simulatorRun", "outputDirectory": "/tmp/out"},
            "input": {"covariates": ["subject.nii.gz"]},
        }

        with mock.patch("vbm.vbm_data.os.path.isfile", return_value=True), mock.patch(
            "vbm.vbm_data.os.access", return_value=False
        ):
            with self.assertRaises(Exception) as context:
                vbm_data.data_parser(args, self.template)

        self.assertIn("Output write permissions denied", str(context.exception))

    def test_args_parser_missing_input_raises_value_error(self):
        with self.assertRaises(ValueError) as context:
            vbm_args.args_parser({}, self.template)
        self.assertIn("missing 'input'", str(context.exception))

    def test_data_parser_missing_sections_raises_value_error(self):
        with self.assertRaises(ValueError) as context:
            vbm_data.data_parser({}, self.template)
        self.assertIn("missing 'state' or 'input'", str(context.exception))


if __name__ == "__main__":
    unittest.main()
