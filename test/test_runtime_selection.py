import json
from pathlib import Path
import tempfile
import unittest

from scripts.select_runtime_products import assemble, filter_installation_registry, plan_selection


class RuntimeSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "runtime"
        for name in ("products", "components", "files"):
            (self.source / "appdata" / name).mkdir(parents=True)
        self.component = "engine_glnxa64 42"
        self.product = "appdata/products/Core.xml"
        (self.source / self.product).write_text(
            "<productData><productNumber>35000</productNumber>"
            "<productName>Core</productName><releaseDescription>R2022b_Update_9</releaseDescription>"
            "<dependsOn><component><name>engine_glnxa64</name><version>42</version>"
            "</component></dependsOn></productData>"
        )
        (self.source / "appdata/components" / (self.component + ".xml")).write_text("<componentData/>")
        self.manifest = self.source / "appdata/files" / (self.component + "_manifest.json")
        self.manifest.write_text(json.dumps(["bin", "bin/glnxa64", "bin/glnxa64/alias.so"]))
        directory = self.source / "bin/glnxa64"
        directory.mkdir(parents=True)
        (directory / "needed.so").write_bytes(b"original binary bytes\x00\xff")
        (directory / "alias.so").symlink_to("needed.so")
        (directory / "unrelated.so").write_bytes(b"unused")
        (self.source / "VersionInfo.xml").write_text("R2022b Update 9")
        (self.source / "matlabruntime_license_agreement.pdf").write_bytes(b"license")
        indexes = {
            "installedProductData.json": {"1000": 1, "35000": 1, "35011": 1},
            "prodcontents.json": {"Core": ["/" + self.product], "Other": ["/appdata/products/Other.xml"]},
            "compcontents.json": {self.component: ["/appdata/components/" + self.component + ".xml"], "other": []},
        }
        for name, value in indexes.items():
            (self.source / "appdata" / name).write_text(json.dumps(value))
        self.requirements = self.root / "requiredMCRProducts.txt"
        self.requirements.write_text("35000\t\n")

    def test_subset_preserves_declared_bytes_symlink_targets_and_licenses(self):
        destination = self.root / "selected/R2022b"
        report = assemble(self.source, self.requirements, destination)
        self.assertEqual(report["products"], {"35000": "Core"})
        self.assertEqual((destination / "bin/glnxa64/alias.so").read_bytes(), b"original binary bytes\x00\xff")
        self.assertTrue((destination / "bin/glnxa64/alias.so").is_symlink())
        self.assertFalse((destination / "bin/glnxa64/unrelated.so").exists())
        self.assertEqual((destination / "matlabruntime_license_agreement.pdf").read_bytes(), b"license")
        self.assertEqual(json.loads((destination / "appdata/installedProductData.json").read_text()), {"1000": 1, "35000": 1})
        self.assertEqual(set(json.loads((destination / "appdata/prodcontents.json").read_text())), {"Core"})
        self.assertEqual(set(json.loads((destination / "appdata/compcontents.json").read_text())), {self.component})

    def test_missing_declared_file_fails_before_creating_destination(self):
        self.manifest.write_text(json.dumps(["bin/missing.so"]))
        destination = self.root / "selected"
        with self.assertRaisesRegex(ValueError, "Missing declared"):
            assemble(self.source, self.requirements, destination)
        self.assertFalse(destination.exists())

    def test_unknown_required_product_fails(self):
        self.requirements.write_text("99999")
        with self.assertRaisesRegex(ValueError, "Missing runtime products"):
            plan_selection(self.source, self.requirements)

    def test_manifest_cannot_traverse_outside_runtime(self):
        self.manifest.write_text(json.dumps(["../outside"]))
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            plan_selection(self.source, self.requirements)

    def test_external_symlink_targets_are_rejected(self):
        outside = self.root / "outside.so"
        outside.write_bytes(b"outside")
        link = self.source / "bin/glnxa64/alias.so"
        link.unlink()
        link.symlink_to(outside)
        with self.assertRaisesRegex(ValueError, "outside installation"):
            plan_selection(self.source, self.requirements)

    def test_existing_destination_cannot_be_overwritten(self):
        destination = self.root / "selected"
        destination.mkdir()
        with self.assertRaisesRegex(ValueError, "must be new"):
            assemble(self.source, self.requirements, destination)

    def test_installation_registry_retains_vendor_numeric_literals(self):
        text = '{"1000":1791598517716,"35000":1.79159851766E+12,"35011":1.791598517665E+12}'
        self.assertEqual(filter_installation_registry(text, {"1000", "35000"}),
                         '{"1000":1791598517716,"35000":1.79159851766E+12}')


if __name__ == "__main__":
    unittest.main()
