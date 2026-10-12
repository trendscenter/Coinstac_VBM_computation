"""Assemble installed MATLAB Runtime products declared by an application's CTF.

Read MathWorks product/component/file manifests; do not guess library dependencies
or edit MATLAB binaries. Intended for a disposable Docker build stage.
"""

import argparse
import json
import os
import re
from pathlib import Path, PurePosixPath
import shutil
import xml.etree.ElementTree as ET


def relative_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError("Unsafe runtime manifest path: " + value)
    return Path(*path.parts)


def filter_installation_registry(text, retained_ids):
    # MATLAB's bootstrap parser is stricter than a general JSON parser. Preserve
    # the installer's exact numeric literals (including scientific notation).
    raw = json.loads(text, parse_float=str, parse_int=str)
    number = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\Z")
    if not isinstance(raw, dict) or not all(isinstance(v, str) and number.fullmatch(v) for v in raw.values()):
        raise ValueError("Unrecognized installedProductData registry")
    return "{" + ",".join(json.dumps(k) + ":" + v for k, v in raw.items() if k in retained_ids) + "}"


def plan_selection(source, requirements):
    required = set(requirements.read_text().split())
    if not required or not all(value.isdigit() for value in required):
        raise ValueError("Expected numeric product IDs in requiredMCRProducts.txt")
    products = {}
    for path in (source / "appdata/products").glob("*.xml"):
        tree = ET.parse(path).getroot()
        product_id = tree.findtext("productNumber")
        if product_id in required or product_id == "1000":
            alternates = tree.find("altDependencies")
            if alternates is not None and len(alternates):
                raise ValueError("Unsupported alternate product dependencies: " + path.name)
            products[path.relative_to(source)] = tree
    found = {tree.findtext("productNumber") for tree in products.values()}
    if required - found:
        raise ValueError("Missing runtime products: " + ", ".join(sorted(required - found)))
    components = set()
    paths = set(products)
    names = {}
    releases = set()
    for tree in products.values():
        names[tree.findtext("productNumber")] = tree.findtext("productName")
        releases.add(tree.findtext("releaseDescription"))
        for component in tree.findall("./dependsOn/component"):
            name = component.findtext("name")
            version = component.findtext("version")
            if not name or not version:
                raise ValueError("Incomplete runtime component definition")
            components.add(name + " " + version)
    if len(releases) != 1:
        raise ValueError("Mixed MATLAB release/update metadata")
    for component in sorted(components):
        manifest = relative_path("appdata/files/" + component + "_manifest.json")
        metadata = relative_path("appdata/components/" + component + ".xml")
        if not (source / manifest).is_file() or not (source / metadata).is_file():
            raise ValueError("Missing component manifest/metadata: " + component)
        members = json.loads((source / manifest).read_text())
        if not isinstance(members, list) or not all(isinstance(p, str) for p in members):
            raise ValueError("Unrecognized file manifest: " + str(manifest))
        paths.update(relative_path(p) for p in members)
        paths.update((manifest, metadata))
    # Installation-generated root files include version, licensing, and notices.
    paths.update(p.relative_to(source) for p in source.iterdir() if p.is_file())
    # Retain symlink targets even when another product owns the target file.
    pending = list(paths)
    while pending:
        relative = pending.pop()
        path = source / relative
        if not path.exists() and not path.is_symlink():
            raise ValueError("Missing declared runtime file: " + str(relative))
        if path.is_symlink():
            target = path.resolve(strict=True)
            try:
                target_relative = target.relative_to(source.resolve())
            except ValueError:
                raise ValueError("Runtime symlink points outside installation: " + str(relative))
            if target_relative not in paths:
                paths.add(target_relative)
                pending.append(target_relative)
    files = [p for p in paths if (source / p).is_file() and not (source / p).is_symlink()]
    report = {
        "method": "installed MathWorks product/component/file manifests",
        "release": next(iter(releases)), "required_product_ids": sorted(required),
        "products": names, "components": len(components),
        "files": len(files), "bytes": sum((source / p).stat().st_size for p in files),
        "paths": sorted(str(p) for p in paths),
    }
    return paths, products, components, report


def assemble(source, requirements, destination):
    paths, products, components, report = plan_selection(source, requirements)
    if destination.exists():
        raise ValueError("Destination must be new: " + str(destination))
    destination.mkdir(parents=True)
    # Create only declared directories, never recursively copy unrelated files.
    for relative in sorted(paths, key=lambda p: (len(p.parts), str(p))):
        src, dst = source / relative, destination / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.is_symlink():
            dst.symlink_to(os.readlink(src))
        elif src.is_dir():
            dst.mkdir(exist_ok=True)
        else:
            shutil.copy2(src, dst)
    metadata = source / "appdata"
    # These installation-generated indexes must describe the subset actually kept.
    retained_ids = set(report["products"]) | {"1000"}
    registry = filter_installation_registry((metadata / "installedProductData.json").read_text(), retained_ids)
    (destination / "appdata/installedProductData.json").write_text(registry)
    indexes = {}
    product_paths = {"/" + str(p) for p in products}
    prodcontents = json.loads((metadata / "prodcontents.json").read_text())
    indexes["prodcontents.json"] = {k: v for k, v in prodcontents.items() if all(p in product_paths for p in v)}
    compcontents = json.loads((metadata / "compcontents.json").read_text())
    indexes["compcontents.json"] = {k: v for k, v in compcontents.items() if k in components}
    for name, value in indexes.items():
        (destination / "appdata" / name).write_text(json.dumps(value, sort_keys=True) + "\n")
    (destination.parent / "runtime-selection.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--requirements", type=Path, required=True)
    parser.add_argument("--destination", type=Path)
    args = parser.parse_args()
    if args.destination:
        report = assemble(args.source, args.requirements, args.destination)
    else:
        report = plan_selection(args.source, args.requirements)[3]
    summary = {k: v for k, v in report.items() if k != "paths"}
    summary["GiB"] = round(report["bytes"] / 1024 ** 3, 3)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
