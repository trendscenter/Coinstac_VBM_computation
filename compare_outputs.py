#!/usr/bin/env python3
"""Compare two COINSTAC VBM output zip files.

Compare decompressed file contents, independently of ZIP timestamps/compression.
All files are strict by default, including MAT files.
"""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile


DEFAULT_IGNORE_PATTERNS = (
    "__MACOSX/*",
    "*/._*",
    "._*",
)


@dataclass(frozen=True)
class EntryInfo:
    sha256: str
    size: int


def _normalize_zip_path(name: str, strip_prefixes: list[str]) -> str:
    normalized = name.replace("\\", "/").lstrip("./")
    for prefix in strip_prefixes:
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]
            break
    return normalized


def _should_ignore(path: str, patterns: tuple[str, ...]) -> bool:
    return any(fnmatch.fnmatch(path, pattern) for pattern in patterns)


def build_manifest(
    zip_path: Path,
    strip_prefixes: list[str],
    ignore_patterns: tuple[str, ...],
) -> dict[str, EntryInfo]:
    manifest: dict[str, EntryInfo] = {}
    with ZipFile(zip_path, "r") as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            normalized = _normalize_zip_path(info.filename, strip_prefixes)
            if not normalized or _should_ignore(normalized, ignore_patterns):
                continue
            if normalized in manifest:
                raise ValueError(f"Duplicate normalized archive path: {normalized}")
            digest = hashlib.sha256()
            with archive.open(info) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            manifest[normalized] = EntryInfo(sha256=digest.hexdigest(), size=info.file_size)
    return manifest


def compare_manifests(
    baseline: dict[str, EntryInfo],
    candidate: dict[str, EntryInfo],
    tolerate_patterns: tuple[str, ...],
) -> tuple[set[str], set[str], list[str], list[str], list[str]]:
    baseline_keys = set(baseline.keys())
    candidate_keys = set(candidate.keys())
    common = baseline_keys & candidate_keys

    only_in_baseline = baseline_keys - candidate_keys
    only_in_candidate = candidate_keys - baseline_keys

    changed_strict: list[str] = []
    changed_tolerated: list[str] = []

    for path in sorted(common):
        if baseline[path] == candidate[path]:
            continue
        if _should_ignore(path, tolerate_patterns):
            changed_tolerated.append(path)
        else:
            changed_strict.append(path)

    matched = sorted(path for path in common if baseline[path] == candidate[path])
    return only_in_baseline, only_in_candidate, changed_strict, changed_tolerated, matched


def _print_list(title: str, items: list[str] | set[str], max_items: int) -> None:
    items_sorted = sorted(items)
    print(f"{title}: {len(items_sorted)}")
    for path in items_sorted[:max_items]:
        print(f"  - {path}")
    if len(items_sorted) > max_items:
        print(f"  ... ({len(items_sorted) - max_items} more)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two VBM output zip archives.")
    parser.add_argument("baseline_zip", type=Path, help="Known-good output zip path.")
    parser.add_argument("candidate_zip", type=Path, help="New output zip path.")
    parser.add_argument(
        "--strip-prefix",
        action="append",
        default=[],
        help="Strip a path prefix from entries before comparison (can repeat).",
    )
    parser.add_argument(
        "--ignore",
        action="append",
        default=[],
        help="Glob path(s) to ignore entirely during comparison.",
    )
    parser.add_argument(
        "--tolerate",
        action="append",
        default=[],
        help="Explicit glob path(s) to report but not fail on (none by default).",
    )
    parser.add_argument(
        "--max-items",
        type=int,
        default=10,
        help="Max paths to print per category.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.baseline_zip.exists():
        print(f"Baseline zip not found: {args.baseline_zip}", file=sys.stderr)
        return 2
    if not args.candidate_zip.exists():
        print(f"Candidate zip not found: {args.candidate_zip}", file=sys.stderr)
        return 2

    ignore_patterns = tuple(DEFAULT_IGNORE_PATTERNS + tuple(args.ignore))
    strip_prefixes = [prefix.strip("/") + "/" for prefix in args.strip_prefix if prefix]
    tolerate_patterns = tuple(args.tolerate)

    baseline_manifest = build_manifest(args.baseline_zip, strip_prefixes, ignore_patterns)
    candidate_manifest = build_manifest(args.candidate_zip, strip_prefixes, ignore_patterns)
    if not baseline_manifest or not candidate_manifest:
        print("Cannot compare empty output archives.", file=sys.stderr)
        return 2

    only_base, only_candidate, changed_strict, changed_tolerated, matched = compare_manifests(
        baseline_manifest,
        candidate_manifest,
        tolerate_patterns,
    )

    print("Comparison summary")
    print(f"- Baseline files: {len(baseline_manifest)}")
    print(f"- Candidate files: {len(candidate_manifest)}")
    print(f"- Exact matches: {len(matched)}")
    print(f"- Strict differences: {len(changed_strict)}")
    print(f"- Tolerated differences: {len(changed_tolerated)}")
    print(f"- Missing from candidate: {len(only_base)}")
    print(f"- Added in candidate: {len(only_candidate)}")
    print()

    _print_list("Missing from candidate", only_base, args.max_items)
    _print_list("Added in candidate", only_candidate, args.max_items)
    _print_list("Changed (strict)", changed_strict, args.max_items)
    _print_list("Changed (tolerated)", changed_tolerated, args.max_items)

    if only_base or only_candidate or changed_strict:
        print("\nResult: FAIL")
        return 1

    print("\nResult: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
