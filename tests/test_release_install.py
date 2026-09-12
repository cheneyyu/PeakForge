from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "peakforge_release_installer",
    Path(__file__).resolve().parents[1] / "scripts/install_release.py",
)
installer = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(installer)


def release_fixture(version="0.2.3"):
    return {
        "tag_name": f"v{version}",
        "draft": False,
        "assets": [
            {
                "name": name,
                "digest": "sha256:" + "a" * 64,
                "browser_download_url": f"https://github.com/cheneyyu/PeakForge/releases/download/v{version}/{name}",
            }
            for name in (
                f"peakforge-{version}.tar.gz",
                f"peakforge-{version}-py3-none-any.whl",
            )
        ],
    }


def archive_fixture(*, link=False, missing=False):
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in ("pyproject.toml", "uv.lock", "../../outside.txt"):
            if missing and name == "uv.lock":
                continue
            content = b"example\n"
            info = tarfile.TarInfo(f"peakforge-0.2.3/{name}")
            info.size = len(content)
            if link and name == "uv.lock":
                info.type = tarfile.SYMTYPE
                info.linkname = "/etc/passwd"
            archive.addfile(info, io.BytesIO(content))
    raw = output.getvalue()
    return raw, {"digest": "sha256:" + hashlib.sha256(raw).hexdigest()}


class ReleaseInstallTests(unittest.TestCase):
    def test_version_is_explicit_and_not_a_branch_or_path(self):
        self.assertEqual(installer.release_version("v0.2.3"), "0.2.3")
        self.assertEqual(installer.release_version("0.2.4rc1"), "0.2.4rc1")
        for value in ("main", "latest", "../0.2.3", "0.2.3;echo no", "v0.2"):
            with (
                self.subTest(value=value),
                self.assertRaises(argparse.ArgumentTypeError),
            ):
                installer.release_version(value)

    def test_selects_release_wheel_and_source_with_digests(self):
        source, wheel = installer.release_assets(release_fixture(), "0.2.3")
        self.assertTrue(source["name"].endswith(".tar.gz"))
        self.assertTrue(wheel["name"].endswith(".whl"))

    def test_rejects_mismatched_release_and_missing_or_untrusted_assets(self):
        mutations = [
            lambda data: data.update(tag_name="v0.2.2"),
            lambda data: data.update(draft=True),
            lambda data: data["assets"].pop(),
            lambda data: data["assets"].append(data["assets"][0]),
            lambda data: data["assets"][0].update(digest=None),
            lambda data: data["assets"][0].update(
                browser_download_url="http://other.invalid/pkg"
            ),
        ]
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                release = release_fixture()
                mutate(release)
                with self.assertRaises(RuntimeError):
                    installer.release_assets(release, "0.2.3")

    def test_extracts_only_lock_and_metadata(self):
        raw, source = archive_fixture()
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "metadata"
            destination.mkdir()
            installer.extract_lock(raw, source, "0.2.3", destination)
            self.assertEqual(
                sorted(p.name for p in destination.iterdir()),
                ["pyproject.toml", "uv.lock"],
            )
            self.assertFalse((Path(tmp) / "outside.txt").exists())

    def test_refuses_bad_checksum_missing_lock_and_links(self):
        for settings in ({"link": True}, {"missing": True}, {}):
            with self.subTest(settings=settings), tempfile.TemporaryDirectory() as tmp:
                raw, source = archive_fixture(**settings)
                if not settings:
                    source["digest"] = "sha256:" + "0" * 64
                with self.assertRaises(RuntimeError):
                    installer.extract_lock(raw, source, "0.2.3", Path(tmp))

    def test_does_not_touch_existing_user_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "keep.txt").write_text("unchanged")
            with self.assertRaises(RuntimeError):
                installer.validate_environment(root, "0.2.3", "3.12")
            self.assertEqual((root / "keep.txt").read_text(), "unchanged")

    def test_rerun_requires_matching_owned_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "pyvenv.cfg").touch()
            (root / "bin").mkdir()
            (root / "bin/python").touch()
            (root / installer.MARKER).write_text(
                json.dumps(
                    {
                        "installer": "peakforge-release",
                        "version": "0.2.3",
                        "python_request": "3.12",
                    }
                )
            )
            installer.validate_environment(root, "0.2.3", "3.12")
            with self.assertRaises(RuntimeError):
                installer.validate_environment(root, "0.2.4", "3.12")
            with self.assertRaises(RuntimeError):
                installer.validate_environment(root, "0.2.3", "3.13")
            link = root / "linked"
            link.symlink_to(root, target_is_directory=True)
            with self.assertRaises(RuntimeError):
                installer.validate_environment(link, "0.2.3", "3.12")

    def test_uv_missing_is_actionable(self):
        with (
            patch.object(installer.importlib.util, "find_spec", return_value=None),
            patch.object(installer.shutil, "which", return_value=None),
        ):
            with self.assertRaisesRegex(RuntimeError, "Install uv first"):
                installer.uv_command()

    def test_changed_release_assets_do_not_modify_an_existing_installation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "pyvenv.cfg").touch()
            (root / "bin").mkdir()
            (root / "bin/python").touch()
            marker = root / installer.MARKER
            original = json.dumps(
                {
                    "installer": "peakforge-release",
                    "version": "0.2.3",
                    "python_request": "3.12",
                    "source_digest": "sha256:" + "0" * 64,
                    "wheel_digest": "sha256:" + "a" * 64,
                    "status": "complete",
                }
            )
            marker.write_text(original)
            with (
                patch.object(installer, "uv_command", return_value=["uv"]),
                patch.object(
                    installer,
                    "fetch",
                    return_value=json.dumps(release_fixture()).encode(),
                ),
                patch.object(installer.subprocess, "run") as run,
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "release assets have changed"
                ):
                    installer.install("0.2.3", root, "3.12")
                run.assert_not_called()
            self.assertEqual(marker.read_text(), original)

if __name__ == "__main__":
    unittest.main()
