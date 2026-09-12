#!/usr/bin/env python3
"""Install a GitHub Release wheel with that release's locked runtime dependencies."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request


REPOSITORY = "cheneyyu/PeakForge"
DEFAULT_VERSION = "0.2.3"
DEFAULT_PYTHON = "3.12"
UV_VERSION = "0.11.7"
MARKER = "peakforge-install.json"


def release_version(value: str) -> str:
    version = value.removeprefix("v")
    if not re.fullmatch(r"\d+\.\d+\.\d+(?:(?:a|b|rc)\d+)?", version):
        raise argparse.ArgumentTypeError(
            "Use a release version such as 0.2.3 or v0.2.3."
        )
    return version


def fetch(url: str) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": "PeakForge-release-installer"}
    )
    with urllib.request.urlopen(request, timeout=45) as response:
        return response.read()


def release_assets(release: dict, version: str) -> tuple[dict, dict]:
    if release.get("tag_name") != f"v{version}" or release.get("draft"):
        raise RuntimeError(
            "The GitHub response does not describe the requested published release."
        )
    selected = []
    for name in (
        f"peakforge-{version}.tar.gz",
        f"peakforge-{version}-py3-none-any.whl",
    ):
        matches = [
            asset for asset in release.get("assets", []) if asset["name"] == name
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"Release v{version} must contain exactly one {name} asset."
            )
        asset = matches[0]
        expected_url = (
            f"https://github.com/{REPOSITORY}/releases/download/v{version}/{name}"
        )
        if asset.get("browser_download_url") != expected_url:
            raise RuntimeError(f"Unexpected download URL for {name}.")
        if not re.fullmatch(r"sha256:[0-9a-f]{64}", asset.get("digest") or ""):
            raise RuntimeError(f"GitHub has not provided a SHA-256 digest for {name}.")
        selected.append(asset)
    return selected[0], selected[1]


def extract_lock(archive: bytes, asset: dict, version: str, destination: Path) -> None:
    digest = "sha256:" + hashlib.sha256(archive).hexdigest()
    if digest != asset["digest"]:
        raise RuntimeError(
            "Release source archive checksum mismatch; installation stopped."
        )
    # Read two specific regular files, never unpack arbitrary archive paths or links.
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as source:
        for name in ("pyproject.toml", "uv.lock"):
            members = [
                m
                for m in source.getmembers()
                if m.name == f"peakforge-{version}/{name}"
            ]
            if (
                len(members) != 1
                or not members[0].isfile()
                or members[0].size > 10_000_000
            ):
                raise RuntimeError(
                    f"Release v{version} does not contain a valid {name}."
                )
            with source.extractfile(members[0]) as handle:
                (destination / name).write_bytes(handle.read())


def uv_command() -> list[str]:
    if importlib.util.find_spec("uv") is not None:
        return [sys.executable, "-m", "uv"]
    executable = shutil.which("uv")
    if executable:
        return [executable]
    raise RuntimeError(
        f"Install uv first (https://docs.astral.sh/uv/getting-started/installation/). "
        f"In Colab: %pip install uv=={UV_VERSION}"
    )


def validate_environment(path: Path, version: str, python: str) -> None:
    if path.is_symlink():
        raise RuntimeError("Choose a new environment directory, not a symbolic link.")
    if not path.exists():
        return
    marker = path / MARKER
    if not marker.is_file():
        raise RuntimeError(f"Refusing to change an existing unowned directory: {path}")
    record = json.loads(marker.read_text())
    if (
        record.get("installer") != "peakforge-release"
        or record.get("version") != version
        or record.get("python_request") != python
    ):
        raise RuntimeError(
            "Use a different --env directory for a different release or Python version."
        )
    if not (path / "pyvenv.cfg").is_file() or not (path / "bin/python").is_file():
        raise RuntimeError(
            "The recorded environment is incomplete; choose a new --env directory."
        )


def install(version: str, environment: Path, python: str) -> Path:
    if sys.platform not in ("linux", "darwin"):
        raise RuntimeError(
            "Use Linux, macOS, or WSL2. Native Windows is not supported by this installer."
        )
    validate_environment(environment, version, python)
    environment = environment.resolve()
    uv = uv_command()
    release = json.loads(
        fetch(f"https://api.github.com/repos/{REPOSITORY}/releases/tags/v{version}")
    )
    source, wheel = release_assets(release, version)
    if environment.exists():
        previous = json.loads((environment / MARKER).read_text())
        if (
            previous.get("source_digest") != source["digest"]
            or previous.get("wheel_digest") != wheel["digest"]
        ):
            raise RuntimeError(
                "The release assets have changed since this environment was created. "
                "Its existing installation has been left unchanged."
            )
    record = {
        "installer": "peakforge-release",
        "version": version,
        "python_request": python,
        "source_url": source["browser_download_url"],
        "source_digest": source["digest"],
        "wheel_url": wheel["browser_download_url"],
        "wheel_digest": wheel["digest"],
        "status": "installing",
    }
    print(f"Installing PeakForge {version} into {environment}", flush=True)
    with tempfile.TemporaryDirectory(prefix="peakforge-release-") as temporary:
        work = Path(temporary)
        extract_lock(fetch(source["browser_download_url"]), source, version, work)
        dependencies = subprocess.check_output(
            [
                *uv,
                "export",
                "--project",
                str(work),
                "--frozen",
                "--offline",
                "--extra",
                "macs3",
                "--no-default-groups",
                "--no-emit-project",
                "--no-annotate",
                "--no-header",
                "--format",
                "requirements.txt",
            ],
            text=True,
            cwd=work,
        )
        requirements = (
            f"# PeakForge {version}; runtime dependencies exported from the release uv.lock.\n"
            + dependencies
            + f"\npeakforge[macs3] @ {wheel['browser_download_url']} \\\n"
            + f"    --hash={wheel['digest']}\n"
        )
        if not environment.exists():
            subprocess.run(
                [
                    *uv,
                    "venv",
                    "--python",
                    python,
                    str(environment),
                ],
                check=True,
                cwd=work,
            )
        (environment / MARKER).write_text(json.dumps(record, indent=2) + "\n")
        requirements_path = environment / "peakforge-requirements.txt"
        requirements_path.write_text(requirements)
        shutil.copyfile(work / "uv.lock", environment / "peakforge-uv.lock")
        interpreter = str(environment / "bin/python")
        subprocess.run(
            [
                *uv,
                "pip",
                "sync",
                "--python",
                interpreter,
                "--require-hashes",
                "--strict",
                str(requirements_path),
            ],
            check=True,
            cwd=work,
        )
        subprocess.run(
            [*uv, "pip", "check", "--python", interpreter], check=True, cwd=work
        )
        # -I and a temporary cwd prevent imports from a local source checkout.
        installed_version = subprocess.check_output(
            [
                interpreter,
                "-I",
                "-c",
                "from importlib.metadata import version; import chipdiff; print(version('peakforge'))",
            ],
            text=True,
            cwd=work,
        ).strip()
        if installed_version != version:
            raise RuntimeError(
                f"Installed version {installed_version} differs from requested {version}."
            )
        process_env = dict(
            os.environ,
            PATH=str(environment / "bin") + os.pathsep + os.environ.get("PATH", ""),
        )
        for executable, argument in (
            ("peakforge", "--help"),
            ("multiBamSummary", "--version"),
            ("macs3", "--version"),
        ):
            subprocess.run(
                [str(environment / "bin" / executable), argument],
                check=True,
                cwd=work,
                env=process_env,
                stdout=subprocess.DEVNULL,
            )
        record.update(
            {
                "status": "complete",
                "installed_version": installed_version,
                "lock_sha256": hashlib.sha256(
                    (work / "uv.lock").read_bytes()
                ).hexdigest(),
                "python_version": subprocess.check_output(
                    [interpreter, "--version"], text=True
                ).strip(),
                "uv_version": subprocess.check_output(
                    [*uv, "--version"], text=True
                ).strip(),
                "samtools_path": shutil.which("samtools"),
            }
        )
        (environment / MARKER).write_text(json.dumps(record, indent=2) + "\n")
    print(f"PeakForge {version}, deepTools, and MACS3 are ready.")
    if not shutil.which("samtools"):
        print(
            "samtools is still required for BAM processing. Install it with your system package manager."
        )
    print(f"Activate: source {shlex.quote(str(environment / 'bin/activate'))}")
    print("Then run: peakforge --help")
    return environment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", type=release_version, default=DEFAULT_VERSION)
    parser.add_argument(
        "--env",
        type=Path,
        help="New environment directory (default: ./peakforge-VERSION)",
    )
    parser.add_argument(
        "--python",
        default=DEFAULT_PYTHON,
        help="Python version or interpreter (default: 3.12)",
    )
    args = parser.parse_args()
    try:
        install(
            args.version, args.env or Path(f"peakforge-{args.version}"), args.python
        )
    except (
        RuntimeError,
        OSError,
        ValueError,
        tarfile.TarError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"Installation failed: {error}", file=sys.stderr)
        print(
            "For GitHub connection/rate-limit errors, retry later or configure HTTPS_PROXY.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
