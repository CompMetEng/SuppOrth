"""Run subprocesses with shared helper binaries and stage logs."""

from __future__ import annotations

import platform
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

from . import paths


class CommandFailed(RuntimeError):
    def __init__(self, argv: Sequence[str], returncode: int, log: Path | None):
        self.argv = list(argv)
        self.returncode = returncode
        self.log = log
        location = f"\n  log: {log}" if log else ""
        super().__init__(
            f"command failed with exit code {returncode}:\n"
            f"  {' '.join(str(a) for a in argv)}{location}"
        )


@dataclass
class Result:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    seconds: float
    log: Path | None = None


def tool_env(
    extra: Mapping[str, str] | None = None,
    *,
    threads: int | None = None,
) -> dict[str, str]:
    """Put shared binaries first on PATH and optionally set the OpenMP thread limit."""
    import os

    env = dict(os.environ)
    bin_dir = str(paths.shared_bin())
    current = env.get("PATH", "")
    entries = [entry for entry in current.split(":") if entry != bin_dir]
    env["PATH"] = ":".join([bin_dir, *entries])
    if threads is not None:
        env["OMP_NUM_THREADS"] = str(threads)
    if extra:
        env.update({k: str(v) for k, v in extra.items()})
    return env


def run(
    argv: Sequence[str],
    *,
    cwd: Path | None = None,
    log: Path | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    capture: bool = True,
    timeout: float | None = None,
    threads: int | None = None,
) -> Result:
    argv = [str(a) for a in argv]
    started = time.monotonic()
    if env is None:
        env = tool_env(threads=threads)

    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("w") as handle:
            handle.write(f"$ {' '.join(argv)}\n")
            handle.write(f"# OMP_NUM_THREADS={env.get('OMP_NUM_THREADS', 'unset')}\n")
            handle.flush()
            completed = subprocess.run(
                argv,
                cwd=str(cwd) if cwd else None,
                env=dict(env),
                stdout=handle,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
        stdout = stderr = ""
    else:
        completed = subprocess.run(
            argv,
            cwd=str(cwd) if cwd else None,
            env=dict(env),
            capture_output=capture,
            text=True,
            timeout=timeout,
        )
        stdout = completed.stdout or ""
        stderr = completed.stderr or ""

    result = Result(
        argv=argv,
        returncode=completed.returncode,
        stdout=stdout,
        stderr=stderr,
        seconds=time.monotonic() - started,
        log=log,
    )
    if check and result.returncode != 0:
        raise CommandFailed(argv, result.returncode, log)
    return result


def which(name: str) -> Path | None:
    """Look for a binary in the shared bin first, then the ambient PATH."""
    found = shutil.which(name, path=str(paths.shared_bin()))
    if found is None:
        found = shutil.which(name)
    return Path(found) if found else None


def probe_version(name: str, args: Sequence[str] = ("--version",)) -> str:
    """Best-effort version string, for provenance rather than for logic."""
    binary = which(name)
    if binary is None:
        return "absent"
    try:
        result = run([binary, *args], check=False, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    text = (result.stdout or "") + (result.stderr or "")
    for line in text.splitlines():
        if line.strip():
            return line.strip()[:200]
    return "unknown"


def platform_key() -> str:
    """Coarse platform identifier used to pick prebuilt release assets."""
    system = platform.system().lower()
    machine = platform.machine().lower()
    if machine in {"aarch64", "arm64"}:
        machine = "arm64"
    if system == "darwin":
        system = "macos"
    return f"{system}-{machine}"


def cpu_count(requested: int | None = None) -> int:
    import os

    try:
        available = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        available = os.cpu_count() or 1
    if requested and requested > 0:
        return requested
    return max(1, available - 1)
