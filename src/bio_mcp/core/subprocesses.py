"""Safe local subprocess helpers for lightweight CLI wrappers."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from os import access, environ, X_OK
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class CommandResult:
    """Captured subprocess result."""

    args: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str


@dataclass(frozen=True)
class CommandStatus:
    """Local command availability and version information."""

    available: bool
    path: str | None
    version: str | None
    resolution_source: str = "missing"
    error: str | None = None


class CommandTimeoutError(RuntimeError):
    """Raised when a local command exceeds its timeout."""


def summarize_stderr(stderr: str, *, max_chars: int = 1_000) -> str:
    """Return a compact stderr summary suitable for JSON tool output."""

    cleaned = "\n".join(line.rstrip() for line in stderr.splitlines() if line.strip())
    if len(cleaned) <= max_chars:
        return cleaned
    return cleaned[: max_chars - 3] + "..."


def run_command(
    args: Sequence[str],
    *,
    cwd: str | Path,
    timeout_sec: int,
) -> CommandResult:
    """Run a local command without a shell and capture text output."""

    try:
        completed = subprocess.run(
            list(args),
            cwd=str(cwd),
            text=True,
            capture_output=True,
            timeout=timeout_sec,
            shell=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise CommandTimeoutError(
            f"Command timed out after {timeout_sec} seconds: {args[0]}"
        ) from exc

    return CommandResult(
        args=tuple(str(arg) for arg in args),
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
    )


def detect_command(
    command: str,
    version_args: Sequence[str],
    *,
    env_var: str | None = None,
    timeout_sec: int = 5,
) -> CommandStatus:
    """Detect a command path and best-effort version string."""

    path: str | None = None
    resolution_source = "missing"
    if env_var:
        env_path = environ.get(env_var)
        if env_path:
            resolution_source = "env_override"
            candidate = Path(env_path)
            if not candidate.is_file():
                return CommandStatus(
                    available=False,
                    path=env_path,
                    version=None,
                    resolution_source=resolution_source,
                    error=(
                        f"Dependency error: {env_var} points to '{env_path}', "
                        "but that file does not exist."
                    ),
                )
            if not access(candidate, X_OK):
                return CommandStatus(
                    available=False,
                    path=env_path,
                    version=None,
                    resolution_source=resolution_source,
                    error=(
                        f"Dependency error: {env_var} points to '{env_path}', "
                        "but it is not executable."
                    ),
                )
            path = str(candidate)

    if path is None:
        path = shutil.which(command)
        if path is not None:
            resolution_source = "PATH"

    if path is None:
        return CommandStatus(
            available=False,
            path=None,
            version=None,
            resolution_source="missing",
            error=f"Dependency error: required binary '{command}' was not found on PATH.",
        )

    try:
        result = run_command([path, *version_args], cwd=Path.cwd(), timeout_sec=timeout_sec)
    except CommandTimeoutError as exc:
        return CommandStatus(
            available=True,
            path=path,
            version=None,
            resolution_source=resolution_source,
            error=str(exc),
        )
    except OSError as exc:
        return CommandStatus(
            available=True,
            path=path,
            version=None,
            resolution_source=resolution_source,
            error=f"Could not run '{command}' to detect version: {exc}",
        )

    combined_output = "\n".join(part for part in (result.stdout, result.stderr) if part.strip())
    version = next((line.strip() for line in combined_output.splitlines() if line.strip()), None)
    error = None if result.returncode == 0 else summarize_stderr(result.stderr) or result.stdout.strip()
    return CommandStatus(
        available=True,
        path=path,
        version=version,
        resolution_source=resolution_source,
        error=error,
    )
