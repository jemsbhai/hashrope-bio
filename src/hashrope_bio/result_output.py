"""Shared result output utilities for hashrope-bio benchmarks.

Every benchmark outputs:
    1. results/{name}.json          — latest run (overwritten each time)
    2. results/{name}_{timestamp}.json — archival copy (never overwritten)

Both files contain identical content including full environment metadata
for reproducibility.
"""

from __future__ import annotations

import datetime
import json
import platform
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any


def get_environment() -> dict[str, Any]:
    """Capture full environment metadata for reproducibility."""
    env = {
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "timestamp_local": datetime.datetime.now().isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "os": platform.system(),
        "os_version": platform.version(),
        "node": platform.node(),
    }

    # CPU name (Windows)
    try:
        import subprocess
        result = subprocess.run(
            ["wmic", "cpu", "get", "name"],
            capture_output=True, text=True, timeout=5
        )
        lines = [l.strip() for l in result.stdout.strip().split("\n") if l.strip() and l.strip() != "Name"]
        if lines:
            env["cpu_name"] = lines[0]
    except Exception:
        pass

    # RAM
    try:
        import subprocess
        result = subprocess.run(
            ["wmic", "computersystem", "get", "totalphysicalmemory"],
            capture_output=True, text=True, timeout=5
        )
        lines = [l.strip() for l in result.stdout.strip().split("\n") if l.strip() and l.strip() != "TotalPhysicalMemory"]
        if lines:
            env["ram_bytes"] = int(lines[0])
            env["ram_gb"] = round(int(lines[0]) / (1024**3), 1)
    except Exception:
        pass

    # hashrope version
    try:
        import hashrope
        env["hashrope_version"] = hashrope.__version__
    except Exception:
        pass

    # hashrope-bio version
    try:
        import hashrope_bio
        env["hashrope_bio_version"] = hashrope_bio.__version__
    except Exception:
        pass

    return env


def save_results(
    name: str,
    output_dir: str | Path,
    data: dict[str, Any] | list,
    config: dict[str, Any] | None = None,
    notes: str | None = None,
) -> tuple[Path, Path]:
    """Save benchmark results as JSON with archival timestamped copy.

    Args:
        name: Benchmark name (e.g., "repeat_construction", "resistance_panel").
              Used as the filename stem.
        output_dir: Directory for output files.
        data: The results data — can be a dict, list of dicts, or list of dataclasses.
        config: Optional benchmark configuration (iterations, parameters, etc.).
        notes: Optional free-text notes about this run.

    Returns:
        (latest_path, archive_path) — paths to the two output files.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    env = get_environment()
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")

    # Normalize data — convert dataclasses to dicts
    normalized = _normalize(data)

    payload = {
        "benchmark": name,
        "environment": env,
        "results": normalized,
    }
    if config is not None:
        payload["config"] = config
    if notes is not None:
        payload["notes"] = notes

    # 1. Latest (overwritten)
    latest_path = output_dir / f"{name}.json"
    _write_json(latest_path, payload)

    # 2. Archival (never overwritten)
    archive_path = output_dir / f"{name}_{timestamp}.json"
    _write_json(archive_path, payload)

    print(f"  Results saved:")
    print(f"    Latest:  {latest_path}")
    print(f"    Archive: {archive_path}")

    return latest_path, archive_path


def _normalize(obj: Any) -> Any:
    """Recursively convert dataclasses, sets, bytes to JSON-serializable types."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return {k: _normalize(v) for k, v in asdict(obj).items()}
    if isinstance(obj, list):
        return [_normalize(item) for item in obj]
    if isinstance(obj, dict):
        return {k: _normalize(v) for k, v in obj.items()}
    if isinstance(obj, set):
        return sorted(_normalize(item) for item in obj)
    if isinstance(obj, bytes):
        return obj.decode("utf-8", errors="replace")
    if isinstance(obj, Path):
        return str(obj)
    return obj


def _write_json(path: Path, data: dict) -> None:
    """Write JSON with consistent formatting."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False, default=str)
