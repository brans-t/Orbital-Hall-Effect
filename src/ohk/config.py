"""Validated, flat JSON inputs for the conductivity command."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

FORMATS = ("csv", "pdf", "png")
RUN_KEYS = {
    "model",
    "t1",
    "t2",
    "mass",
    "phi",
    "nearest_neighbour_distance",
    "n_k",
    "mu",
    "eta",
    "variant",
    "hr",
    "lattice",
    "fermi",
    "fermi_range",
    "convergence",
    "output_dir",
    "formats",
}


def load_run_config(path):
    """Paths in a config are relative to that file, independent of the working directory."""
    path = Path(path).resolve()
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(data, dict):
        raise ValueError("Config must be a JSON object")
    unknown = set(data) - RUN_KEYS
    if unknown:
        raise ValueError(f"Unknown config keys: {', '.join(sorted(unknown))}")
    for key in ("hr", "lattice", "output_dir"):
        if key in data:
            if not isinstance(data[key], str) or not data[key]:
                raise ValueError(f"{key} must be a nonempty path string")
            data[key] = (path.parent / data[key]).resolve()
    return data


def validate_run(args):
    for key in ("t1", "t2", "mass", "phi", "nearest_neighbour_distance", "mu", "eta"):
        value = getattr(args, key)
        if value is None and key == "eta":
            continue
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ValueError(f"{key} must be a finite number")
    if args.nearest_neighbour_distance <= 0:
        raise ValueError("nearest_neighbour_distance must be positive")
    if args.eta is not None and args.eta <= 0:
        raise ValueError("eta must be positive")
    if isinstance(args.n_k, bool) or not isinstance(args.n_k, int) or args.n_k < 1:
        raise ValueError("n_k must be a positive integer")
    if args.model not in ("haldane", "wannier90"):
        raise ValueError("model must be haldane or wannier90")
    if args.variant not in ("covariant", "spectral", "pezo", "aca"):
        raise ValueError("variant must be covariant, spectral, pezo or aca")
    for key in ("fermi", "convergence", "fermi_range"):
        value = getattr(args, key)
        if value is not None:
            if not isinstance(value, list) or not value:
                raise ValueError(f"{key} must be a nonempty array")
            if any(isinstance(x, bool) or not isinstance(x, (int, float)) for x in value):
                raise ValueError(f"{key} must contain numbers")
            if not np.isfinite(value).all():
                raise ValueError(f"{key} must contain finite numbers")
    if args.convergence is not None and any(int(n) != n or n < 1 for n in args.convergence):
        raise ValueError("convergence must contain positive integers")
    if args.fermi is not None and args.fermi_range is not None:
        raise ValueError("Choose either fermi or fermi_range")
    if args.fermi_range is not None:
        if len(args.fermi_range) != 3:
            raise ValueError("fermi_range must be [start, stop, count]")
        start, stop, count = args.fermi_range
        if start >= stop or int(count) != count or count < 2:
            raise ValueError("fermi_range requires start < stop and an integer count >= 2")


def output_formats(value):
    if value is None:
        return ["csv", "pdf"]
    if not isinstance(value, list) or not value or any(x not in FORMATS for x in value):
        raise ValueError("formats must be a nonempty array containing csv, pdf or png")
    return list(dict.fromkeys(value))


def resolved_input(args):
    """A reusable JSON config with absolute input/output paths."""
    data = {
        key: getattr(args, key) for key in sorted(RUN_KEYS) if getattr(args, key, None) is not None
    }
    data["model"] = "wannier90" if args.hr is not None else args.model
    for key in ("hr", "lattice", "output_dir"):
        if key in data:
            data[key] = str(Path(data[key]).resolve())
    return data
