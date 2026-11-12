"""Pure-data reports, including lossless JSON encoding of complex matrices."""

from __future__ import annotations
from dataclasses import asdict, dataclass, field, is_dataclass
from pathlib import Path
from typing import Any
import json
import numpy as np


def _json_value(value):
    if is_dataclass(value):
        return _json_value(asdict(value))
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, complex):
        return {"real": _json_value(value.real), "imag": _json_value(value.imag)}
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


@dataclass
class Report:
    title: str = ""
    quantities: dict[str, Any] = field(default_factory=dict)
    checks: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def add(self, name, value):
        self.quantities[name] = value
        return self

    def warn(self, message):
        self.warnings.append(message)
        return self

    def failed_checks(self, tol=1e-8):
        return {
            name: value
            for name, value in self.checks.items()
            if not np.isfinite(value) or abs(value) > tol
        }

    def to_dict(self):
        return _json_value(asdict(self))

    def to_json(self, path=None):
        serialized = json.dumps(self.to_dict(), indent=2, ensure_ascii=False, allow_nan=False)
        if path is not None:
            Path(path).write_text(serialized + "\n", encoding="utf-8")
        return serialized

    def summary(self):
        lines = [self.title]
        for key, value in self.quantities.items():
            if isinstance(value, np.ndarray):
                lines.append(f"  {key}: array{value.shape}")
            else:
                lines.append(f"  {key}: {value}")
        for key, value in self.checks.items():
            lines.append(
                f"  [{'FAIL' if key in self.failed_checks() else 'PASS'}] {key}: {value:.3g}"
            )
        lines.extend(f"  Warning: {message}" for message in self.warnings)
        return "\n".join(lines)
