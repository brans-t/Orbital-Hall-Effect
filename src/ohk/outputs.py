"""CSV and headless figure exports from computed data; no physics calculations here."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def require_plotting(formats):
    """Keep CSV-only use independent of Matplotlib; fail before expensive calculations."""
    if set(formats) & {"pdf", "png"}:
        try:
            import matplotlib  # noqa: F401
        except ImportError as exc:
            raise ValueError(
                'PDF/PNG output requires Matplotlib: python -m pip install -e ".[plot]"; '
                "or use --formats csv"
            ) from exc


def _write_csv(path, columns, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def _new_figure(*args, **kwargs):
    # A private, noninteractive canvas avoids changing the caller's global pyplot backend.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(*args, **kwargs)
    FigureCanvasAgg(fig)
    return fig


def _save_figure(fig, directory, stem, formats):
    import matplotlib as mpl

    files = []
    with mpl.rc_context({"pdf.fonttype": 42}):
        for extension in formats:
            if extension in ("pdf", "png"):
                path = directory / f"{stem}.{extension}"
                fig.savefig(path, dpi=180, facecolor="white")
                files.append(path)
    fig.clear()
    return files


def _finish(report, directory, files, inputs=None):
    if inputs is not None:
        path = directory / "input.json"
        path.write_text(json.dumps(inputs, indent=2) + "\n", encoding="utf-8")
        files.append(path)
    path = directory / "report.json"
    files.append(path)
    report.metadata["outputs"] = [str(p.resolve()) for p in files]
    report.metadata["plot_processing"] = "Raw computed values; no smoothing or interpolation"
    report.to_json(path)
    return files


def export_conductivity(report, directory, formats, inputs=None):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    q, files = report.quantities, []
    mus = np.asarray(q.get("fermi_energies", [report.metadata["mu"]]))
    sigma = np.asarray(q.get("fermi_sigma", [q["sigma_oh"]]))
    if "csv" in formats:
        path = directory / "conductivity.csv"
        _write_csv(path, ["mu", "sigma_oh_internal"], zip(mus, sigma))
        files.append(path)
    if set(formats) & {"pdf", "png"}:
        fig = _new_figure(figsize=(7.2, 4.6), layout="constrained")
        ax = fig.subplots()
        order = np.argsort(mus, kind="stable")
        ax.plot(mus[order], sigma[order], "o-", color="#2358AF", ms=2.5, lw=1.2)
        ax.axhline(0, color="0.65", lw=0.7)
        ax.set(
            xlabel="Fermi energy (model energy unit)",
            ylabel=f"OHC coefficient [internal {report.metadata['units']}]",
            title=f"{report.metadata['model']} / {report.metadata['variant']}",
        )
        ax.grid(alpha=0.18)
        files.extend(_save_figure(fig, directory, "conductivity", formats))
    if "convergence" in q:
        study = q["convergence"]
        if "csv" in formats:
            path = directory / "convergence.csv"
            _write_csv(path, ["n_k", "sigma_oh_internal"], zip(study["grid_sizes"], study["sigma"]))
            files.append(path)
        if set(formats) & {"pdf", "png"}:
            fig = _new_figure(figsize=(7.2, 4.6), layout="constrained")
            ax = fig.subplots()
            ax.plot(study["grid_sizes"], study["sigma"], "s--", color="#00775A", ms=4)
            ax.set(
                xlabel="Grid points per periodic direction",
                ylabel=f"OHC coefficient [internal {report.metadata['units']}]",
                title="Mesh study (convergence is not automatically assumed)",
            )
            files.extend(_save_figure(fig, directory, "convergence", formats))
    return _finish(report, directory, files, inputs)


def export_oam_path(report, directory, formats, n_bands, dimension, ticks, labels):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    q, files = report.quantities, []
    rows = np.asarray(q["rows"])
    if "csv" in formats:
        path = directory / "oam.csv"
        _write_csv(path, q["columns"], rows)
        files.append(path)
    if set(formats) & {"pdf", "png"}:
        fig = _new_figure(figsize=(7.2, 8.8), layout="constrained")
        axes = fig.subplots(3, 1, sharex=True)
        start = 1 + dimension
        energies = rows[:, start : start + n_bands]
        values = rows[:, start + n_bands :].reshape(-1, n_bands, n_bands, 2)
        matrices = values[..., 0] + 1j * values[..., 1]
        diagonal = np.diagonal(matrices, axis1=1, axis2=2)
        off = matrices.copy()
        for n in range(n_bands):
            axes[0].plot(
                rows[:, 0], energies[:, n], lw=1, ls=("-", "--", "-.")[n % 3], label=f"Band {n}"
            )
            axes[1].plot(
                rows[:, 0],
                diagonal[:, n].real,
                lw=1,
                ls=("-", "--", "-.")[n % 3],
                label=f"Band {n}",
            )
            off[:, n, n] = 0
        axes[2].plot(rows[:, 0], np.linalg.norm(off, axis=(1, 2)), color="#2358AF")
        axes[0].set(ylabel="Band energy (model unit)", title="All bands along the momentum path")
        axes[1].set(ylabel="Diagonal OAM (internal)")
        axes[2].set(
            ylabel="Off-diagonal OAM\nFrobenius norm (internal)",
            xlabel="Momentum path (Cartesian distance)",
        )
        axes[2].set_xticks(ticks, labels)
        if n_bands <= 8:
            for ax in axes[:2]:
                ax.legend(fontsize=8, ncols=2)
        for ax in axes:
            ax.grid(alpha=0.18)
        files.extend(_save_figure(fig, directory, "oam", formats))
    return _finish(report, directory, files)
