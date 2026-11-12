"""Recompute the three main Haldane figures, source tables and an audit manifest.

Run: python scripts/reproduce_paper.py --output output/paper-reproduction
Requires the plot extra. The default is a 240x240 Bloch mesh and an independent 60x60
sparse real-space torus. No author data are digitized, fitted, smoothed or rescaled to fit.
The eV/Angstrom convention is an explicit, unconfirmed physical-scale assumption.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import ohk
from ohk import HaldaneModel, EigenBasis, CovariantProjector, full_matrix
from ohk.units import UnitSystem
from ohk.validation.reproduction import bloch_response, finite_pbc_response, uniform_response

COLORS = ["#222222", "#00775A", "#2358AF", "#B51E48"]
STYLES = ["-", "--", "-.", ":"]
PHASES = [0.0, 0.15, 0.25, 0.5]
SOURCE = "https://doi.org/10.1103/pd3s-q2z7"


def paper_model(phi_over_pi):
    # Paper positive phi has C=+1, while this library's positive phi has C=-1.
    return HaldaneModel(t1=1.0, t2=0.15, M=0.4, phi=-phi_over_pi * np.pi)


def write_csv(path, header, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)


def export(fig, directory, name):
    for extension in ("png", "svg"):
        fig.savefig(directory / f"{name}.{extension}", dpi=300, facecolor="white")
    plt.close(fig)


def figure_one(directory, finite_size, scale):
    coordinate = np.linspace(-1, 1, 801)
    model = paper_model(0)
    corner = np.array([1 / 3, 2 / 3]) @ model.spec.reciprocal_vectors()
    kpoints = -coordinate[:, None] * corner
    # Every marker is an allowed momentum: increments (1,2)/N along Gamma-K.
    count = finite_size // 3
    marker_x = np.linspace(-1, 1, 2 * count + 1)
    marker_frac = -marker_x[:, None] * np.array([1 / 3, 2 / 3])
    indices = np.rint(marker_frac * finite_size).astype(int) % finite_size
    flat_indices = np.ravel_multi_index(indices.T, (finite_size, finite_size))
    fig, axes = plt.subplots(
        3, 1, figsize=(7.4, 9.4), layout="constrained", gridspec_kw={"height_ratios": [1, 1, 0.85]}
    )
    rows, finite_rows, maxima = [], [], {}
    maximum_error = 0.0
    for phase, color, style in zip(PHASES, COLORS, STYLES):
        model = paper_model(phase)
        dense = bloch_response(model, kpoints)
        finite = finite_pbc_response(model, finite_size)
        direct = finite.oam[flat_indices]
        expected = bloch_response(model, -marker_x[:, None] * corner).oam
        error = float(np.max(abs(abs(direct) - abs(expected))))
        maximum_error = max(maximum_error, error)
        off = abs(dense.oam[:, 0, 1]) * scale
        diagonal = dense.oam[:, 0, 0].real * scale
        off_f = abs(direct[:, 0, 1]) * scale
        diag_f = direct[:, 0, 0].real * scale
        label = rf"$\phi={phase:g}\pi$" if phase else r"$\phi=0$"
        axes[0].plot(coordinate, off, color=color, ls=style, label=label)
        axes[0].plot(marker_x, off_f, "o", color=color, mfc="white", ms=3.6, mew=0.8)
        for ax in axes[1:]:
            ax.plot(coordinate, diagonal, color=color, ls=style)
            ax.plot(marker_x, diag_f, "o", color=color, mfc="white", ms=3.6, mew=0.8)
        maxima[str(phase)] = {
            "max_abs_L12_hbar": float(off.max()),
            "L11_range_hbar": [float(diagonal.min()), float(diagonal.max())],
        }
        rows.extend(
            zip(
                np.repeat(phase, len(coordinate)),
                coordinate,
                kpoints[:, 0],
                kpoints[:, 1],
                dense.energies[:, 0],
                dense.energies[:, 1],
                abs(dense.oam[:, 0, 1]),
                dense.oam[:, 0, 0].real,
                off,
                diagonal,
            )
        )
        finite_rows.extend(zip(np.repeat(phase, len(marker_x)), marker_x, off_f, diag_f))
    axes[0].set(ylabel=r"$|L_{12}|/\hbar$", title="(a) Off-diagonal OAM")
    axes[0].legend(ncol=2, loc="upper center", frameon=False)
    axes[1].set(
        ylabel=r"$L_{11}/\hbar$", ylim=(-1, 1), title="(b) Diagonal OAM: paper display window"
    )
    axes[1].text(
        0.5,
        0.05,
        "Valley extrema exceed this window; all values appear in (c).",
        transform=axes[1].transAxes,
        ha="center",
        va="bottom",
        fontsize=8,
    )
    axes[2].set_yscale("symlog", linthresh=0.1, linscale=0.6)
    axes[2].set(
        ylabel=r"$L_{11}/\hbar$", title="(c) Full range (symmetric log; linear within +/-0.1)"
    )
    for ax in axes:
        ax.set(xlim=(-1, 1), xticks=[-1, 0, 1], xticklabels=["K", r"$\Gamma$", "K'"])
        ax.axvline(0, color="0.65", ls="--", lw=0.7)
        ax.axhline(0, color="0.75", lw=0.6, zorder=0)
    fig.suptitle("Figure 1 reproduction: orbital angular momentum", fontsize=14)
    fig.supxlabel(
        f"Lines: projected Bloch formula; circles: independent {finite_size} x {finite_size} PBC torus"
    )
    export(fig, directory, "fig1_oam")
    write_csv(
        directory / "fig1_bloch.csv",
        [
            "paper_phi_over_pi",
            "path_coordinate",
            "kx",
            "ky",
            "eps0",
            "eps1",
            "abs_L12_internal",
            "L11_internal",
            "abs_L12_hbar",
            "L11_hbar",
        ],
        rows,
    )
    write_csv(
        directory / "fig1_finite.csv",
        ["paper_phi_over_pi", "path_coordinate", "abs_L12_hbar", "L11_hbar"],
        finite_rows,
    )
    return {"max_finite_bloch_abs_error_internal": maximum_error, "extrema": maxima}


def figure_two(directory, n_dense, finite_size):
    critical = HaldaneModel().phase_boundary() / np.pi
    singular = np.array([-1 + critical, -critical, critical, 1 - critical])
    regular = np.linspace(-1, 1, 201)
    phases = np.sort(np.unique(np.r_[regular, singular]))
    curves = {finite_size: [], n_dense // 2: [], n_dense: []}
    rows = []
    started = time.perf_counter()
    for index, phase in enumerate(phases):
        if np.any(np.isclose(phase, singular, atol=1e-12, rtol=0)):
            for values in curves.values():
                values.append(float("nan"))
            for n in curves:
                rows.append([phase, n, "undefined_gap_closing", np.nan, np.nan, np.nan])
            continue
        model = paper_model(phase)
        for n in curves:
            response = (
                finite_pbc_response(model, n) if n == finite_size else uniform_response(model, n)
            )
            value = response.sigma()  # Explicitly one occupied band, not fixed mu=0.
            curves[n].append(value)
            rows.append(
                [
                    phase,
                    n,
                    "finite_PBC" if n == finite_size else "Bloch",
                    -3 * model.t2 * np.cos(phase * np.pi),
                    value,
                    response.min_gap,
                ]
            )
        if index % 25 == 0:
            print(
                f"Phase sweep {index + 1}/{len(phases)}; {time.perf_counter() - started:.1f} s",
                flush=True,
            )
    write_csv(
        directory / "fig2_phase.csv",
        [
            "paper_phi_over_pi",
            "mesh_n",
            "method",
            "gap_center_mu",
            "sigma_internal_e_over_2pi",
            "minimum_gap",
        ],
        rows,
    )
    return phases, curves, singular


def plot_two(directory, phases, curves, singular, n_dense, finite_size, scale):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5), layout="constrained")
    marker_mask = np.isclose(phases * 20, np.rint(phases * 20), atol=1e-9)
    for ax in axes:
        for n, color, style in [(n_dense, "#184C93", "-"), (n_dense // 2, "#697585", "--")]:
            ax.plot(
                phases,
                -np.asarray(curves[n]) * scale,
                color=color,
                ls=style,
                label=f"Bloch {n} x {n}",
            )
        ax.plot(
            phases[marker_mask],
            -np.asarray(curves[finite_size])[marker_mask] * scale,
            "o",
            ms=4,
            mfc="#F0BE43",
            mec="#6A4710",
            mew=0.6,
            label=f"Real-space PBC {finite_size} x {finite_size}",
        )
        for x in singular:
            ax.axvline(x, color="0.5", ls=":", lw=0.8)
        ax.set(
            xlim=(-1, 1),
            xlabel=r"Paper $\phi/\pi$",
            ylabel=r"$-\sigma_{xy}^{\rm OH}\;[e/(2\pi)]$",
            xticks=[-1, -0.5, 0, 0.5, 1],
        )
    axes[0].set(ylim=(0, 1.3), title="(a) Paper display window")
    axes[0].legend(loc="lower center", fontsize=8, frameon=False)
    # Show every computed finite-torus sample in the full-range panel, including large
    # mesh artifacts near the transition. The larger markers above are the paper-style subset.
    axes[1].plot(
        phases,
        -np.asarray(curves[finite_size]) * scale,
        ".",
        color="#A36300",
        ms=3,
        label=f"PBC {finite_size}: all phase samples",
        zorder=2,
    )
    axes[1].legend(loc="upper center", fontsize=8, frameon=False)
    axes[1].set_yscale("log")
    axes[1].set(title="(b) Full computed range (log scale)")
    fig.suptitle("Figure 2 reproduction: phase sweep at half filling", fontsize=14)
    fig.supxlabel(
        "Dashed vertical lines: exact gap closings (undefined); mesh dependence is retained.",
        fontsize=9,
    )
    export(fig, directory, "fig2_phase")


def figure_three(directory, n_dense, finite_size, scale, energy_ev=1.0):
    model = paper_model(0.5)
    fine = uniform_response(model, n_dense)
    finite = finite_pbc_response(model, finite_size)
    ef = np.linspace(-2, 2, 801)
    fine_sigma, finite_sigma = fine.fermi_scan(ef), finite.fermi_scan(ef)
    physical_ef = ef * energy_ev
    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    ax.plot(physical_ef, -fine_sigma * scale, color="#184C93", label=f"Bloch {n_dense} x {n_dense}")
    ax.plot(
        physical_ef[::16],
        -finite_sigma[::16] * scale,
        "o",
        mfc="#F0BE43",
        mec="#6A4710",
        ms=4.2,
        mew=0.6,
        label=f"Independent PBC {finite_size} x {finite_size}",
    )
    ax.axvline(0, color="0.5", ls="--", lw=0.8)
    ax.set(
        xlim=(-2 * energy_ev, 2 * energy_ev),
        ylim=(0, max(0.31, float(np.max(-fine_sigma * scale)) * 1.1)),
        xlabel=r"$E_F$ [eV]",
        ylabel=r"$-\sigma_{xy}^{\rm OH}\;[e/(2\pi)]$",
        title=r"Figure 3 reproduction: $\phi=\pi/2$",
    )
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    fig.supxlabel("No broadening or smoothing; exact zero-temperature occupations.", fontsize=9)
    export(fig, directory, "fig3_fermi")
    write_csv(
        directory / "fig3_fermi.csv",
        [
            "fermi_energy_eV",
            "bloch_sigma_internal",
            "finite_sigma_internal",
            "bloch_minus_sigma_e_over_2pi",
            "finite_minus_sigma_e_over_2pi",
        ],
        zip(physical_ef, fine_sigma, finite_sigma, -fine_sigma * scale, -finite_sigma * scale),
    )
    return {
        "bloch_plateau": -fine.sigma(mu=0) * scale,
        "finite_plateau": -finite.sigma(mu=0) * scale,
        "max_curve_mesh_difference": float(np.max(abs(fine_sigma - finite_sigma))) * scale,
    }


def convergence(directory, scale, n_dense):
    sizes = sorted(set([30, 60, 120, n_dense, 2 * n_dense]))
    rows, errors = [], []
    for phase in [0.0, 0.1, 0.15, 0.2, 0.25, 0.5]:
        last = None
        for n in sizes:
            response = uniform_response(paper_model(phase), n)
            value = -response.sigma() * scale
            rows.append([phase, n, value, np.nan if last is None else value - last])
            last = value
        if phase == 0.5:
            errors.append(abs(rows[-1][3]))
    write_csv(
        directory / "convergence.csv",
        ["paper_phi_over_pi", "mesh_n", "minus_sigma_e_over_2pi", "change_from_previous"],
        rows,
    )
    return {"plateau_last_mesh_change": errors[0], "sizes": sizes}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("output/paper-reproduction"))
    parser.add_argument(
        "--n-k",
        type=int,
        default=240,
        help="Fine Bloch mesh; multiple of 6, greater than twice finite-size",
    )
    parser.add_argument("--finite-size", type=int, default=60)
    parser.add_argument("--energy-ev", type=float, default=1.0)
    parser.add_argument("--length-angstrom", type=float, default=1.0)
    args = parser.parse_args(argv)
    if (
        args.finite_size % 3
        or args.finite_size < 3
        or args.n_k % 6
        or args.n_k <= 2 * args.finite_size
    ):
        parser.error(
            "finite-size must be >=3 and divisible by 3; n-k must be divisible by 6 and >2*finite-size"
        )
    args.output.mkdir(parents=True, exist_ok=True)
    units = UnitSystem(energy=args.energy_ev, length=args.length_angstrom)
    started = time.perf_counter()
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.7,
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
    }
    with plt.rc_context(style):
        one = figure_one(args.output, args.finite_size, units.oam_scale)
        phases, curves, singular = figure_two(args.output, args.n_k, args.finite_size)
        plot_two(args.output, phases, curves, singular, args.n_k, args.finite_size, units.oam_scale)
        three = figure_three(
            args.output, args.n_k, args.finite_size, units.oam_scale, args.energy_ev
        )
        conv = convergence(args.output, units.oam_scale, args.n_k)
    # Validate one sampled point against the established single-point public API.
    model, k = paper_model(0.25), np.array([0.7, 1.1])
    point = bloch_response(model, k[None])
    reference = full_matrix(CovariantProjector(), EigenBasis.at(model, k))[0]
    manifest = {
        "source": SOURCE,
        "software_version": ohk.__version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_pdf_sha256": hashlib.sha256(
            Path("docs/Orbital Hall effect from the Wannier function approach.pdf").read_bytes()
        ).hexdigest(),
        "physical_scale_assumption": {
            "energy_eV": args.energy_ev,
            "length_angstrom": args.length_angstrom,
            "oam_factor": units.oam_scale,
            "confirmed_by_paper_text": False,
        },
        "phase_mapping": "library_phi = -paper_phi; paper positive phi has lowest-band C=+1",
        "occupation_fig2": "lowest band occupied at every k (half filling); mu follows gap center",
        "occupation_fig3": "eps < E_F at zero temperature",
        "broadening": 0.0,
        "mesh_bloch": args.n_k,
        "mesh_finite": args.finite_size,
        "phase_boundaries_over_pi": singular.tolist(),
        "fig1": one,
        "fig3": three,
        "convergence": conv,
        "single_point_crosscheck": float(np.max(abs(point.oam[0] - reference))),
        "transformations": [
            "paper flux sign mapped to library convention",
            "physical unit restoration",
            "plot minus sigma as in paper",
            "no fitting, smoothing or digitization",
        ],
        "singular_policy": "Exact gap-closing phase samples retained as NaN gaps in CSV and curves",
        "display_windows": "Fig1b and Fig2a match paper bounds; full ranges are shown in adjacent panels",
        "runtime_seconds": time.perf_counter() - started,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "fig1": one,
                "fig3": three,
                "convergence": conv,
                "runtime_seconds": manifest["runtime_seconds"],
            },
            indent=2,
        ),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
