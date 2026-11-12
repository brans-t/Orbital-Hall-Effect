"""Thin, reproducible command-line wrappers around L0-L4."""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
from ohk.eigen import EigenBasis
from ohk.model.haldane import HaldaneModel
from ohk.oam.assemble import full_matrix
from ohk.report import Report
from ohk.transport.grids import BZGrid
from ohk.transport.kubo import orbital_hall_conductivity, fermi_scan, convergence_study


def _model_options(parser):
    parser.add_argument("--model", choices=["haldane", "wannier90"], default="haldane")
    parser.add_argument("--t1", type=float, default=1.0)
    parser.add_argument("--t2", type=float, default=0.15)
    parser.add_argument("--mass", "--M", type=float, default=0.4)
    parser.add_argument("--phi", type=float, default=np.pi / 4, help="Hopping phase in radians")
    parser.add_argument("--hr", type=Path, help="Wannier90 hr.dat; implies --model wannier90")
    parser.add_argument(
        "--lattice", type=Path, help="Whitespace-delimited lattice matrix, rows in model lengths"
    )
    parser.add_argument(
        "--variant", choices=["covariant", "spectral", "pezo", "aca"], default="covariant"
    )


def _build_parser_run():
    p = argparse.ArgumentParser(prog="ohc-run", description="Compute orbital Hall conductivity")
    _model_options(p)
    p.add_argument("--n-k", type=int, default=40)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--eta", type=float, help="Opt-in broadening for OAM and Kubo denominators")
    p.add_argument("--fermi", nargs="+", type=float, help="Additional chemical potentials to scan")
    p.add_argument("--convergence", nargs="+", type=int, help="Additional mesh sizes")
    p.add_argument("--output", type=Path, help="JSON report path")
    return p


def _build_parser_scan():
    p = argparse.ArgumentParser(
        prog="oam-scan", description="Export OAM along a reciprocal-space path"
    )
    _model_options(p)
    p.add_argument("--points-per-segment", type=int, default=50)
    p.add_argument(
        "--vertices", help="Semicolon-separated fractional vertices, e.g. '0,0;0.333,0.667'"
    )
    p.add_argument("--output", type=Path, default=Path("oam-scan.csv"))
    return p


def _build_parser_validate():
    p = argparse.ArgumentParser(
        prog="ohk-validate", description="Run independent consistency checks"
    )
    p.add_argument("--n-k", type=int, default=30)
    p.add_argument("--output", type=Path, help="JSON validation report")
    return p


def _model_variant(args):
    if args.hr is not None or args.model == "wannier90":
        if args.hr is None or args.lattice is None:
            raise ValueError("Wannier90 requires --hr and --lattice")
        from ohk.model.wannier90 import Wannier90Model

        model = Wannier90Model.from_hr_dat(args.hr, np.loadtxt(args.lattice))
    else:
        model = HaldaneModel(t1=args.t1, t2=args.t2, M=args.mass, phi=args.phi)
    from ohk.oam import CovariantProjector, SpectralSum, PezoTruncated
    from ohk.oam.atomic import AtomCentered, ACASpec

    if args.variant == "aca":
        if not isinstance(model, HaldaneModel):
            raise ValueError(
                "CLI ACA supports the known s-orbital Haldane model; use Python ACASpec for other models"
            )
        variant = AtomCentered(ACASpec((0, 0)))
    else:
        variant = {"covariant": CovariantProjector, "spectral": SpectralSum, "pezo": PezoTruncated}[
            args.variant
        ]()
    return model, variant


def main(argv=None):
    args = _build_parser_run().parse_args(argv)
    try:
        model, variant = _model_variant(args)
        grid = BZGrid.uniform(model.lattice(), args.n_k)
        kwargs = {} if args.eta is None else {"subband_regularization": True, "eta": args.eta}
        result = orbital_hall_conductivity(model, variant, grid, mu=args.mu, **kwargs)
        report = Report(
            "Orbital Hall conductivity",
            metadata={
                "model": model.spec.name,
                "parameters": model.spec.metadata,
                "variant": variant.name,
                "n_k": args.n_k,
                "mu": args.mu,
                "units": result.units,
                "regularization_eta": result.regularization_eta,
            },
        )
        report.add("sigma_oh", result.sigma).add("min_gap", result.min_gap)
        report.add("degenerate_fraction", result.degenerate_fraction)
        if result.degenerate_fraction > 0.001:
            report.warn("Near-degenerate points detected; inspect mesh convergence")
        if args.fermi is not None:
            report.add("fermi_energies", args.fermi)
            report.add("fermi_sigma", fermi_scan(model, variant, grid, args.fermi, **kwargs))
        if args.convergence is not None:
            report.add(
                "convergence",
                convergence_study(model, variant, args.convergence, mu=args.mu, **kwargs),
            )
        if args.output:
            report.to_json(args.output)
        print(report.summary())
        return 0
    except (ValueError, TypeError, OSError) as exc:
        print(f"ohc-run: {exc}", file=sys.stderr)
        return 2


def main_scan(argv=None):
    args = _build_parser_scan().parse_args(argv)
    try:
        model, variant = _model_variant(args)
        if args.vertices:
            vertices = [
                [float(x) for x in vertex.split(",")] for vertex in args.vertices.split(";")
            ]
            labels = [str(i) for i in range(len(vertices))]
        elif model.spec.periodic_dim == 2:
            vertices, labels = [[0, 0], [1 / 3, 2 / 3], [2 / 3, 1 / 3]], ["Gamma", "K", "Kprime"]
        else:
            vertices, labels = [[0, 0, 0], [0.5, 0, 0]], ["Gamma", "X"]
        grid = BZGrid.path(model.lattice(), vertices, labels, args.points_per_segment)
        columns = ["distance"] + [f"k{j}" for j in range(model.spec.periodic_dim)]
        columns += [f"eps_{n}" for n in range(model.num_bands)]
        columns += [
            f"L_{n}_{m}_{part}"
            for n in range(model.num_bands)
            for m in range(model.num_bands)
            for part in ("real", "imag")
        ]
        rows, distance, previous = [], 0.0, None
        for k in grid.kpoints:
            if previous is not None:
                distance += np.linalg.norm(k - previous)
            previous = k
            basis = EigenBasis.at(model, k)
            L = full_matrix(variant, basis)[0]
            components = np.stack([L.real, L.imag], axis=-1).ravel().tolist()
            rows.append([distance, *k, *basis.eps, *components])
        with args.output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerows(rows)
        print(f"Wrote {len(rows)} momentum points to {args.output}")
        return 0
    except (ValueError, TypeError, OSError) as exc:
        print(f"oam-scan: {exc}", file=sys.stderr)
        return 2


def main_validate(argv=None):
    args = _build_parser_validate().parse_args(argv)
    try:
        from ohk.validation.benchmarks import run_haldane_benchmark

        measured = run_haldane_benchmark(n_k=args.n_k)
        report = Report("Orbital Hall kit validation", checks=measured.pop("checks"))
        report.quantities.update(measured)
        if args.output:
            report.to_json(args.output)
        print(report.summary())
        return 0 if measured["passed"] else 1
    except (ValueError, TypeError, OSError) as exc:
        print(f"ohk-validate: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
