"""Thin, reproducible command-line wrappers around L0-L4."""

from __future__ import annotations
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
from ohk.eigen import EigenBasis
from ohk.config import load_run_config, output_formats, resolved_input, validate_run
from ohk.model.haldane import HaldaneModel
from ohk.oam.assemble import full_matrix
from ohk.report import Report
from ohk.transport.grids import BZGrid
from ohk.transport.kubo import orbital_hall_conductivity, fermi_scan, convergence_study
from ohk.outputs import export_conductivity, export_oam_path, require_plotting


def _model_options(parser):
    parser.add_argument("--model", choices=["haldane", "wannier90"], default="haldane")
    parser.add_argument("--t1", type=float, default=1.0)
    parser.add_argument("--t2", type=float, default=0.15)
    parser.add_argument("--mass", "--M", type=float, default=0.4)
    parser.add_argument("--phi", type=float, default=np.pi / 4, help="Hopping phase in radians")
    parser.add_argument("--nearest-neighbour-distance", type=float, default=1.0)
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
    p.add_argument("--config", type=Path, help="JSON config; command-line values override its keys")
    p.add_argument("--n-k", type=int, default=40)
    p.add_argument("--mu", type=float, default=0.0)
    p.add_argument("--eta", type=float, help="Opt-in broadening for OAM and Kubo denominators")
    scan = p.add_mutually_exclusive_group()
    scan.add_argument(
        "--fermi", nargs="+", type=float, help="Additional chemical potentials to scan"
    )
    scan.add_argument(
        "--fermi-range",
        nargs=3,
        type=float,
        metavar=("START", "STOP", "COUNT"),
        help="Inclusive linear Fermi scan",
    )
    p.add_argument("--convergence", nargs="+", type=int, help="Additional mesh sizes")
    p.add_argument("--output-dir", type=Path, help="Result directory (default: output/ohc)")
    p.add_argument(
        "--formats",
        nargs="+",
        choices=["csv", "pdf", "png"],
        help="Exports (default: csv pdf); report.json and input.json are always saved",
    )
    p.add_argument("--output", type=Path, help="Legacy JSON-only report path; prefer --output-dir")
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
    p.add_argument("--output-dir", type=Path, default=Path("output/oam"))
    p.add_argument("--formats", nargs="+", choices=["csv", "pdf", "png"])
    p.add_argument("--output", type=Path, help="Legacy CSV-only path; prefer --output-dir")
    return p


def _build_parser_validate():
    p = argparse.ArgumentParser(
        prog="ohk-validate", description="Run independent consistency checks"
    )
    p.add_argument("--n-k", type=int, default=30)
    p.add_argument(
        "--output",
        type=Path,
        default=Path("output/validation/report.json"),
        help="JSON validation report (default: output/validation/report.json)",
    )
    return p


def _model_variant(args):
    if args.hr is not None or args.model == "wannier90":
        if args.hr is None or args.lattice is None:
            raise ValueError("Wannier90 requires --hr and --lattice")
        from ohk.model.wannier90 import Wannier90Model

        model = Wannier90Model.from_hr_dat(args.hr, np.loadtxt(args.lattice))
    else:
        model = HaldaneModel(
            t1=args.t1,
            t2=args.t2,
            M=args.mass,
            phi=args.phi,
            nearest_neighbour_distance=args.nearest_neighbour_distance,
        )
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
    parser = _build_parser_run()
    try:
        argv = list(sys.argv[1:] if argv is None else argv)
        bootstrap = argparse.ArgumentParser(add_help=False)
        bootstrap.add_argument("--config", type=Path)
        config_path = bootstrap.parse_known_args(argv)[0].config
        defaults = load_run_config(config_path) if config_path is not None else {}
        args = parser.parse_args(argv, namespace=argparse.Namespace(**defaults))
        # Explicit scan flags also override a different scan specification in the config.
        flags = {arg.split("=", 1)[0] for arg in argv}
        if "--fermi" in flags:
            args.fermi_range = None
        elif "--fermi-range" in flags:
            args.fermi = None
        validate_run(args)
        if args.output is not None and (args.formats is not None or args.output_dir is not None):
            raise ValueError("--output is JSON-only; use --output-dir with --formats")
        args.formats = output_formats(args.formats)
        args.output_dir = (args.output_dir or Path("output/ohc")).resolve()
        if args.output is None:
            require_plotting(args.formats)
        if args.fermi_range is not None:
            start, stop, count = args.fermi_range
            mus = np.linspace(start, stop, int(count)).tolist()
        else:
            mus = args.fermi
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
                "config": str(config_path.resolve()) if config_path else None,
                "lattice_vectors": model.lattice(),
                "num_bands": model.num_bands,
                "unit_note": "Internal hbar=m_e=e=1; eV/Angstrom models need UnitSystem conversion",
            },
        )
        report.add("sigma_oh", result.sigma).add("min_gap", result.min_gap)
        report.add("degenerate_fraction", result.degenerate_fraction)
        if result.degenerate_fraction > 0.001:
            report.warn("Near-degenerate points detected; inspect mesh convergence")
        if mus is not None:
            report.add("fermi_energies", np.asarray(mus))
            report.add("fermi_sigma", fermi_scan(model, variant, grid, mus, **kwargs))
        if args.convergence is not None:
            report.add(
                "convergence",
                convergence_study(model, variant, args.convergence, mu=args.mu, **kwargs),
            )
        if args.output:
            report.to_json(args.output)
            files = [args.output]
        else:
            files = export_conductivity(report, args.output_dir, args.formats, resolved_input(args))
        print(report.summary())
        _print_outputs(files)
        return 0
    except (ValueError, TypeError, OSError) as exc:
        print(f"ohc-run: {exc}", file=sys.stderr)
        return 2


def main_scan(argv=None):
    args = _build_parser_scan().parse_args(argv)
    try:
        if args.output is not None and args.formats is not None:
            raise ValueError("--output is CSV-only; use --output-dir with --formats")
        formats = output_formats(args.formats)
        if args.output is None:
            require_plotting(formats)
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
        rows, distance, previous, distances = [], 0.0, None, []
        for k in grid.kpoints:
            if previous is not None:
                distance += np.linalg.norm(k - previous)
            previous = k
            distances.append(distance)
            basis = EigenBasis.at(model, k)
            L = full_matrix(variant, basis)[0]
            components = np.stack([L.real, L.imag], axis=-1).ravel().tolist()
            rows.append([distance, *k, *basis.eps, *components])
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(columns)
                writer.writerows(rows)
            files = [args.output]
        else:
            ticks = [distances[i * args.points_per_segment] for i in range(len(vertices))]
            report = Report(
                "OAM momentum path",
                metadata={
                    "model": model.spec.name,
                    "parameters": model.spec.metadata,
                    "variant": variant.name,
                    "vertices_fractional": vertices,
                    "labels": labels,
                    "lattice_vectors": model.lattice(),
                    "num_bands": model.num_bands,
                    "unit_note": "Energies and OAM are in internal model units",
                },
            )
            report.add("columns", columns).add("rows", rows)
            files = export_oam_path(
                report,
                args.output_dir,
                formats,
                model.num_bands,
                model.spec.periodic_dim,
                ticks,
                labels,
            )
        print(f"Computed {len(rows)} momentum points")
        _print_outputs(files)
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
        if args.output:
            _print_outputs([args.output])
        return 0 if measured["passed"] else 1
    except (ValueError, TypeError, OSError) as exc:
        print(f"ohk-validate: {exc}", file=sys.stderr)
        return 2


def _print_outputs(files):
    print("\nOutput files:")
    for path in files:
        print(f"  {Path(path).resolve()}")


if __name__ == "__main__":
    sys.exit(main())
