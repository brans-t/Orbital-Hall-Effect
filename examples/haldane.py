"""Run a small Haldane benchmark, Fermi scan and convergence study."""

import argparse
from pathlib import Path
import numpy as np
from ohk import HaldaneModel, CovariantProjector, BZGrid, orbital_hall_conductivity, fermi_scan
from ohk.transport import convergence_study
from ohk.report import Report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-k", type=int, default=30)
    parser.add_argument("--output", type=Path, default=Path("output/haldane.json"))
    args = parser.parse_args()
    model, variant = HaldaneModel(phi=np.pi / 2), CovariantProjector()
    grid = BZGrid.uniform(model.lattice(), args.n_k)
    result = orbital_hall_conductivity(model, variant, grid)
    mus = np.linspace(-3, 3, 61)
    report = Report("Haldane OHC at phi=pi/2", metadata={"n_k": args.n_k, "units": result.units})
    report.add("sigma_oh", result.sigma).add("fermi_energies", mus)
    report.add("fermi_sigma", fermi_scan(model, variant, grid, mus))
    report.add("convergence", convergence_study(model, variant, [12, 24, 36]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report.to_json(args.output)
    print(report.summary())


if __name__ == "__main__":
    main()
