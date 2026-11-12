"""End-to-end commands and complex JSON report round trips."""

import csv
import json
import numpy as np
from ohk.cli import main, main_scan, main_validate
from ohk.report import Report


def test_report_serializes_complex_without_nan(tmp_path):
    report = Report("test", metadata={"z": np.complex128(1 + 2j)})
    assert report.add("matrix", np.array([[1, 2j]])) is report
    report.add("nonfinite", np.inf)
    report.checks.update({"okay": 1e-12, "bad": float("nan")})
    report.warn("example")
    path = tmp_path / "report.json"
    report.to_json(path)
    value = json.loads(path.read_text(encoding="utf-8"))
    assert value["quantities"]["matrix"][0][1] == {"real": 0, "imag": 2}
    assert value["quantities"]["nonfinite"] is None
    assert list(report.failed_checks()) == ["bad"]
    assert "FAIL" in report.summary()


def test_cli_conductivity_and_fermi_scan(tmp_path):
    output = tmp_path / "run.json"
    assert (
        main(
            [
                "--n-k",
                "6",
                "--fermi",
                "-1",
                "0",
                "1",
                "--convergence",
                "3",
                "6",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    value = json.loads(output.read_text(encoding="utf-8"))
    assert value["quantities"]["sigma_oh"] < 0
    assert len(value["quantities"]["fermi_sigma"]) == 3
    assert value["metadata"]["units"] == "e/(2*pi)"


def test_cli_scan_exports_complex_elements(tmp_path):
    output = tmp_path / "scan.csv"
    assert main_scan(["--points-per-segment", "4", "--output", str(output)]) == 0
    with output.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 9
    assert "L_0_1_real" in rows[0] and "L_0_1_imag" in rows[0]
    assert float(rows[-1]["distance"]) > 0


def test_cli_full_validation(tmp_path):
    output = tmp_path / "validation.json"
    assert main_validate(["--n-k", "12", "--output", str(output)]) == 0
    value = json.loads(output.read_text(encoding="utf-8"))
    assert value["quantities"]["passed"] is True
    assert max(value["checks"].values()) < 1e-8


def test_cli_wannier_model(tmp_path):
    from ohk.model.wannier90 import Wannier90Model
    from ohk.model import ModelSpec

    hr = tmp_path / "toy_hr.dat"
    lattice = tmp_path / "lattice.txt"
    model = Wannier90Model(
        ModelSpec(2, np.eye(2)), np.array([np.diag([-1.0, 1.0])]), np.zeros((1, 2))
    )
    model.to_hr_dat(hr)
    np.savetxt(lattice, np.eye(2))
    output = tmp_path / "ohc.json"
    assert (
        main(["--hr", str(hr), "--lattice", str(lattice), "--n-k", "3", "--output", str(output)])
        == 0
    )
    assert json.loads(output.read_text())["quantities"]["sigma_oh"] == 0


def test_cli_errors_have_nonzero_status(capsys):
    assert main(["--n-k", "0"]) == 2
    assert "positive integer" in capsys.readouterr().err
    assert main(["--model", "wannier90"]) == 2
    assert main_scan(["--vertices", "garbage"]) == 2
