"""User-facing configuration, output formats and reproducibility boundaries."""

import builtins
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from ohk import BZGrid, CovariantProjector, HaldaneModel, orbital_hall_conductivity
from ohk.cli import main, main_scan


def config_file(tmp_path, **overrides):
    path = tmp_path / "inputs" / "model.json"
    path.parent.mkdir(exist_ok=True)
    data = dict(
        model="haldane",
        t1=1.0,
        t2=0.0,
        mass=0.4,
        phi=0.0,
        n_k=6,
        mu=0.0,
        fermi_range=[-0.2, 0.2, 5],
        output_dir="../results",
        formats=["csv"],
    )
    data.update(overrides)
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_config_paths_overrides_csv_and_replay(tmp_path, monkeypatch, capsys):
    config = config_file(tmp_path)
    monkeypatch.chdir(tmp_path.parent)
    original = config.read_bytes()
    assert main(["--config", str(config), "--n-k", "9"]) == 0
    folder = tmp_path / "results"
    report = json.loads((folder / "report.json").read_text())
    assert report["metadata"]["n_k"] == 9
    assert report["metadata"]["regularization_eta"] == 0
    assert config.read_bytes() == original
    with (folder / "conductivity.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 5
    assert [float(row["mu"]) for row in rows] == pytest.approx(np.linspace(-0.2, 0.2, 5))
    expected = orbital_hall_conductivity(
        HaldaneModel(t2=0, M=0.4, phi=0),
        CovariantProjector(),
        BZGrid.uniform(HaldaneModel().lattice(), 9),
    ).sigma
    assert float(rows[2]["sigma_oh_internal"]) == pytest.approx(expected)
    paths = report["metadata"]["outputs"]
    assert all(Path(p).is_file() for p in paths)
    assert str((folder / "conductivity.csv").resolve()) in capsys.readouterr().out
    assert not list(folder.glob("*.pdf"))
    before = (folder / "conductivity.csv").read_bytes()
    assert main(["--config", str(folder / "input.json")]) == 0
    assert (folder / "conductivity.csv").read_bytes() == before


def test_csv_only_never_imports_matplotlib(tmp_path, monkeypatch):
    original = builtins.__import__

    def guard(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise AssertionError("CSV-only export imported optional plotting dependency")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guard)
    assert main(["--n-k", "3", "--formats", "csv", "--output-dir", str(tmp_path)]) == 0


def test_missing_plot_dependency_has_actionable_error(tmp_path, monkeypatch, capsys):
    original = builtins.__import__

    def guard(name, *args, **kwargs):
        if name.startswith("matplotlib"):
            raise ImportError("optional dependency unavailable")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guard)
    assert main(["--n-k", "3", "--formats", "pdf", "--output-dir", str(tmp_path / "out")]) == 2
    assert "--formats csv" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize(
    "change",
    [
        {"unknown_option": 1},
        {"n_k": True},
        {"mu": "zero"},
        {"eta": 0},
        {"formats": ["svg"]},
        {"formats": "csv"},
        {"fermi_range": [0, 1, 2.5]},
        {"fermi_range": [1, 0, 3]},
        {"fermi_range": [0, 1]},
        {"convergence": [0, 3]},
        {"nearest_neighbour_distance": -1},
        {"model": "invalid"},
        {"variant": "invalid"},
    ],
)
def test_bad_config_fails_without_output(tmp_path, change):
    path = config_file(tmp_path, **change)
    assert main(["--config", str(path)]) == 2
    assert not (tmp_path / "results").exists()


def test_explicit_fermi_replaces_config_range(tmp_path):
    path = config_file(tmp_path)
    assert main(["--config", str(path), "--fermi", "0", "1"]) == 0
    data = json.loads((tmp_path / "results/report.json").read_text())
    assert data["quantities"]["fermi_energies"] == [0, 1]


def test_pdf_png_and_convergence_exports(tmp_path):
    pytest.importorskip("matplotlib")
    folder = tmp_path / "nested/output"
    assert (
        main(
            [
                "--n-k",
                "6",
                "--fermi-range",
                "-1",
                "1",
                "7",
                "--convergence",
                "3",
                "6",
                "--formats",
                "csv",
                "pdf",
                "png",
                "--output-dir",
                str(folder),
            ]
        )
        == 0
    )
    for stem in ("conductivity", "convergence"):
        assert (folder / f"{stem}.pdf").read_bytes().startswith(b"%PDF-")
        assert (folder / f"{stem}.png").read_bytes().startswith(b"\x89PNG")
        assert (folder / f"{stem}.csv").is_file()
    data = json.loads((folder / "report.json").read_text())
    assert len(data["quantities"]["fermi_sigma"]) == 7
    assert len(data["metadata"]["outputs"]) == 8


def test_pdf_only_preserves_source_data_in_report(tmp_path):
    pytest.importorskip("matplotlib")
    assert (
        main_scan(["--points-per-segment", "4", "--formats", "pdf", "--output-dir", str(tmp_path)])
        == 0
    )
    assert (tmp_path / "oam.pdf").read_bytes().startswith(b"%PDF-")
    assert not (tmp_path / "oam.csv").exists()
    report = json.loads((tmp_path / "report.json").read_text())
    assert len(report["quantities"]["rows"]) == 9
    assert "L_0_1_imag" in report["quantities"]["columns"]


def test_legacy_nested_output(tmp_path):
    path = tmp_path / "nested/legacy.json"
    assert main(["--n-k", "3", "--output", str(path)]) == 0
    assert json.loads(path.read_text())["quantities"]["sigma_oh"] < 0
