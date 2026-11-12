"""Run the bundled gapped-graphene input and print where to find the figures/data.

From the repository root: python examples/first_run.py
Optional overrides: python examples/first_run.py --formats csv --n-k 24
The package must first be installed: python -m pip install -e ".[plot]"
"""

from pathlib import Path
import sys

from ohk.cli import main


if __name__ == "__main__":
    config = Path(__file__).resolve().parents[1] / "inputs" / "graphene.json"
    raise SystemExit(main(["--config", str(config), *sys.argv[1:]]))
