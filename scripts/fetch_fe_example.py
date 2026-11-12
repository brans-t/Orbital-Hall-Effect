"""Fetch a pinned public WannierTools Fe dataset into .local, verifying both checksums.

The third-party dataset remains outside version control. No archive paths are extracted.
Run: python scripts/fetch_fe_example.py
Source: https://github.com/quanshengwu/wannier_tools/tree/362171ce35450118828c6960a882e3a6d08c73f5/examples/Fe
"""

from pathlib import Path
import hashlib
import io
import json
import tarfile
from urllib.request import urlopen
import numpy as np

COMMIT = "362171ce35450118828c6960a882e3a6d08c73f5"
BASE = f"https://raw.githubusercontent.com/quanshengwu/wannier_tools/{COMMIT}/examples/Fe"
ARCHIVE_SHA = "388792dc706600f1b63e57e1c0a8ba06eb32b60d43e24a3c7138d548faee269f"
HR_SHA = "80bc1e558e0095532613e2b829ccded2b3e00652abd1f2da15cf47588a2bcd12"


def main():
    root = Path(".local/fe")
    root.mkdir(parents=True, exist_ok=True)
    archive = root / "wannier90_hr.dat.tar.gz"
    if not archive.exists():
        with urlopen(BASE + "/wannier90_hr.dat.tar.gz", timeout=60) as response:
            archive.write_bytes(response.read())
    data = archive.read_bytes()
    if hashlib.sha256(data).hexdigest() != ARCHIVE_SHA:
        raise ValueError("Fe archive checksum mismatch")
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and m.name == "wannier90_hr.dat"]
        if len(members) != 1:
            raise ValueError("Archive must contain exactly the expected hr.dat member")
        hr = tar.extractfile(members[0]).read()
    if hashlib.sha256(hr).hexdigest() != HR_SHA:
        raise ValueError("Fe hr.dat checksum mismatch")
    (root / "wannier90_hr.dat").write_bytes(hr)
    # Lattice and Fermi level from the pinned wt.in-ahe.
    np.savetxt(
        root / "lattice.txt",
        [[1.435, 1.435, -1.435], [-1.435, 1.435, 1.435], [1.435, -1.435, 1.435]],
    )
    provenance = {
        "source": BASE,
        "commit": COMMIT,
        "archive_sha256": ARCHIVE_SHA,
        "hr_sha256": HR_SHA,
        "fermi_energy_eV": 9.1839,
        "lattice_source": BASE + "/wt.in-ahe",
        "note": "External validation data; not redistributed as part of ohk",
    }
    (root / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Fe data ready at {root}; use mu=9.1839 and lattice.txt")


if __name__ == "__main__":
    main()
