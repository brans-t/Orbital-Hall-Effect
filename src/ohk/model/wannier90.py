"""Native Wannier90 hr reader/writer, with degeneracy weights and analytic derivatives.

Uses the periodic cell gauge and Cartesian k. The hr file alone contains no position-operator
matrix elements: results are those of the Hamiltonian-only, diagonal-position approximation.
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
from ohk.model.base import TightBindingModel, ModelSpec


class Wannier90Model(TightBindingModel):
    def __init__(self, spec, hr, r_vectors, degeneracies=None):
        self.spec = spec
        self.hr = np.array(hr, dtype=complex, copy=True)
        self.r_vectors = np.array(r_vectors, dtype=float, copy=True)
        if self.hr.ndim != 3 or self.hr.shape[1:] != (spec.num_bands,) * 2:
            raise ValueError("hr must have shape (nrpts, num_bands, num_bands)")
        n = len(self.hr)
        if n == 0 or self.r_vectors.shape != (n, spec.periodic_dim):
            raise ValueError("r_vectors must have shape (nrpts, periodic_dim)")
        self.degeneracies = np.ones(n, int) if degeneracies is None else np.asarray(degeneracies)
        if self.degeneracies.shape != (n,) or np.any(self.degeneracies <= 0):
            raise ValueError("degeneracies must contain nrpts positive integers")
        if not np.isfinite(self.degeneracies).all() or np.any(
            self.degeneracies != np.rint(self.degeneracies)
        ):
            raise ValueError("degeneracies must be positive integers")
        if not np.isfinite(self.hr).all() or not np.isfinite(self.r_vectors).all():
            raise ValueError("Hamiltonian and R vectors must be finite")
        self.degeneracies = self.degeneracies.astype(int)
        self._weighted = self.hr / self.degeneracies[:, None, None]

    @classmethod
    def from_hr_dat(cls, path, lattice_vectors=None):
        if lattice_vectors is None:
            raise ValueError("lattice_vectors are required; hr.dat does not contain them")
        path = Path(path)
        try:
            with path.open(encoding="utf-8-sig") as stream:
                comment = stream.readline().strip()
                n_bands, nrpts = int(stream.readline()), int(stream.readline())
                if n_bands <= 0 or nrpts <= 0:
                    raise ValueError("num_wann and nrpts must be positive")
                degeneracies = []
                while len(degeneracies) < nrpts:
                    line = stream.readline()
                    if not line:
                        raise ValueError("truncated degeneracy list")
                    degeneracies.extend(map(int, line.split()))
                if len(degeneracies) != nrpts:
                    raise ValueError("incorrect degeneracy count")
                blocks = np.zeros((nrpts, n_bands, n_bands), complex)
                cells = []
                for block in range(nrpts):
                    seen, r_current = set(), None
                    for _ in range(n_bands**2):
                        fields = stream.readline().split()
                        if len(fields) != 7:
                            raise ValueError("expected seven fields per hr matrix element")
                        r = tuple(map(int, fields[:3]))
                        i, j = int(fields[3]) - 1, int(fields[4]) - 1
                        if not 0 <= i < n_bands or not 0 <= j < n_bands or (i, j) in seen:
                            raise ValueError("invalid or duplicate orbital indices")
                        if r_current is not None and r != r_current:
                            raise ValueError("inconsistent R within a block")
                        seen.add((i, j))
                        r_current = r
                        blocks[block, i, j] = complex(
                            *[float(x.replace("D", "E").replace("d", "e")) for x in fields[5:]]
                        )
                    cells.append(r_current)
                if any(line.strip() for line in stream):
                    raise ValueError("unexpected trailing matrix elements")
        except (IndexError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid Wannier90 hr file {path}: {exc}") from exc
        a = np.asarray(lattice_vectors, dtype=float)
        spec = ModelSpec(
            n_bands,
            a,
            periodic_dim=len(a),
            name=path.stem,
            metadata={
                "source": str(path),
                "comment": comment,
                "position_approximation": "Hamiltonian-only cell gauge",
            },
        )
        cells = np.asarray(cells)
        if spec.periodic_dim == 2 and np.any(cells[:, 2] != 0):
            raise ValueError("2D lattice cannot represent nonzero third R components")
        if len(set(map(tuple, cells))) != nrpts:
            raise ValueError("duplicate R blocks")
        return cls(spec, blocks, cells[:, : spec.periodic_dim] @ a, degeneracies)

    @classmethod
    def from_pythtb(cls, w90):
        """Adapt a legacy PythTB tb_model (spinless and fully periodic).

        A w90 object is also accepted if its model() method needs no extra arguments.
        """
        obj = w90.model() if callable(getattr(w90, "model", None)) else w90
        required = ("_lat", "_orb", "_hoppings", "_site_energies", "_nspin", "_dim_k")
        if not all(hasattr(obj, name) for name in required):
            raise TypeError("Expected a legacy PythTB tb_model or w90 object")
        if obj._nspin != 1 or obj._dim_k != len(obj._lat):
            raise ValueError("Only spinless, fully periodic PythTB models are supported")
        from ohk.model.analytic import AnalyticModel, Hopping

        spec = ModelSpec(
            len(obj._orb),
            np.asarray(obj._lat),
            periodic_dim=obj._dim_k,
            name="PythTB",
            metadata={"source": "PythTB legacy adapter"},
        )
        # Return the common model contract while retaining physical orbital embedding.
        return AnalyticModel(
            spec,
            [
                Hopping(int(i), int(j), tuple(map(int, r)), complex(t))
                for t, i, j, r in obj._hoppings
            ],
            onsite=np.asarray(obj._site_energies),
            orbital_positions=np.asarray(obj._orb) @ spec.lattice_vectors,
        )

    def H(self, k):
        phase = np.exp(1j * (self.r_vectors @ self._momentum(k)))
        return np.einsum("r,rij->ij", phase, self._weighted)

    def dH(self, k):
        phase = np.exp(1j * (self.r_vectors @ self._momentum(k)))
        return np.einsum("r,rd,rij->dij", 1j * phase, self.r_vectors, self._weighted)

    def to_hr_dat(self, path, comment=""):
        cells = self.r_vectors @ np.linalg.inv(self.lattice())
        if not np.allclose(cells, np.rint(cells), atol=1e-9, rtol=0):
            raise ValueError("R vectors are not integer lattice translations")
        cells = np.rint(cells).astype(int)
        cells = np.pad(cells, ((0, 0), (0, 3 - self.spec.periodic_dim)))
        if "\n" in comment or "\r" in comment:
            raise ValueError("comment must be a single line")
        lines = [comment or "Written by ohk", str(self.num_bands), str(len(self.hr))]
        for start in range(0, len(self.hr), 15):
            lines.append(" ".join(str(x) for x in self.degeneracies[start : start + 15]))
        for r, block in zip(cells, self.hr):
            for j in range(self.num_bands):
                for i in range(self.num_bands):
                    t = block[i, j]
                    lines.append(
                        f"{r[0]} {r[1]} {r[2]} {i + 1} {j + 1} {t.real:.17e} {t.imag:.17e}"
                    )
        Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
