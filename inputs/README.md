# Model inputs

Run `graphene.json` first. It is self-contained and needs no downloads:

```bash
ohc-run --config inputs/graphene.json
```

This uses the existing Haldane implementation with `t2=0`: gapped nearest-neighbor
graphene, hopping 2.8 eV, onsite energies +/-0.25 eV, nearest-neighbor distance
1.42 Angstrom. The band gap is 0.5 eV. All output numbers remain in the program's
internal units; see the root README for physical conversion.

Config paths are relative to the JSON file. Consequently `../output/first-run`
means the repository's `output/first-run/` even if invoked from another directory.
Command-line `--output-dir` paths are relative to the working directory instead.

`wannier90.example.json` is a template, not bundled material data. To use it:

1. Place your `wannier90_hr.dat` and `lattice.txt` under `inputs/wannier90/`.
2. Set the correct Fermi energy (`mu`) and mesh (`n_k`) in your copy of the config.
3. Run `ohc-run --config inputs/wannier90.example.json`.

The lattice is a 2x2 or 3x3 whitespace-delimited matrix with primitive lattice
vectors as rows, in the length unit consistent with the model. The hr file
contains hopping matrices, not lattice vectors or full position matrices.
Private material inputs under `inputs/wannier90/` are ignored by Git.

Supported JSON keys are documented in the root README. Unknown keys and invalid
types are rejected. CLI flags override the corresponding JSON values.
