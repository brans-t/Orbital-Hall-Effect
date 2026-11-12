"""Symmetrized orbital current J_a = {v_a,L_z}/2."""

import numpy as np


def orbital_current(L, v, direction=0):
    L, v = np.asarray(L), np.asarray(v)
    if L.ndim != 2 or L.shape[0] != L.shape[1] or v.ndim != 3 or v.shape[1:] != L.shape:
        raise ValueError("L and v must have matching band dimensions")
    if not 0 <= direction < len(v):
        raise ValueError("invalid current direction")
    return 0.5 * (v[direction] @ L + L @ v[direction])


def orbital_current_tensor(L, v):
    return np.stack([orbital_current(L, v, j) for j in range(len(v))])


def check_hermitian_current(J, tol=1e-10):
    J = np.asarray(J)
    return float(np.linalg.norm(J - J.conj().swapaxes(-1, -2)))
