"""Abstract OAM definition (L2).

An OAM variant maps an :class:`~ohk.eigen.EigenBasis` to a complex matrix ``L[n, n']`` in the
Bloch band basis, in units of hbar.  Everything else — assembly of the diagonal, Hermiticity
enforcement, comparison across variants — is handled outside the variant.

Contract
--------
* Return shape is exactly ``(n_bands, n_bands)``.
* The result must be Hermitian up to numerical tolerance.  A variant may return a matrix that
  violates this (e.g. a deliberately naive one); :func:`ohk.oam.assemble.full_matrix` will catch
  it and report rather than silently symmetrize.
* Under ``U -> U D``, the result transforms as ``L -> D^dag L D``. Diagonal
  expectations and Kubo contractions are invariant. The historical flag
  :attr:`gauge_invariant` refers to those physical observables, not entrywise matrix equality.

Why a class and not a function
------------------------------
Variants carry metadata (``name``, ``gauge_invariant``, the formula they implement, their
reference) that flows into reports and test failure messages.  That metadata is what makes a
side-by-side comparison of conflicting formulas interpretable.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from ohk.eigen import EigenBasis


class OAMVariant(ABC):
    """Base class for an OAM matrix-element definition.

    Attributes
    ----------
    name : str
        Short identifier used in reports and test ids.
    reference : str
        Citation for the definition.
    gauge_invariant : bool
        Historical name: whether the variant has gauge-invariant observables and a
        gauge-covariant matrix under independent phase rotations of the Bloch states.
    diagonal : bool
        Whether the variant supplies the band-diagonal elements.  ``False`` means the diagonal is
        filled in by :func:`ohk.oam.assemble.full_matrix` from the closed form.
    """

    name: str = ""
    reference: str = ""
    gauge_invariant: bool = True
    diagonal: bool = False

    def __init_subclass__(cls, **kwargs) -> None:
        """Guard against shadowing the declaration flags with methods of the same name.

        ``diagonal`` and ``gauge_invariant`` are *data* consumed by the assembler and by the
        validator.  A subclass that also defines a method called ``diagonal`` silently replaces
        the flag with a function object, which is truthy and therefore makes the assembler skip
        filling the diagonal.  That failure is invisible until the numbers are wrong, so it is
        checked at class creation time instead.
        """
        super().__init_subclass__(**kwargs)
        for flag in ("diagonal", "gauge_invariant"):
            value = getattr(cls, flag, None)
            if not isinstance(value, bool):
                raise TypeError(
                    f"{cls.__name__}.{flag} must be a bool, got {type(value).__name__}; "
                    f"rename the colliding method (see ohk.oam.base.OAMVariant docstring)"
                )

    @abstractmethod
    def matrix(self, basis: EigenBasis) -> np.ndarray:
        """Return ``L[n, n']`` as an ``(n_bands, n_bands)`` complex array (units of hbar)."""

    # ------------------------------------------------------------------ helpers
    def __repr__(self) -> str:
        return f"<{type(self).__name__} name={self.name!r} gauge_invariant={self.gauge_invariant}>"

    def __call__(self, basis: EigenBasis) -> np.ndarray:
        return self.matrix(basis)
