"""Dimensionless potentials nu(x) on x in [0, 1] (energies scaled so the well floor is -1)."""
import numpy as np


def infinite_well(x: np.ndarray) -> np.ndarray:
    """Flat floor nu = -1 between hard walls at x = 0 and x = 1."""
    return -np.ones_like(x, dtype=float)


def harmonic(x: np.ndarray) -> np.ndarray:
    """Parabolic well nu = 8 (x - 1/2)^2 - 1, centred in the box."""
    return 8.0 * (x - 0.5) ** 2 - 1.0
