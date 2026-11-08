"""Closed-form results used to validate the solver."""
import numpy as np


def well_eigenvalue(n, gamma2: float):
    """Infinite well: k = n pi  =>  eps_n = n^2 pi^2 / gamma^2 - 1."""
    n = np.asarray(n, dtype=float)
    return n**2 * np.pi**2 / gamma2 - 1.0


def harmonic_eigenvalue(n, gamma2: float):
    """psi'' + gamma^2 (E - 8 y^2) psi = 0 has gamma^2 E_n = (2n - 1) sqrt(8 gamma^2).

    So eps_n = (2n - 1) sqrt(8) / gamma - 1 (n = 1, 2, ...); valid while the
    state stays far below the box walls, where nu = 1.
    """
    n = np.asarray(n, dtype=float)
    return (2 * n - 1) * np.sqrt(8.0 / gamma2) - 1.0


def well_uncertainty(n):
    """Infinite well on [0, 1]: Delta x Delta p = sqrt(n^2 pi^2 / 12 - 1/2)."""
    n = np.asarray(n, dtype=float)
    return np.sqrt(n**2 * np.pi**2 / 12.0 - 0.5)


def harmonic_uncertainty(n):
    """Harmonic oscillator: Delta x Delta p = n - 1/2 (n = 1, 2, ...)."""
    return np.asarray(n, dtype=float) - 0.5
