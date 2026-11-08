"""Eigenvalue search: march eps upward until psi(x_end) changes sign, then bisect.

Each sign change of the end-point value brackets one eigenvalue; reversing and
halving the step on every crossing converges to it geometrically.
"""
import numpy as np

from schrodinger.numerov import Potential, shoot


def _endpoint(eps: float, potential: Potential, gamma2: float, x: np.ndarray) -> float:
    return shoot(eps, potential, gamma2, x)[-1]


def find_eigenvalue(
    eps_start: float,
    potential: Potential,
    gamma2: float,
    x: np.ndarray,
    step: float = 0.04,
    tol: float = 1e-10,
) -> float:
    """First eigenvalue above eps_start (step must be smaller than the level spacing)."""
    eps = eps_start
    f_prev = _endpoint(eps, potential, gamma2, x)
    while abs(step) >= tol:
        eps += step
        f_cur = _endpoint(eps, potential, gamma2, x)
        if f_prev * f_cur < 0:
            step /= -2.0
        f_prev = f_cur
    return eps


def eigenvalues(
    n: int,
    potential: Potential,
    gamma2: float,
    x: np.ndarray,
    eps_start: float = -1.0,
    step: float = 0.04,
    tol: float = 1e-10,
    gap: float | None = None,
) -> np.ndarray:
    """Lowest n eigenvalues, each search starting just above the previous root."""
    gap = 1.0 / gamma2 if gap is None else gap
    out = []
    eps = eps_start
    for _ in range(n):
        e = find_eigenvalue(eps, potential, gamma2, x, step=step, tol=tol)
        out.append(e)
        eps = e + gap
    return np.array(out)
