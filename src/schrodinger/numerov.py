"""Numerov integration of psi'' + k^2(x) psi = 0, with k^2 = gamma^2 (eps - nu(x)).

Numerov's three-term recurrence is O(h^6) locally and O(h^4) globally:

    (1 + h^2 k_{i+1}^2 / 12) psi_{i+1}
        = 2 (1 - 5 h^2 k_i^2 / 12) psi_i - (1 + h^2 k_{i-1}^2 / 12) psi_{i-1}
"""
from collections.abc import Callable

import numpy as np

Potential = Callable[[np.ndarray], np.ndarray]


def grid(n_points: int = 1000, length: float = 1.0) -> np.ndarray:
    """Uniform grid including both walls."""
    return np.linspace(0.0, length, n_points)


def shoot(
    eps: float,
    potential: Potential,
    gamma2: float,
    x: np.ndarray,
    psi0: float = 0.0,
    psi1: float = 1e-4,
) -> np.ndarray:
    """Integrate from the left wall; psi(x_end) = 0 iff eps is an eigenvalue.

    psi1 only fixes the (arbitrary) amplitude: the equation is linear.
    """
    h = x[1] - x[0]
    k2 = gamma2 * (eps - potential(x))
    f = 1.0 + h * h * k2 / 12.0  # Numerov weights
    psi = np.empty_like(x, dtype=float)
    psi[0], psi[1] = psi0, psi1
    # Plain floats: the recurrence is inherently sequential, and scalar maths
    # is several times faster than indexing numpy arrays in the loop.
    fl = f.tolist()
    p_prev, p_cur = psi0, psi1
    for i in range(1, len(x) - 1):
        p_next = ((12.0 - 10.0 * fl[i]) * p_cur - fl[i - 1] * p_prev) / fl[i + 1]
        psi[i + 1] = p_next
        p_prev, p_cur = p_cur, p_next
    return psi
