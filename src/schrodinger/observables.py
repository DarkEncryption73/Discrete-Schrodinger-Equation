"""Normalisation and position/momentum spreads (hbar = 1, p = -i d/dx)."""
import numpy as np
from scipy.integrate import simpson


def normalize(psi: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Scale psi so that the integral of psi^2 is 1 (Simpson's rule)."""
    return psi / np.sqrt(simpson(psi**2, x=x))


def delta_x(psi: np.ndarray, x: np.ndarray) -> float:
    """sqrt(<x^2> - <x>^2) for a normalised real psi."""
    mean = simpson(x * psi**2, x=x)
    mean2 = simpson(x**2 * psi**2, x=x)
    return float(np.sqrt(mean2 - mean**2))


def delta_p(psi: np.ndarray, x: np.ndarray) -> float:
    """sqrt(<p^2>), since <p> = 0 for a real bound state; <p^2> = -<psi|psi''>."""
    h = x[1] - x[0]
    d2psi = (psi[:-2] - 2.0 * psi[1:-1] + psi[2:]) / h**2
    return float(np.sqrt(-simpson(psi[1:-1] * d2psi, x=x[1:-1])))


def uncertainty_product(psi: np.ndarray, x: np.ndarray) -> float:
    """Delta x * Delta p of the normalised state."""
    psi = normalize(psi, x)
    return delta_x(psi, x) * delta_p(psi, x)
