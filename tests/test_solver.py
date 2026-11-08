"""Validate the Numerov shooting solver against closed-form results."""
import numpy as np
import pytest

from schrodinger import (
    eigenvalues, find_eigenvalue, grid, harmonic, infinite_well, normalize,
    shoot, uncertainty_product,
)
from schrodinger.analytic import (
    harmonic_eigenvalue, harmonic_uncertainty, well_eigenvalue, well_uncertainty,
)
from scipy.integrate import simpson

X = grid(1000)
N_IDX = np.arange(1, 11)


@pytest.fixture(scope="module")
def well_eps():
    return eigenvalues(10, infinite_well, 200, X)


@pytest.fixture(scope="module")
def harm_eps():
    return eigenvalues(6, harmonic, 1000, X, step=0.05, gap=1 / 200)


def test_well_eigenvalues_match_analytic(well_eps):
    assert np.allclose(well_eps, well_eigenvalue(N_IDX, 200), rtol=1e-7, atol=0)


def test_harmonic_eigenvalues_match_analytic(harm_eps):
    # Low states sit deep below the walls (nu = 1), so the box barely perturbs them.
    assert np.allclose(harm_eps[:5], harmonic_eigenvalue(np.arange(1, 6), 1000), atol=5e-5)


def test_box_walls_shift_higher_harmonic_states(harm_eps):
    # Higher states reach further out, so the hard walls push them up more.
    shift = harm_eps - harmonic_eigenvalue(np.arange(1, 7), 1000)
    assert np.all(np.diff(np.abs(shift)) > 0)


def test_harmonic_levels_equally_spaced(harm_eps):
    spacing = np.diff(harm_eps[:4])
    assert np.allclose(spacing, 2 * np.sqrt(8 / 1000), rtol=1e-4)


def test_eigenfunction_satisfies_right_wall(well_eps):
    psi = normalize(shoot(well_eps[0], infinite_well, 200, X), X)
    assert abs(psi[-1]) < 1e-6


def test_off_eigenvalue_misses_wall(well_eps):
    psi = normalize(shoot(well_eps[0] + 0.1, infinite_well, 200, X), X)
    assert abs(psi[-1]) > 1e-2


def test_nth_state_has_n_minus_1_nodes(well_eps):
    for n, eps in enumerate(well_eps, start=1):
        psi = shoot(eps, infinite_well, 200, X)[1:-2]
        assert np.count_nonzero(np.diff(np.sign(psi))) == n - 1


def test_normalisation():
    psi = normalize(shoot(-0.95, infinite_well, 200, X), X)
    assert simpson(psi**2, x=X) == pytest.approx(1.0, abs=1e-12)


def test_well_uncertainty_products(well_eps):
    u = [uncertainty_product(shoot(e, infinite_well, 200, X), X) for e in well_eps]
    assert np.allclose(u, well_uncertainty(N_IDX), rtol=1e-4)


def test_harmonic_uncertainty_products(harm_eps):
    u = [uncertainty_product(shoot(e, harmonic, 1000, X), X) for e in harm_eps[:4]]
    assert np.allclose(u, harmonic_uncertainty(np.arange(1, 5)), rtol=1e-3)


def test_heisenberg_bound_holds(well_eps, harm_eps):
    for eps, pot, g2 in [(well_eps[0], infinite_well, 200), (harm_eps[0], harmonic, 1000)]:
        assert uncertainty_product(shoot(eps, pot, g2, X), X) >= 0.5 - 1e-4


def test_numerov_is_fourth_order():
    """Eigenvalue error must fall as h^4: fitted log-log slope ~ 4."""
    ns = np.array([41, 81, 161])
    exact = well_eigenvalue(3, 200)
    errs = [
        abs(find_eigenvalue(exact - 0.05, infinite_well, 200, grid(n), step=0.01, tol=1e-14) - exact)
        for n in ns
    ]
    slope = np.polyfit(np.log(1 / (ns - 1)), np.log(errs), 1)[0]
    assert 3.7 < slope < 4.3
