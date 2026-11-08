"""Regenerate every figure in results/ from the solver package.

Run from the repo root:  python scripts/make_figures.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from schrodinger import (
    eigenvalues, find_eigenvalue, grid, harmonic, infinite_well, normalize, shoot,
    uncertainty_product,
)
from schrodinger.analytic import (
    harmonic_uncertainty, well_eigenvalue, well_uncertainty,
)

OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)
X = grid(1000)
G_WELL, G_HARM = 200, 1000
HARM = dict(step=0.05, gap=1 / 200)


def wavefunctions():
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    cases = [(infinite_well, G_WELL, {}, r"Infinite well, $\nu(x) = -1$"),
             (harmonic, G_HARM, HARM, r"Harmonic, $\nu(x) = 8(x - 1/2)^2 - 1$")]
    for ax, (pot, g2, kw, title) in zip(axs, cases):
        for n, e in enumerate(eigenvalues(4, pot, g2, X, **kw), start=1):
            ax.plot(X, normalize(shoot(e, pot, g2, X), X), label=rf"$\psi_{n}$")
        ax.set_title(title)
        ax.set_xlabel(r"$x$")
        ax.legend()
    axs[0].set_ylabel(r"$\psi_n(x)$")
    fig.tight_layout()
    fig.savefig(OUT / "wavefunctions.png", dpi=150)


def convergence():
    ns = np.array([21, 41, 81, 161, 321])
    exact = well_eigenvalue(3, G_WELL)
    errs = [abs(find_eigenvalue(exact - 0.05, infinite_well, G_WELL, grid(n), step=0.01, tol=1e-15) - exact)
            for n in ns]
    h = 1 / (ns - 1)
    slope = np.polyfit(np.log(h), np.log(errs), 1)[0]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(h, errs, "ro", label="Numerov error")
    ax.loglog(h, errs[0] * (h / h[0]) ** 4, "k--", label=r"$\propto h^4$")
    ax.set_xlabel(r"grid spacing $h$")
    ax.set_ylabel(r"$|\epsilon_3 - \epsilon_3^{\mathrm{exact}}|$")
    ax.set_title(f"Eigenvalue convergence (fitted order {slope:.2f})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "convergence.png", dpi=150)
    return slope


def uncertainties():
    n = np.arange(1, 11)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    cases = [(infinite_well, G_WELL, {}, well_uncertainty(n), "Infinite well"),
             (harmonic, G_HARM, HARM, harmonic_uncertainty(n), "Harmonic")]
    for ax, (pot, g2, kw, exact, title) in zip(axs, cases):
        num = [uncertainty_product(shoot(e, pot, g2, X), X) for e in eigenvalues(10, pot, g2, X, **kw)]
        ax.plot(n, exact, "k_", ms=14, label="analytic")
        ax.plot(n, num, "ro", label="Numerov")
        ax.axhline(0.5, color="gray", ls=":")
        ax.set_title(title)
        ax.set_xlabel(r"state index $n$")
        ax.set_xticks(n)
        ax.legend()
    axs[0].set_ylabel(r"$\Delta x\,\Delta p$")
    fig.tight_layout()
    fig.savefig(OUT / "uncertainty.png", dpi=150)


def level_spacing():
    eps = eigenvalues(101, harmonic, G_HARM, X, **HARM)
    n = np.arange(1, 101)
    d = np.diff(eps)
    mask = n >= 10
    slope = np.polyfit(np.log(n[mask]), np.log(d[mask]), 1)[0]
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(n, d, "ro", ms=3, label="Numerov")
    ax.axhline(2 * np.sqrt(8 / G_HARM), color="gray", ls=":", label="oscillator spacing")
    ax.set_xlabel(r"state index $n$")
    ax.set_ylabel(r"$\epsilon_{n+1} - \epsilon_n$")
    ax.set_title(f"Level spacing, slope for $n \\geq 10$: {slope:.3f}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "level_spacing.png", dpi=150)
    return slope


def precision_limits():
    """Error vs. search tolerance, and the truncation/rounding trade-off in N."""
    exact = well_eigenvalue(1, G_WELL)
    tols = np.logspace(-4, -16, 25)
    err_tol = [abs(find_eigenvalue(-1, infinite_well, G_WELL, X, tol=t) - exact) for t in tols]
    ns = np.array([51, 101, 201, 401, 801, 1601, 3201])
    err_n = [abs(find_eigenvalue(-1, infinite_well, G_WELL, grid(n), tol=1e-15) - exact) for n in ns]
    h = 1 / (ns - 1)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    axs[0].loglog(tols, err_tol, "ro", label="error")
    axs[0].loglog(tols, 2 * tols, "k--", label=r"$2\,\mathrm{tol}$ bound")
    axs[0].set_xlabel("search tolerance")
    axs[0].set_ylabel(r"$|\epsilon_1 - \epsilon_1^{\mathrm{exact}}|$")
    axs[0].set_title(r"Tolerance sweep ($N = 1000$)")
    axs[0].legend()
    axs[1].loglog(ns, err_n, "ro", label="error (tol $10^{-15}$)")
    axs[1].loglog(ns[:4], err_n[1] * (h[:4] / h[1]) ** 4, "k--", label=r"truncation $\propto h^4$")
    axs[1].set_xlabel(r"grid points $N$")
    axs[1].set_title("Truncation vs. rounding")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "precision_limits.png", dpi=150)
    return ns[int(np.argmin(err_n))], min(err_n)


if __name__ == "__main__":
    wavefunctions()
    print("convergence order:", round(convergence(), 3))
    uncertainties()
    print("level-spacing slope:", round(level_spacing(), 3))
    print("best N, error:", precision_limits())
