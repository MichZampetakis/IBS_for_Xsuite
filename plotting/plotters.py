"""
Utilities for the plots
"""
import matplotlib.pylab as plt
import numpy as np

PARAMS = {
    "font.family": "serif",
    "font.size": 20,
    "axes.titlesize": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "legend.fontsize": 15,
    "figure.titlesize": 20,
}


def plot_emittances_and_momentum_spread_evolution(
    kinetic_epsx: np.ndarray,
    kinetic_epsy: np.ndarray,
    kinetic_sig_delta: np.ndarray,
    simple_epsx: np.ndarray,
    simple_epsy: np.ndarray,
    simple_sig_delta: np.ndarray,
    analytical_epsx: np.ndarray,
    analytical_epsy: np.ndarray,
    analytical_sig_delta: np.ndarray,
) -> plt.Figure:
    figure, (epsx, epsy, sigdelta) = plt.subplots(1, 3, figsize=(16, 5))

    # ----- Horizontal Emittance ----- #
    epsx.plot(1e10 * kinetic_epsx, label="Kinetic")
    epsx.plot(1e10 * simple_epsx, label="Simple")
    epsx.plot(1e10 * analytical_epsx, c="magenta", label="Analytical")
    epsx.set_ylabel(r"$\varepsilon_x$ [$10^{-10}$m]")
    epsx.set_xlabel("Turns")

    # ----- Vertical Emittance ----- #
    epsy.plot(1e13 * kinetic_epsy)
    epsy.plot(1e13 * simple_epsy)
    epsy.plot(1e13 * analytical_epsy, c="magenta")
    epsy.set_ylabel(r"$\varepsilon_y$ [$10^{-13}$m]")
    epsy.set_xlabel("Turns")

    # ----- Longitudinal Emittance ----- #
    sigdelta.plot(1e3 * kinetic_sig_delta)
    sigdelta.plot(1e3 * simple_sig_delta)
    sigdelta.plot(1e3 * analytical_sig_delta, c="magenta")
    sigdelta.set_ylabel(r"$\sigma_{\delta}$ [$10^{-3}$]")
    sigdelta.set_xlabel("Turns")

    # ----- Final Touches ----- #
    epsx.legend(loc="upper left")
    plt.tight_layout()
    return figure
