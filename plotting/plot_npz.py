"""
Script to plot my own new data which was saved to npz files.
Can be a first quick comparision when making modifications.
"""
import matplotlib.pylab as plt
import numpy as np
from plotters import PARAMS, plot_emittances_and_momentum_spread_evolution

plt.rcParams.update(PARAMS)

# ----- Load Data ----- #
with np.load("../outputs/xsuite_kinetic.npz") as data:
    kinetic_epsx = data["epsilon_x"]
    kinetic_epsy = data["epsilon_y"]
    kinetic_sig_delta = data["sig_delta"]

with np.load("../outputs/xsuite_simple.npz") as data:
    simple_epsx = data["epsilon_x"]
    simple_epsy = data["epsilon_y"]
    simple_sig_delta = data["sig_delta"]

with np.load("../outputs/xsuite_analytical.npz") as data:
    analytical_epsx = data["epsilon_x"]
    analytical_epsy = data["epsilon_y"]
    analytical_sig_delta = data["sig_delta"]

figure = plot_emittances_and_momentum_spread_evolution(
    kinetic_epsx,
    kinetic_epsy,
    kinetic_sig_delta,
    simple_epsx,
    simple_epsy,
    simple_sig_delta,
    analytical_epsx,
    analytical_epsy,
    analytical_sig_delta,
)
plt.tight_layout()
plt.show()
