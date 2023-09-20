"""
Script to plot Michalis' data which was saved to parquet.
Can be compared to when making modifications.
"""
import matplotlib.pylab as plt
import pandas as pd
from plotters import PARAMS, plot_emittances_and_momentum_spread_evolution

plt.rcParams.update(PARAMS)

# ----- Load Data ----- #
kinetic = pd.read_parquet("../outputs/xsuite_kinetic.parquet")
simple = pd.read_parquet("../outputs/xsuite_simple.parquet")
analytical = pd.read_parquet("../outputs/xsuite_analytical.parquet")

figure = plot_emittances_and_momentum_spread_evolution(
    kinetic.eps_x.to_numpy(),
    kinetic.eps_y.to_numpy(),
    kinetic.sig_delta.to_numpy(),
    simple.eps_x.to_numpy(),
    simple.eps_y.to_numpy(),
    simple.sig_delta.to_numpy(),
    analytical.eps_x.to_numpy(),
    analytical.eps_y.to_numpy(),
    analytical.sig_delta.to_numpy(),
)

plt.tight_layout()
plt.show()
