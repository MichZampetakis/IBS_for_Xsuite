"""
Script to plot Michalis' data which was saved to parquet.
Can be compared to when making modifications.
"""
import matplotlib.pylab as plt
import pandas as pd

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 20,
        "axes.titlesize": 20,
        "axes.labelsize": 20,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "legend.fontsize": 15,
        "figure.titlesize": 20,
    }
)

# ----- Load Data ----- #
kinetic = pd.read_parquet("xsuite_kinetic.parquet")
simple = pd.read_parquet("xsuite_simple.parquet")
analytical = pd.read_parquet("xsuite_analytical.parquet")

figure, (epsx, epsy, sigdelta) = plt.subplots(1, 3, figsize=(16, 5))

# ----- Horizontal Emittance ----- #
epsx.plot(1e10 * kinetic["eps_x"].to_numpy(), label="Kinetic")
epsx.plot(1e10 * simple["eps_x"].to_numpy(), label="Simple")
epsx.plot(1e10 * analytical["eps_x"].to_numpy(), c="magenta", label="Analytical")
epsx.set_ylabel(r"$\varepsilon_x$ [$10^{-10}$m]")
epsx.set_xlabel("Turns")

# ----- Vertical Emittance ----- #
epsy.plot(1e13 * kinetic["eps_y"].to_numpy())
epsy.plot(1e13 * simple["eps_y"].to_numpy())
epsy.plot(1e13 * analytical["eps_y"].to_numpy(), c="magenta")
epsy.set_ylabel(r"$\varepsilon_y$ [$10^{-13}$m]")
epsy.set_xlabel("Turns")

# ----- Longitudinal Emittance ----- #
sigdelta.plot(1e3 * kinetic["sig_delta"].to_numpy())
sigdelta.plot(1e3 * simple["sig_delta"].to_numpy())
sigdelta.plot(1e3 * analytical["sig_delta"].to_numpy(), c="magenta")
sigdelta.set_ylabel(r"$\sigma_{\delta}$ [$10^{-3}$]")
sigdelta.set_xlabel("Turns")

# ----- Final Touches ----- #
epsx.legend(loc="upper left")
plt.tight_layout()
plt.show()
