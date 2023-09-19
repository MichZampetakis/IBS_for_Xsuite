import numpy as np
import matplotlib.pylab as plt
import pandas as pd

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rcParams["font.family"] = "serif"
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

kinetic    = pd.read_parquet("xsuite_kinetic.parquet")
simple     = pd.read_parquet("xsuite_simple.parquet")
analytical = pd.read_parquet("xsuite_analytical.parquet")

f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (16,5))

# ax1.plot(nag[0], 'r')
plt.sca(ax1)
plt.plot(kinetic['eps_x'].values, c='b', label='kinetic')
plt.plot(simple['eps_x'].values, c='g', label='simple kick')
plt.plot(analytical['eps_x'].values, c='k', label='analytical')
plt.legend(fontsize=12)

plt.sca(ax2)
plt.plot(kinetic['eps_y'].values, c='b')
plt.plot(simple['eps_y'].values, c='g')
plt.plot(analytical['eps_y'].values, c='k')

plt.sca(ax3)
plt.plot(kinetic['sig_delta'].values*1e3, c='b')
plt.plot(simple['sig_delta'].values*1e3, c='g')
plt.plot(analytical['sig_delta'].values*1e3, c='k')

ax1.set_ylabel(r'$\varepsilon_x$ [m]')
ax1.set_xlabel('Turns')

ax2.set_ylabel(r'$\varepsilon_y$ [m]')
ax2.set_xlabel('Turns')

ax3.set_ylabel(r'$\sigma_{\delta}$ [$10^{-3}$]')
ax3.set_xlabel('Turns')

plt.tight_layout()

plt.savefig('comparison_parquet.png', dpi = 400)
plt.show()
