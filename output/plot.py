import numpy as np
import matplotlib.pylab as plt

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

kin = np.loadtxt('xsuite_kinetic.txt')
sim = np.loadtxt('xsuite_simple.txt')
# nag = np.loadtxt('nagaitsev.txt')


f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (16,5))

# ax1.plot(nag[0], 'r')
ax1.plot(kin.T[0], 'b')
ax1.plot(sim.T[0], 'g')

ax2.plot(kin.T[1], 'b')
ax2.plot(sim.T[1], 'g')

ax3.plot(kin.T[2] * 1e3, 'b')
ax3.plot(sim.T[2] * 1e3, 'g')


ax1.set_ylabel(r'$\varepsilon_x$ [m]')
ax1.set_xlabel('Turns')

ax2.set_ylabel(r'$\varepsilon_y$ [m]')
ax2.set_xlabel('Turns')

ax3.set_ylabel(r'$\sigma_{\delta}$ [$10^{-3}$]')
ax3.set_xlabel('Turns')

plt.tight_layout()

plt.savefig('comparison.png', dpi = 400)
plt.show()