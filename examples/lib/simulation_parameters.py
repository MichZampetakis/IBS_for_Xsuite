import numpy as np

from scipy.constants import c, e, epsilon_0, m_p, physical_constants

# TODO: this is not called elsewhere, can go into some constants module in pyibs

# r0 = physical_constants["classical electron radius"][0]                     #Electron Classical Radius [m]
E0p = physical_constants["proton mass energy equivalent in MeV"][0] * 1e-3  # Proton Rest Mass [Gev]

Am = 208.0  # Atomic Mass (207.2)
Z = 54.0  # Number of Charges
Nb = 2e8  # Bunch Population
h = 2  # Harmonic Number

E0i = 193.6999  # Rest energy in GeV
Ek = 0.0042  # Kinetic Energy  [GeV]
En = E0i + Ek * Am  # Ion Total Energy [GeV]
P0 = np.sqrt(En**2 - E0i**2) / c  # Reference momentum [GeV/c]
mi = (E0i * m_p) / E0p  # Ions Mass

gamma_rel = En / E0i
beta_rel = np.sqrt(1 - 1 / gamma_rel**2)

gamma_tr = 2.838351
a_p = 1 / gamma_tr**2
SlipF = (1 / gamma_tr**2) - (1 / gamma_rel**2)  # Slip Factor

V0 = 1.1e-6  # RF-Voltage [GV]
U0 = 0.0  # Energy Loss / Turn [Gev]
p_increment = U0 * 1e9 * 1.602 * 1e-19 / c  # kg*m/s //turn

r0 = (Z * e) ** 2 / (4 * np.pi * epsilon_0 * c**2 * mi)  # Ion Classical Radius

# ! -------------------------------------------- !

turns = 16000
IBS_step = 50

# ! ---------- Equilibrium ------------ !
emit_x0 = 2.95814809e-6
emit_y0 = 2.95814809e-6
Nemit_x0 = emit_x0 * beta_rel * gamma_rel
Nemit_y0 = emit_y0 * beta_rel * gamma_rel
bl0 = 4.256474951
# Sigma_E0 = 1.05834e-05
