import numpy as np

# TODO: this is not called elsewhere, can go into some utils module in pyibs


def drift(eta, C0, z, delta):
    z1 = z - eta * delta * C0
    delta1 = delta
    return z1, delta1


def RF_map(P0, E0, beta_0, f_RF, c, lag, V_RF, En, z, delta):
    P = (1 + delta) * P0  # GeV/c
    gamma = np.sqrt(1 + (P * c / E0) ** 2)
    beta = np.sqrt(1 - 1 / gamma**2)
    old_rvv = beta / beta_0

    # RF energy kick
    k = 2 * np.pi * f_RF / c
    tau = z / old_rvv / beta_0
    phase = lag * np.pi / 180.0 - k * tau
    energy = V_RF * np.sin(phase)

    deltabeta_0 = delta * beta_0
    ptaubeta_0 = np.sqrt(deltabeta_0**2 + 2 * deltabeta_0 * beta_0 + 1) - 1
    ptaubeta_0 += energy / En  # energy kick
    ptau = ptaubeta_0 / beta_0

    delta1 = np.sqrt(ptau**2 + 2 * ptau / beta_0 + 1) - 1
    rvv = (1 + delta1) / (1 + ptaubeta_0)
    z1 = z * rvv / old_rvv
    return z1, delta1
