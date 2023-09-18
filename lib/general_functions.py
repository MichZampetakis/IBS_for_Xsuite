import numpy as np


def Hi(beta, alpha, eta, eta_d):
    """Go over with Michalis to figure out what this does."""
    return (1 / beta) * (eta**2 + (beta * eta_d + alpha * eta) ** 2)


def Phi(beta, alpha, eta, eta_d):
    """Go over with Michalis to figure out what this does."""
    return eta_d + alpha * eta / beta


def mean2(numb):
    """Go over with Michalis to figure out what this does."""
    return np.mean((numb - np.mean(numb)) ** 2)


def mean3(numbx, numbpx):
    """Go over with Michalis to figure out what this does."""
    return np.mean((numbx - np.mean(numbx)) * (numbpx - np.mean(numbpx)))


def emittance(x, px):
    """Go over with Michalis to figure out what this does."""
    return np.sqrt(mean2(x) * mean2(px) - mean3(x, px) ** 2)


def Sigma(beta, emit, eta, Sigma_M):
    """Go over with Michalis to figure out what this does."""
    return np.sqrt(beta * emit + (eta * Sigma_M) ** 2)


def BunchLength(
    Circumferance, Harmonic_Num, Energy_total, SlipF, Sigma_E, beta_rel, RF_Voltage, Energy_loss, Z
):
    """
    Go over with Michalis to figure out what this does.
    ~~~ from Wiedermanns book ~~~
    """
    return (
        Sigma_E
        * Circumferance
        * np.sqrt(
            abs(SlipF)
            * Energy_total
            / (2 * np.pi * beta_rel * Harmonic_Num * np.sqrt(Z**2 * RF_Voltage**2 - Energy_loss**2))
        )
    )


def EnergySpread(Circumferance, Harmonic_Num, Energy_total, SlipF, BL, beta_rel, RF_Voltage, Energy_loss, Z):
    """
    Go over with Michalis to figure out what this does.
    ~~~ from Wiedermanns book ~~~
    """
    return BL / (
        Circumferance
        * np.sqrt(
            abs(SlipF)
            * Energy_total
            / (2 * np.pi * beta_rel * Harmonic_Num * np.sqrt(Z**2 * RF_Voltage**2 - Energy_loss**2))
        )
    )


def ion_BunchLength(
    Circumference, Harmonic_Num, Energy_total, SlipF, Sigma_E, beta_rel, RF_Voltage, Energy_loss, Z
):
    """Go over with Michalis to figure out what this does."""
    return (
        Circumference
        / (2.0 * np.pi * Harmonic_Num)
        * np.arccos(
            1
            - (Sigma_E**2 * Energy_total * abs(SlipF) * Harmonic_Num * np.pi)
            / (beta_rel**2 * Z * RF_Voltage)
        )
    )


def ion_EnergySpread(
    Circumference, Harmonic_Num, Energy_total, SlipF, BL, beta_rel, RF_Voltage, Energy_loss, Z
):
    """Go over with Michalis to figure out what this does."""
    tau_phi = 2 * np.pi * Harmonic_Num * BL / Circumference  # bunch length in rad
    return np.sqrt(
        beta_rel**2
        * Z
        * RF_Voltage
        * (-(np.cos(tau_phi) - 1))
        / (Energy_total * abs(SlipF) * Harmonic_Num * np.pi)
    )
