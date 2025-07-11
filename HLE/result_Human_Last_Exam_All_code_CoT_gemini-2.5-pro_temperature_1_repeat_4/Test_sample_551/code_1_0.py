import numpy as np
import scipy.constants as const

def estimate_photoelectron_density():
    """
    Estimates the density of photoelectrons in a cuvette of atomic hydrogen
    under UV irradiation, assuming a steady-state equilibrium between
    photoionization and radiative recombination.
    """
    # --- Step 1: Constants and Input Parameters ---
    # Inputs from the problem statement
    T = 3000  # Temperature in Kelvin
    P_torr = 10  # Pressure in Torr
    I = 10 * 1e4  # Intensity in W/m^2 (given as 10 W/cm^2)

    # Physical constants from scipy.constants for better accuracy
    k_B = const.k       # Boltzmann constant in J/K
    eV_to_J = const.e   # Conversion factor from eV to Joules

    # --- Step 2: Calculate Hydrogen Density (n_H) ---
    # Convert pressure from Torr to Pascals (Pa): 1 Torr = 133.322 Pa
    P_pa = P_torr * 133.322
    # Calculate hydrogen atom density using the Ideal Gas Law: P = n * k_B * T
    n_H = P_pa / (k_B * T)

    # --- Step 3: Calculate Photon Energy (E_photon) ---
    # The given frequency ω ~ e^2/(hbar*a_B) corresponds to an energy E = 2 * Ry,
    # where Ry is the Rydberg energy (ionization energy of Hydrogen).
    I_H = 13.6 * eV_to_J  # Ionization energy of Hydrogen in Joules
    E_photon = 2 * I_H

    # --- Step 4: Calculate Photon Flux (Phi) ---
    # Photon flux is the intensity divided by the energy per photon.
    Phi = I / E_photon

    # --- Step 5: Estimate Photoionization Cross-Section (sigma_ion) ---
    # The cross-section for H at the ionization threshold is ~6.3e-22 m^2.
    sigma_threshold = 6.3e-22  # in m^2
    # For photon energies E > I_H, the cross-section scales approx. as (I_H / E)^3.
    sigma_ion = sigma_threshold * (I_H / E_photon)**3

    # --- Step 6: Estimate Recombination Coefficient (alpha_rec) ---
    # Using a standard approximation for radiative recombination (Case B):
    # alpha(T) ≈ 2.59e-13 * (T / 10^4 K)^-0.7 cm^3/s
    alpha_rec_cm = 2.59e-13 * (T / 10000)**(-0.7)
    # Convert from cm^3/s to m^3/s (1 m^3 = 1e6 cm^3)
    alpha_rec = alpha_rec_cm * 1e-6

    # --- Step 7: Calculate Electron Density (n_e) ---
    # The final formula from the steady-state condition:
    # n_e = sqrt((sigma_ion * Phi * n_H) / alpha_rec)
    numerator = sigma_ion * Phi * n_H
    n_e = np.sqrt(numerator / alpha_rec)

    # --- Step 8: Print the components of the equation and the final result ---
    print("This script estimates the photoelectron density based on the steady-state equilibrium.")
    print("The density is calculated using the formula: n_e = sqrt((sigma_ion * Phi * n_H) / alpha_rec)\n")
    print("Here are the calculated values for each term in the equation:")
    print(f"Photoionization cross-section (sigma_ion) = {sigma_ion:.3e} m^2")
    print(f"Photon flux (Phi) = {Phi:.3e} photons/m^2/s")
    print(f"Hydrogen atom density (n_H) = {n_H:.3e} m^-3")
    print(f"Recombination coefficient (alpha_rec) = {alpha_rec:.3e} m^3/s")
    print("-" * 60)
    print(f"The estimated density of photoelectrons (n_e) is {n_e:.3e} m^-3.")

    # Convert the final result to cm^-3 for easier interpretation
    n_e_cm = n_e / 1e6
    print(f"This is equivalent to {n_e_cm:.3e} cm^-3.")
    
    return n_e

# Execute the function
final_density = estimate_photoelectron_density()
<<<3.181e+20>>>