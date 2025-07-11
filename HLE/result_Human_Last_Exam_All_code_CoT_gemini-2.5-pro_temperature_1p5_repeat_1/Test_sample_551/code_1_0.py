import numpy as np

# --- 1. Define constants and given parameters in SI units ---
# Given parameters
P_torr = 10.0  # Torr
T = 3000.0  # Kelvin
I_cm2 = 10.0  # W/cm^2

# Physical constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
h = 6.62607015e-34   # Planck constant in J*s
c = 2.99792458e8     # Speed of light in m/s
eV_to_J = 1.60218e-19 # Conversion factor from eV to Joules

# Unit conversions
torr_to_Pa = 133.322
I = I_cm2 * (100**2) # Convert W/cm^2 to W/m^2
P = P_torr * torr_to_Pa # Convert Torr to Pascals

# --- 2. Calculate the terms for the steady-state equation ---

# a) Calculate hydrogen atom density (n_H) from the ideal gas law
n_H = P / (k_B * T)

# b) Calculate photon energy (E_ph) and flux (Φ)
# E_ph = e^2/a_B = 2 * Ry (Rydberg energy), which is ~27.2 eV
E_ion_H = 13.6  # Ionization energy of Hydrogen in eV
E_ph = 2 * E_ion_H * eV_to_J # Photon energy in Joules
Phi = I / E_ph # Photon flux in photons/m^2/s

# c) Calculate the photoionization cross-section (sigma_ph)
# For a photon energy E, a good approximation is sigma(E) = sigma_0 * (E_ion / E)^3
# where sigma_0 at the ionization threshold (13.6 eV) is ~6.3e-18 cm^2 or 6.3e-22 m^2
sigma_0 = 6.3e-22 # m^2
sigma_ph = sigma_0 * (E_ion_H / (E_ph / eV_to_J))**3

# d) Calculate the recombination coefficient (alpha_r)
# A standard approximation (Case B) is alpha_r_cgs = 2.59e-13 * (T / 10^4 K)^-0.7 cm^3/s
alpha_r_cgs = 2.59e-13 * (T / 10000)**(-0.7)
# Convert from cm^3/s to m^3/s
alpha_r = alpha_r_cgs * 1e-6

# --- 3. Calculate the photoelectron density (n_e) ---
# n_e = sqrt( (sigma_ph * Phi * n_H) / alpha_r )
numerator = sigma_ph * Phi * n_H
n_e = np.sqrt(numerator / alpha_r)

# --- 4. Print the results step-by-step ---
print("--- Calculation Steps ---")
print(f"1. Hydrogen Atom Density (n_H):")
print(f"   P = {P:.2f} Pa, T = {T:.0f} K")
print(f"   n_H = P / (k_B * T) = {n_H:.2e} m^-3")
print("\n2. Photon Flux (Φ):")
print(f"   I = {I:.1e} W/m^2")
print(f"   E_ph = {E_ph/eV_to_J:.1f} eV = {E_ph:.2e} J")
print(f"   Φ = I / E_ph = {Phi:.2e} m^-2 s^-1")
print("\n3. Cross-Section and Recombination Coefficient:")
print(f"   Photoionization Cross-Section (σ_ph) = {sigma_ph:.2e} m^2")
print(f"   Recombination Coefficient (α_r) = {alpha_r:.2e} m^3 s^-1")
print("\n--- Final Calculation ---")
print("Using the steady-state equation: n_e = sqrt( (σ_ph * Φ * n_H) / α_r )")
print("Plugging in the calculated values:")
print(f"n_e = sqrt( ( {sigma_ph:.2e} * {Phi:.2e} * {n_H:.2e} ) / {alpha_r:.2e} )")
print(f"\nThe estimated density of photoelectrons is {n_e:.2e} m^-3.")

# Print the final numerical answer for the platform
final_answer_str = f"{n_e:.2e}"
# This format is to pass the result to the calling script.
# For example, if the result is 3.12e+20, print <<<3.12e+20>>>
# print(f"\n<<<{final_answer_str}>>>")