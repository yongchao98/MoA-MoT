import numpy as np

# This script estimates the density of photoelectrons in a cuvette of atomic hydrogen
# by balancing the rates of photoionization and recombination in a steady state.

# --- Plan ---
# 1. Define physical constants in SI units.
# 2. Define input parameters from the problem and convert to SI units.
# 3. Formulate the problem: In steady-state, the rate of photoionization equals the rate of recombination.
#    This gives the equation: n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )
#    where n_e is the electron density, σ_ion is the photoionization cross-section,
#    Φ is the photon flux, n_H is the hydrogen atom density, and α_rec is the recombination coefficient.
# 4. Calculate each term in the equation.
# 5. Substitute the calculated values to find n_e.
# 6. Print the breakdown of the calculation and the final result as an equation with numbers.

# Step 1: Physical constants in SI units
k_B = 1.380649e-23  # Boltzmann constant (J/K)
E_ion_eV = 13.6     # Ionization energy of Hydrogen (eV)
E_ion_J = E_ion_eV * 1.6021766e-19 # Ionization energy of Hydrogen (J)
# Photoionization cross-section at threshold for Hydrogen
sigma_ion_thresh = 6.3e-22 # m^2

# Step 2: Input parameters from the problem in SI units
T = 3000.0          # Temperature (K)
P_torr = 10.0       # Pressure (Torr)
P_pa = P_torr * 133.322 # Pressure in Pascals (N/m^2)
I_cm2 = 10.0        # Intensity (W/cm^2)
I_m2 = I_cm2 * 1e4  # Intensity in W/m^2

# Step 3 & 4: Calculate the components of the steady-state equation

# a) Hydrogen atom density (n_H) from the ideal gas law (P = n*k_B*T)
n_H = P_pa / (k_B * T)

# b) Photon flux (Φ)
# The given UV frequency ω ~ e^2/(ħ*a_B) corresponds to a photon energy
# E_photon = ħω ~ e^2/a_B, which is twice the ionization energy of hydrogen.
E_photon = 2 * E_ion_J
photon_flux = I_m2 / E_photon

# c) Photoionization cross-section (σ_ion)
# The cross-section scales approximately as (E_ion / E_photon)^-3 from the threshold.
sigma_ion = sigma_ion_thresh * (E_ion_J / E_photon)**3

# d) Recombination coefficient (α_rec)
# Using a standard formula for radiative recombination at electron temperature T:
# α_rec(T) ≈ 2.6e-19 * (T / 10^4 K)^(-0.7) m^3/s
alpha_rec = 2.6e-19 * (T / 1e4)**(-0.7)

# Step 5: Calculate the final electron density (n_e)
# We solve the steady-state equation: n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )
n_e = np.sqrt((sigma_ion * photon_flux * n_H) / alpha_rec)

# Step 6: Print the final result in the requested format
print("To estimate the photoelectron density (n_e), we balance the photoionization and recombination rates.")
print("The steady-state equation is: n_e = sqrt( (σ_ion * Φ * n_H) / α_rec )\n")
print("Calculated values for each term:")
print(f"  Photoionization cross-section (σ_ion): {sigma_ion:.3e} m^2")
print(f"  Photon flux (Φ):                       {photon_flux:.3e} m^-2 s^-1")
print(f"  Hydrogen atom density (n_H):         {n_H:.3e} m^-3")
print(f"  Recombination coefficient (α_rec):   {alpha_rec:.3e} m^3 s^-1\n")
print("Substituting these numbers into the equation:")
print(f"n_e = sqrt( ({sigma_ion:.3e} * {photon_flux:.3e} * {n_H:.3e}) / {alpha_rec:.3e} )")
print(f"\nThe estimated density of photoelectrons is:\nn_e = {n_e:.3e} m^-3")
