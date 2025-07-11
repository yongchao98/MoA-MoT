import math

# --- Step 1: Define Constants and Input Parameters in SI units ---
# Physical Constants
k_B = 1.380649e-23  # Boltzmann constant in J/K
Torr_to_Pa = 133.322 # Conversion factor for pressure
e_charge = 1.602e-19   # Elementary charge in C
H_ionization_energy_eV = 13.6  # Ionization energy of Hydrogen in eV
H_ionization_energy_J = H_ionization_energy_eV * e_charge # in Joules
sigma0_cm2 = 6.3e-18 # Peak photoionization cross-section in cm^2
sigma0_m2 = sigma0_cm2 * 1e-4 # Peak photoionization cross-section in m^2

# Input Parameters from the problem
T_K = 3000.0  # Temperature in Kelvin
P_Torr = 10.0  # Pressure in Torr
I_W_cm2 = 10.0 # Intensity in W/cm^2

# Convert inputs to SI units
P_Pa = P_Torr * Torr_to_Pa  # Pressure in Pascals (N/m^2)
I_W_m2 = I_W_cm2 * (100**2) # Intensity in W/m^2

print("### Photoelectron Density Estimation ###\n")
print(f"This script calculates the steady-state photoelectron density (n_e) using the equation:")
print(f"n_e = sqrt( (n_H * σ_ph * Φ) / α(T) )\n")
print("--- Calculating each term of the equation ---\n")

# --- Step 2: Calculate the density of Hydrogen atoms (n_H) ---
# Using the Ideal Gas Law: P = n_H * k_B * T  =>  n_H = P / (k_B * T)
n_H = P_Pa / (k_B * T_K)
print(f"1. Hydrogen Atom Density (n_H):")
print(f"   - P = {P_Pa:.2f} Pa")
print(f"   - T = {T_K:.0f} K")
print(f"   - n_H = P / (k_B * T) = {n_H:.3e} atoms/m^3\n")


# --- Step 3: Calculate the Photon Flux (Φ) ---
# The radiation frequency ω ~ e^2/(ħ*a_B) corresponds to the Hartree energy,
# which is 2 * Ry ~ 27.2 eV.
photon_energy_eV = 2 * H_ionization_energy_eV
photon_energy_J = photon_energy_eV * e_charge
photon_flux = I_W_m2 / photon_energy_J
print(f"2. Photon Flux (Φ):")
print(f"   - I = {I_W_m2:.1e} W/m^2")
print(f"   - Photon Energy (E_γ) = {photon_energy_eV:.1f} eV = {photon_energy_J:.3e} J")
print(f"   - Φ = I / E_γ = {photon_flux:.3e} photons/m^2/s\n")


# --- Step 4: Calculate the Photoionization Cross-Section (σ_ph) ---
# Using the approximation σ(E) ≈ σ₀ * (E_ion / E)^(7/2)
energy_ratio = photon_energy_J / H_ionization_energy_J
sigma_ph = sigma0_m2 * (1 / energy_ratio)**(3.5)
print(f"3. Photoionization Cross-Section (σ_ph):")
print(f"   - E_ion / E_γ = 1 / {energy_ratio:.1f}")
print(f"   - σ_ph ≈ σ₀ * (E_ion / E_γ)^(7/2) = {sigma_ph:.3e} m^2\n")


# --- Step 5: Calculate the Recombination Coefficient (α(T)) ---
# Using the approximation α(T) ≈ 2.6e-13 * (T / 10⁴ K)^(-0.7) cm³/s
alpha_recomb_cm3_s = 2.6e-13 * (T_K / 10000)**(-0.7)
alpha_recomb_m3_s = alpha_recomb_cm3_s * 1e-6 # Convert from cm³/s to m³/s
print(f"4. Recombination Coefficient (α(T)):")
print(f"   - T / 10⁴K = {T_K/10000:.2f}")
print(f"   - α(T) ≈ 2.6e-13 * ({T_K/10000:.2f})**(-0.7) cm³/s = {alpha_recomb_m3_s:.3e} m³/s\n")


# --- Step 6: Final Calculation for Electron Density (n_e) ---
# n_e = sqrt( (n_H * σ_ph * Φ) / α(T) )
numerator = n_H * sigma_ph * photon_flux
n_e_squared = numerator / alpha_recomb_m3_s
n_e = math.sqrt(n_e_squared)

print("--- Final Calculation ---")
print("Plugging the values into the final equation:\n")
print(f"n_e = sqrt( (n_H * σ_ph * Φ) / α(T) )")
print(f"    = sqrt( ({n_H:.3e} * {sigma_ph:.3e} * {photon_flux:.3e}) / {alpha_recomb_m3_s:.3e} )")
print(f"    = sqrt( {numerator:.3e} / {alpha_recomb_m3_s:.3e} )")
print(f"    = sqrt( {n_e_squared:.3e} )")
print(f"\nEstimated Photoelectron Density (n_e) = {n_e:.3e} electrons/m^3")
final_answer = f"{n_e:.3e}"