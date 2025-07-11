import math

# Step 1: Define constants and problem parameters in SI units
P_torr = 10.0  # Torr
T_K = 3000.0  # Kelvin
I_W_cm2 = 10.0  # W/cm^2

# Physical Constants (SI units)
TORR_TO_PA = 133.322  # Pa per Torr
K_B = 1.380649e-23     # Boltzmann constant, J/K
E_CHARGE = 1.602177e-19 # Elementary charge, C
I_H_eV = 13.6          # Ionization energy of Hydrogen, eV

# Convert parameters to SI units
P_Pa = P_torr * TORR_TO_PA               # Pressure in Pascals
I_W_m2 = I_W_cm2 * 10000.0             # Intensity in W/m^2
I_H_J = I_H_eV * E_CHARGE              # Ionization energy in Joules

# Step 2: Calculate the terms needed for the steady-state equation

# a) Calculate the density of hydrogen atoms (n_H) using the ideal gas law
n_H = P_Pa / (K_B * T_K)

# b) Calculate the energy of the incident photons (E_photon)
# From the problem statement, E_photon = 2 * I_H
E_photon = 2 * I_H_J

# c) Calculate the photon flux (F)
F = I_W_m2 / E_photon

# d) Estimate the photoionization cross-section (sigma_ph) for H at E=2*I_H
# The value from literature at this energy is ~1.0e-18 cm^2.
sigma_ph_cm2 = 1.0e-18
sigma_ph_m2 = sigma_ph_cm2 * 1e-4 # Convert cm^2 to m^2

# e) Estimate the radiative recombination coefficient (alpha_r)
# We assume electron temperature T_e = T_gas = 3000 K.
# A standard approximation is alpha_r ~ 5.3e-13 cm^3/s at 3000 K.
alpha_r_cm3_s = 5.3e-13
alpha_r_m3_s = alpha_r_cm3_s * 1e-6 # Convert cm^3 to m^3

# Step 3: Solve the steady-state equation for electron density (n_e)
# n_e = sqrt((sigma_ph * F * n_H) / alpha_r)
numerator = sigma_ph_m2 * F * n_H
n_e_squared = numerator / alpha_r_m3_s
n_e = math.sqrt(n_e_squared)

# Step 4: Print the results, including the final equation with values
print("--- Calculated Values for the Steady-State Equation ---")
print(f"Hydrogen Density (n_H): {n_H:.3e} m^-3")
print(f"Photon Flux (F): {F:.3e} photons/m^2/s")
print(f"Photoionization Cross-section (sigma_ph): {sigma_ph_m2:.3e} m^2")
print(f"Recombination Coefficient (alpha_r): {alpha_r_m3_s:.3e} m^3/s")
print("\n--- Final Equation ---")
print("n_e = sqrt((sigma_ph * F * n_H) / alpha_r)")
print(f"n_e = sqrt(({sigma_ph_m2:.3e} m^2 * {F:.3e} m^-2 s^-1 * {n_H:.3e} m^-3) / {alpha_r_m3_s:.3e} m^3 s^-1)")

print("\n--- Final Result ---")
print(f"The estimated density of photoelectrons is: {n_e:.3e} m^-3")

# Extract the numerical value for the final answer format
final_answer = f"{n_e:.2e}"
# print(f"\n<<<{final_answer}>>>")
print(f"\n<<<{n_e:.3e}>>>")