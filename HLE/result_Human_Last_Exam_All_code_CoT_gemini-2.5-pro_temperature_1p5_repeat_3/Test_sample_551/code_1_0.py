import math

# --- 1. Define constants and initial parameters ---
# Physical constants in SI units
k_B = 1.380649e-23  # Boltzmann constant in J/K
e_charge = 1.60217663e-19  # Elementary charge in Coulombs

# Given parameters
T_K = 3000.0  # Temperature in Kelvin
P_torr = 10.0  # Pressure in Torr
I_W_cm2 = 10.0 # Intensity in W/cm^2

# Conversion factors
Pa_per_Torr = 133.322  # Pascals per Torr
m2_per_cm2 = 1e-4      # m^2 per cm^2
m3_per_cm3 = 1e-6      # m^3 per cm^3

# --- 2. Convert all parameters to SI units ---
P_Pa = P_torr * Pa_per_Torr
I_W_m2 = I_W_cm2 / m2_per_cm2

# --- 3. Calculate the density of hydrogen atoms (n_H) using the Ideal Gas Law ---
# P = n_H * k_B * T  =>  n_H = P / (k_B * T)
n_H = P_Pa / (k_B * T_K)

# --- 4. Define and convert physical parameters for the rate equation ---
# Photoionization cross-section for H at threshold (a standard literature value)
sigma_photo_cm2 = 6.3e-18
sigma_photo_m2 = sigma_photo_cm2 * m2_per_cm2

# Photon energy (E_ph), assumed to be the ionization energy of Hydrogen
E_ion_eV = 13.6
E_ph_J = E_ion_eV * e_charge

# Recombination coefficient (alpha) for H+. Using the empirical formula:
# alpha_cm3_s = 2.6e-13 * (T / 10^4 K)^(-0.7)
alpha_cm3_s = 2.6e-13 * (T_K / 10000.0)**(-0.7)
alpha_m3_s = alpha_cm3_s * m3_per_cm3

# --- 5. Solve for the photoelectron density (n_e) ---
# The steady-state equation is: n_e = sqrt( (n_H * sigma_photo * I) / (alpha * E_ph) )
# We calculate the term inside the square root first
term_in_sqrt = (n_H * sigma_photo_m2 * I_W_m2) / (alpha_m3_s * E_ph_J)
n_e_m3 = math.sqrt(term_in_sqrt)

# Convert the final result to electrons per cm^3 for convenience
n_e_cm3 = n_e_m3 * m3_per_cm3

# --- 6. Print the results step-by-step ---
print("--- Calculation of Photoelectron Density ---")
print("The steady-state is reached when the photoionization rate equals the recombination rate.")
print("The governing equation is: n_e = sqrt( (n_H * \u03C3_photo * I) / (\u03B1 * E_ph) )")
print("\n--- Values Used in the Equation (in SI units) ---")
print(f"1. Hydrogen Atom Density (n_H): {n_H:.3e} m\u207b\u00b3")
print(f"2. Photoionization Cross-Section (\u03C3_photo): {sigma_photo_m2:.3e} m\u00b2")
print(f"3. Radiation Intensity (I): {I_W_m2:.3e} W/m\u00b2")
print(f"4. Recombination Coefficient (\u03B1): {alpha_m3_s:.3e} m\u00b3/s")
print(f"5. Photon Energy (E_ph): {E_ph_J:.3e} J")
print("\n--- Final Result ---")
print(f"The estimated density of photoelectrons is {n_e_m3:.3e} electrons/m\u00b3.")
print(f"In more common units, this is {n_e_cm3:.3e} electrons/cm\u00b3.")
