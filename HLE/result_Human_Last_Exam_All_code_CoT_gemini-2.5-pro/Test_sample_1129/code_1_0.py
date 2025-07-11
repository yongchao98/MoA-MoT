import math

# --- Define Constants and Parameters ---

# Lawson parameter (L) for D-T fusion, used as an optimistic value for this hypothetical reaction.
# Units: keV * s / m^3
lawson_parameter_L = 3e21

# Particle confinement time (tau), given as 'optimistic'.
# Units: seconds (s)
confinement_time_tau = 1.0

# Temperature in Celsius, given in the problem.
temp_celsius = 20.0

# Side length of the cube-shaped reaction chamber.
# Units: meters (m)
chamber_side_length_m = 0.1

# --- Step-by-step Calculation ---

# 1. Convert temperature from Celsius to Kelvin.
temp_kelvin = temp_celsius + 273.15

# 2. Convert temperature from Kelvin to kiloelectronvolts (keV).
# Boltzmann constant (kB) in eV/K.
boltzmann_constant_eV_K = 8.617333e-5
temp_eV = temp_kelvin * boltzmann_constant_eV_K
temp_keV = temp_eV / 1000

# 3. Calculate the volume of the reaction chamber.
# Units: cubic meters (m^3)
volume_m3 = chamber_side_length_m ** 3

# 4. Calculate the minimum required particle density (n_min) using the Lawson criterion.
# n_min = L / (T * tau)
# Units: atoms / m^3
min_density_n = lawson_parameter_L / (temp_keV * confinement_time_tau)

# 5. Calculate the minimum number of atoms (N_min) required.
# N_min = n_min * V
min_atoms_N = min_density_n * volume_m3

# --- Final Output ---
# The prompt requires printing the final equation with all numbers.
# We will format the output to show the full calculation for the number of atoms (N).
# N = n * V = (L / (T * tau)) * V
print("This script calculates the minimum number of atoms (N) for a sustained fusion reaction based on the Lawson Criterion.")
print("\n--- Calculation ---")
print(f"Minimum Atoms (N) = (Lawson Parameter / (Temperature_keV * Confinement Time)) * Volume")
print(f"N = ({lawson_parameter_L:.2e} keV·s/m³) / ({temp_keV:.4e} keV * {confinement_time_tau:.1f} s) * {volume_m3:.3f} m³")
print(f"N = {min_atoms_N:.4e} atoms")

# The final numerical answer in the required format.
# We use scientific notation for the large number.
final_answer = f"{min_atoms_N:.4e}"
print(f"\n<<<>>>")
print(f"The minimum number of atoms required is approximately {final_answer}.")
print(f"<<<{final_answer}>>>")