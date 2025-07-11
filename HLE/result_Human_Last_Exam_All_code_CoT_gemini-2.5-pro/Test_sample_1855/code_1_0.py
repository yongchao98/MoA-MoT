import math

# Step 1: Define physical constants and given values.
# Physical constants are rounded to three decimal places as requested.
# All other calculations will retain full precision.

# Proton mass in GeV
m_p = 0.938
# Boltzmann's constant in GeV/K
k = 8.617e-14
# Delta baryon mass in GeV (given)
m_delta = 1.233
# CMB temperature in Kelvin (given)
T = 2.73
# Base temperature for perturbation in Kelvin (given)
T0 = 2.7
# Proportionality coefficient for the main energy term (given)
C_prop = 2.71
# Perturbative coefficient in GeV/K^2 (given as O(10^-10))
C_pert = 1.0e-10

print("--- Step 1: Define Constants and Given Values ---")
print(f"Proton mass (m_p): {m_p} GeV")
print(f"Boltzmann's constant (k): {k:.3e} GeV/K")
print(f"Delta baryon mass (m_delta): {m_delta} GeV")
print(f"CMB temperature (T): {T} K")
print("\n")

# Step 2: Calculate the mean energy of a CMB photon (E_gamma).
print("--- Step 2: Calculate Mean Photon Energy (E_gamma) ---")
# The model for photon energy is E_gamma = C_prop * k * T + C_pert * (T - T0)^2
delta_T = T - T0
main_energy_term = C_prop * k * T
pert_energy_term = C_pert * (delta_T)**2
E_gamma = main_energy_term + pert_energy_term

print(f"The formula for E_gamma is: (C_prop * k * T) + (C_pert * (T - T0)^2)")
print(f"E_gamma = ({C_prop} * {k:.3e} * {T}) + ({C_pert:.1e} * ({T} - {T0})^2)")
print(f"E_gamma = {main_energy_term:.5e} GeV + {pert_energy_term:.5e} GeV")
print(f"E_gamma = {E_gamma:.5e} GeV")
print("\n")


# Step 3: Calculate the threshold proton energy (E_p).
print("--- Step 3: Calculate Proton Threshold Energy (E_p) ---")
# The formula is E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p = numerator / denominator

print("The formula for E_p is: (m_delta^2 - m_p^2) / (4 * E_gamma)")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.5e})")
print(f"E_p = ({m_delta**2:.5f} - {m_p**2:.5f}) / {4 * E_gamma:.5e}")
print(f"E_p = {numerator:.5f} / {denominator:.5e}")
print(f"Calculated E_p = {E_p:.5e} GeV")
print("\n")

# Step 4: Format and present the final answer.
print("--- Step 4: Final Answer ---")
# Format the result to scientific notation with three decimal places.
final_answer_formatted = f"{E_p:.3e}"
print(f"The threshold energy for the proton is {final_answer_formatted} GeV.")

print(f"<<<{final_answer_formatted}>>>")