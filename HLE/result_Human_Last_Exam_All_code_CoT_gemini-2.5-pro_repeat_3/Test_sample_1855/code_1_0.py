import math

# Plan:
# 1. Define all physical constants and given values, rounding constants as required.
# 2. Calculate the total energy of a CMB photon (E_gamma) by summing the primary energy and the perturbative correction.
# 3. Calculate the proton threshold energy (E_p) using the relativistic kinematic formula.
# 4. Print the formula with the numerical values substituted and the final result in the specified format.

# Step 1: Define Constants and Parameters
# Physical constants are rounded to three decimal places as specified.
m_delta = 1.233      # Mass of Delta baryon in GeV
m_p = 0.938          # Mass of proton in GeV (rounded from 0.938272)
k = 8.617e-14        # Boltzmann's constant in GeV/K (rounded from 8.617333e-14)

# Given experimental parameters
T = 2.73             # CMB temperature in Kelvin
T_ref = 2.7          # Reference temperature in Kelvin
prop_coeff_factor = 2.71 # Dimensionless factor for the proportionality constant
pert_coeff = 1.0e-10   # Perturbative coefficient in GeV/K^2

# Step 2: Calculate the mean energy of the CMB photon (E_gamma)
# Primary energy term
E_gamma_primary = prop_coeff_factor * k * T
# Perturbative correction term
delta_T = T - T_ref
E_gamma_pert = pert_coeff * (delta_T**2)
# Total photon energy
E_gamma = E_gamma_primary + E_gamma_pert

# Step 3: Calculate the threshold energy for the proton (E_p)
# The formula for threshold energy E_p is (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p = numerator / denominator

# Step 4: Print the results as requested
print("The formula for the proton threshold energy (E_p) is:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print("\nSubstituting the values:")
# Print the equation with the calculated numbers
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma})")
print(f"E_p = ({m_delta**2} - {m_p**2}) / ({4*E_gamma})")


# Print the final result in scientific notation, rounded to three decimal places
print(f"\nThe calculated threshold energy for the proton is:")
print(f"E_p = {E_p:.3e} GeV")

final_answer = f"{E_p:.3e}"
print(f"<<<{final_answer}>>>")