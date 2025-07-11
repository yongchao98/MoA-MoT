import math

# This script calculates the proton threshold energy for the reaction p + gamma -> Delta.
# It first calculates the energy of the CMB photon and then uses that to find the proton energy.

# Step 1: Define physical constants and given values.
# Physical constants are rounded to three decimal places as per the instructions.
m_delta_gev = 1.233      # Mass of Delta baryon in GeV
m_p_gev = 0.938          # Mass of proton in GeV (0.938272... rounded to 3 decimal places)
k_boltzmann_gev_k = 8.617e-14  # Boltzmann's constant in GeV/K (8.617333...e-14 rounded to 3 decimal places)

# Given experimental parameters
T_k = 2.73              # Temperature in Kelvin
T_ref_k = 2.7           # Reference temperature in Kelvin for the perturbation
c1_coeff = 2.71         # Proportionality coefficient for CMB energy
c2_coeff = 1.0e-10      # Perturbative coefficient in GeV/K^2

# Step 2: Calculate the average energy of a CMB photon (E_gamma) at T = 2.73 K.
print("Step 1: Calculating the average CMB photon energy (E_gamma)")
delta_T = T_k - T_ref_k
# The formula is E_gamma = (2.71 * k * T) + C_pert * (T - 2.7)^2
E_gamma_primary = c1_coeff * k_boltzmann_gev_k * T_k
E_gamma_perturbation = c2_coeff * (delta_T**2)
E_gamma_gev = E_gamma_primary + E_gamma_perturbation
print(f"The energy of a CMB photon is calculated as E_gamma = ({c1_coeff} * {k_boltzmann_gev_k} * {T_k}) + ({c2_coeff} * ({T_k} - {T_ref_k})^2)")
print(f"E_gamma = {E_gamma_gev} GeV\n")


# Step 3: Calculate the proton threshold energy (E_p).
print("Step 2: Calculating the proton threshold energy (E_p)")
print("The threshold energy is calculated using the formula: E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")

# Calculate the numerator and denominator for the final equation
m_delta_sq = m_delta_gev**2
m_p_sq = m_p_gev**2
numerator = m_delta_sq - m_p_sq
denominator = 4 * E_gamma_gev

# Calculate the final result
E_p_thresh_gev = numerator / denominator

# Step 4: Display the final equation with all numbers and the result.
print("\nSubstituting the numerical values into the equation:")
print(f"E_p = ({m_delta_gev}^2 - {m_p_gev}^2) / (4 * {E_gamma_gev})")
print(f"E_p = ({m_delta_sq} - {m_p_sq}) / ({denominator})")
print(f"E_p = {numerator} / {denominator}")

# Format the final answer into scientific notation with three decimal places.
final_answer_sci = f"{E_p_thresh_gev:.3e}"
print(f"\nThe resulting proton threshold energy is {final_answer_sci} GeV.")

# Final answer in the specified format.
print(f"<<<{final_answer_sci}>>>")