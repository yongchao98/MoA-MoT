import math

# Step 1: Define constants and given values
# Physical constants are rounded to three decimal places as requested.
# All other calculations retain full precision based on these constants.

# Given masses in GeV
m_delta = 1.233  # Delta baryon mass
m_p = 0.938      # Proton mass (PDG value ~0.938272 GeV, rounded)

# Temperature values in Kelvin
T = 2.73   # Current CMB temperature
T_0 = 2.7    # Base temperature for perturbation

# Boltzmann's constant in GeV/K, rounded to three decimal places
# The full value is k ≈ 8.617333 x 10^-14 GeV/K
k = 8.617e-14  # GeV/K

# Coefficients for photon energy calculation from the problem description
c1_coeff = 2.71       # Coefficient for the linear term
c2_coeff = 1.0e-10    # Coefficient for the perturbative term [GeV/K^2]

print("This script calculates the threshold energy for a proton reacting with a CMB photon.")
print("-" * 70)

# Step 2: Calculate the mean CMB photon energy (E_gamma)
delta_T = T - T_0
linear_term = c1_coeff * k * T
perturbative_term = c2_coeff * (delta_T**2)
E_gamma = linear_term + perturbative_term

print("1. Calculation of the mean CMB photon energy (E_gamma):")
print(f"   The model for photon energy is: E_gamma = ({c1_coeff} * k * T) + (c2_coeff * (T - T_0)^2)")
print(f"   Plugging in the values:")
print(f"   E_gamma = ({c1_coeff} * {k:.3e} GeV/K * {T} K) + ({c2_coeff:.1e} GeV/K^2 * ({T} K - {T_0} K)^2)")
print(f"   E_gamma = {linear_term:.5e} GeV + {perturbative_term:.5e} GeV")
print(f"   E_gamma = {E_gamma:.5e} GeV\n")

# Step 3: Calculate the proton threshold energy (E_th)
m_delta_sq = m_delta**2
m_p_sq = m_p**2
m_diff_sq = m_delta_sq - m_p_sq

# We use the full formula for E_th for completeness.
term1 = m_diff_sq / (4 * E_gamma)
term2 = (m_p_sq * E_gamma) / m_diff_sq
E_th = term1 + term2

print("2. Calculation of the proton threshold energy (E_th):")
print(f"   The formula for threshold energy is: E_th = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print(f"   Plugging in the values:")
print(f"   E_th = ({m_delta}^2 GeV^2 - {m_p}^2 GeV^2) / (4 * {E_gamma:.5e} GeV)")
print(f"   E_th = ({m_delta_sq:.5f} - {m_p_sq:.5f}) GeV^2 / ({4 * E_gamma:.5e} GeV)")
print(f"   E_th = {m_diff_sq:.5f} GeV^2 / {4 * E_gamma:.5e} GeV")
# The second term of the more precise formula is negligible, but we show its calculation.
# print(f"   E_th_term1 = {term1:.5e} GeV")
# print(f"   E_th_term2 = {term2:.5e} GeV")
print(f"   E_th = {E_th:.5e} GeV\n")

# Step 4: Format and print the final answer
# The final answer is required in scientific notation, rounded to three decimal places.
final_answer_sci = f"{E_th:.3e}"
print("-" * 70)
print(f"The final calculated average threshold energy for the proton is {final_answer_sci} GeV.")

# Final answer in the specified format for automated grading.
print(f"<<<{final_answer_sci}>>>")