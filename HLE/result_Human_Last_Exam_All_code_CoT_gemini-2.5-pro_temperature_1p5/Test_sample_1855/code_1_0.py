import math

# Plan:
# 1. Define the physical constants, rounded as specified.
# 2. Calculate the mean CMB photon energy with the perturbative correction.
# 3. Calculate the proton threshold energy using the relativistic invariant formula.
# 4. Print the breakdown of the equation and the final formatted answer.

# Step 1: Define Constants
# The mass of the delta baryon is given in GeV.
m_delta = 1.233  # GeV
# The mass of the proton is rounded to three decimal places.
m_p = 0.938  # GeV
# Boltzmann's constant in GeV/K, rounded to three decimal places for the significand.
k = 8.617e-14  # GeV/K

# Temperature and perturbative correction parameters
T = 2.73  # K
T_base = 2.7  # K
c1_factor = 2.71
c2_coeff = 1.0e-10  # GeV/K^2

# Step 2: Calculate the mean CMB photon energy (E_gamma)
# E_gamma = (2.71 * k * T) + (1e-10) * (delta_T)^2
delta_T = T - T_base
E_gamma = c1_factor * k * T + c2_coeff * (delta_T**2)

# Step 3: Calculate the proton threshold energy (E_p)
# The formula is E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p = numerator / denominator

# Step 4: Output the equation breakdown and the final result
print("The threshold energy for the proton (E_p) is found using the formula derived from relativistic kinematics:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)\n")

print("The values used in the calculation are:")
print(f"m_delta (Δ baryon mass) = {m_delta} GeV")
print(f"m_p (proton mass)         = {m_p} GeV")
print(f"E_gamma (photon energy)   = {E_gamma:.6e} GeV\n")

print("Substituting these values into the equation:")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma})")
print(f"E_p = ({numerator}) / ({denominator})")
print(f"E_p = {E_p} GeV\n")

# Format the final answer into scientific notation with three decimal places
formatted_E_p = f"{E_p:.3e}"
print("Final Answer:")
print(f"The average threshold energy for the proton is {formatted_E_p} GeV.")
