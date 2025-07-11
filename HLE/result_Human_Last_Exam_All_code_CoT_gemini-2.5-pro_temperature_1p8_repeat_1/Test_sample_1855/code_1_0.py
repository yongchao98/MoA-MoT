import math

# Step 1: Define physical constants and problem parameters.
# Constants are rounded to three decimal places as requested.
m_delta = 1.233      # Mass of the Delta baryon in GeV
m_p = 0.938          # Mass of the proton in GeV (rounded from 938.272 MeV)
k = 8.617e-14        # Boltzmann's constant in GeV/K (rounded from 8.61733...)

# Temperature parameters and coefficients for the photon energy model
T = 2.73             # CMB temperature in Kelvin
T0 = 2.7             # Reference temperature in Kelvin
c1_k_coeff = 2.71      # Coefficient for the linear energy term
c2_pert_coeff = 1.0e-10 # Coefficient for the perturbative energy term in GeV/K^2

# Step 2: Calculate the average energy of a CMB photon (E_gamma).
# The model includes a linear term and a perturbative quadratic term.
delta_T = T - T0
linear_term = c1_k_coeff * k * T
perturbative_term = c2_pert_coeff * (delta_T)**2
E_gamma = linear_term + perturbative_term

# Step 3: Calculate the proton's threshold energy (E_p_thresh).
# This is derived from the Lorentz-invariant s = (p_p + p_gamma)^2 = m_delta^2
# In a head-on collision, this simplifies to E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
m_delta_sq = m_delta**2
m_p_sq = m_p**2
numerator = m_delta_sq - m_p_sq
denominator = 4 * E_gamma
E_p_thresh = numerator / denominator

# Step 4: Print the calculation and the final answer.
# The problem asks to show the equation with the numbers plugged in.
print("Equation for the proton threshold energy:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print("\nPlugging in the values:")
print(f"E_p = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma})")
print(f"E_p = ({m_delta_sq} - {m_p_sq}) / ({denominator})")
print(f"E_p = ({numerator}) / ({denominator})")

# Print the final result in scientific notation to three decimal places.
print(f"\nThe average threshold energy for the proton is {E_p_thresh:.3e} GeV.")

# Final answer in the specified format
final_answer_string = f"{E_p_thresh:.3e}"
print(f"\n<<<{final_answer_string}>>>")