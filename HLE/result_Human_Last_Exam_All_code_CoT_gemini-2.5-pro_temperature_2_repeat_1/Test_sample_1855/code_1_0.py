import math

# --- 1. Define Constants and Given Values ---
# Given values from the problem
m_delta = 1.233  # Mass of the Delta baryon in GeV
T = 2.73         # Temperature of the CMB in K
T0 = 2.7         # Base temperature in K
c1_coeff = 2.71  # Dimensionless coefficient for the linear energy term
c2 = 1e-10       # Perturbative coefficient in GeV/K^2

# Physical constants, rounded to three decimal places as specified
# Boltzmann constant in GeV/K
k = 8.617e-14
# Proton mass in GeV
m_p = 0.938

# --- 2. Calculate Mean Photon Energy (E_gamma) ---
# Calculate the temperature difference
delta_T = T - T0

# Calculate the mean photon energy using the provided model
# E_gamma = (coefficient * k * T) + (perturbative_coeff * delta_T^2)
E_gamma_linear_term = c1_coeff * k * T
E_gamma_pert_term = c2 * (delta_T**2)
E_gamma = E_gamma_linear_term + E_gamma_pert_term

print("--- Step 1: Calculating the average photon energy (E_gamma) ---")
print(f"The formula used is: E_gamma = ({c1_coeff} * k * T) + ({c2} * (T - {T0})^2)")
print("Using values:")
print(f"  k = {k:.3e} GeV/K")
print(f"  T = {T} K")
print(f"  delta_T = {delta_T:.2f} K")
print(f"The calculated average photon energy is: E_gamma = {E_gamma:.6e} GeV\n")


# --- 3. Calculate Proton Threshold Energy (E_p_th) ---
# Use the formula for threshold energy in an ultra-relativistic case:
# E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p_th = numerator / denominator

print("--- Step 2: Calculating the proton threshold energy (E_p_th) ---")
print("The formula used is: E_p_th = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print("Equation with final numerical values:")
print(f"  E_p_th = ({m_delta}^2 - {m_p}^2) / (4 * {E_gamma:.6e})")
print(f"  E_p_th = ({m_delta**2:.6f} - {m_p**2:.6f}) / {denominator:.6e}")
print(f"  E_p_th = {numerator:.6f} / {denominator:.6e}\n")


# --- 4. Final Result ---
# Format the final answer in scientific notation with three decimal places.
final_answer = "{:.3e}".format(E_p_th)
print("--- Final Result ---")
print(f"The average threshold energy for the proton is: {final_answer} GeV")


print(f'<<<{final_answer}>>>')