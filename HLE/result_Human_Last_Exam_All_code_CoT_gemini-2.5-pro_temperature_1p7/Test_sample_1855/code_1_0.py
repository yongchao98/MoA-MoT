import math

# Step 1: Define the physical constants and given values.
# Masses are given in GeV, and rounded to three decimal places as requested.
m_delta = 1.233  # Mass of Delta baryon in GeV
m_p = 0.938      # Mass of proton in GeV

# Temperature values in Kelvin (K)
T = 2.73
T_ref = 2.7

# Boltzmann constant in GeV/K.
# k = 8.617333262... x 10^-14 GeV/K
# We use a precise value for the calculation.
k_boltzmann = 8.617333e-14 # GeV/K

# Proportionality coefficient factor for the main energy term
prop_coeff_factor = 2.71

# Perturbative coefficient in GeV/K^2
pert_coeff = 1.0e-10

# Step 2: Calculate the CMB photon energy (E_gamma).
delta_T = T - T_ref
main_energy_term = prop_coeff_factor * k_boltzmann * T
pert_energy_term = pert_coeff * (delta_T**2)
E_gamma = main_energy_term + pert_energy_term

# Step 3: Calculate the threshold proton energy (E_p) using the derived formula.
# E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)
numerator = m_delta**2 - m_p**2
denominator = 4 * E_gamma
E_p = numerator / denominator

# Step 4: Print the calculation steps and the final answer.
print("The threshold proton energy E_p is calculated using the formula:")
print("E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)\n")

print("The required values are:")
print(f"  m_delta = {m_delta:.3f} GeV")
print(f"  m_p = {m_p:.3f} GeV")
print(f"  T = {T:.2f} K")
print(f"  E_gamma = ({prop_coeff_factor} * {k_boltzmann:.4e} * {T}) + ({pert_coeff:.1e} * ({T} - {T_ref})**2)")
print(f"  E_gamma = {E_gamma:.4e} GeV\n")

print("Plugging the values into the formula:")
print(f"E_p = ({m_delta:.3f}^2 - {m_p:.3f}^2) / (4 * {E_gamma:.4e})")
print(f"E_p = ({numerator:.6f}) / ({denominator:.6e})")
print(f"E_p = {E_p:.4e} GeV\n")

# Format the final answer in scientific notation with three decimal places.
final_answer_str = f"{E_p:.3e}"
print(f"The final threshold energy for the proton, rounded to three decimal places, is: {final_answer_str} GeV")

# Final answer in the required format
final_answer_value = f"{E_p:.3e}".replace('e+','e')
print(f"<<<{final_answer_value}>>>")
