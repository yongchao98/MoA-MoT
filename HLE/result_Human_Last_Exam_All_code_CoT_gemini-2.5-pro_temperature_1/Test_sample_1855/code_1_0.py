import math

# Define the physical constants, rounded to three decimal places as requested.
# Masses are in GeV.
m_delta = 1.233  # Delta baryon mass
m_proton = 0.938  # Proton mass

# Boltzmann constant in GeV/K. The standard value is ~8.617333e-14 GeV/K.
k_boltzmann = 8.617e-14

# Define temperature parameters in Kelvin.
T = 2.73  # Current CMB temperature
T0 = 2.7   # Base CMB temperature

# Define coefficients from the energy model.
prop_coeff = 2.71  # Proportionality coefficient
pert_coeff = 1.0e-10  # Perturbative coefficient in GeV*K^-2

# --- Step 1: Calculate the mean energy of a CMB photon (E_gamma) ---
print("Step 1: Calculate the mean photon energy E_gamma.")
delta_T = T - T0
e_gamma_linear = prop_coeff * k_boltzmann * T
e_gamma_pert = pert_coeff * (delta_T**2)
e_gamma = e_gamma_linear + e_gamma_pert

print("The formula for the mean photon energy is: E_gamma = 2.71 * k * T + C * (T - T0)^2")
print("Plugging in the numbers:")
print(f"E_gamma = {prop_coeff} * {k_boltzmann:.3e} * {T} + {pert_coeff:.1e} * ({T} - {T0})^2")
print(f"E_gamma = {prop_coeff} * {k_boltzmann:.3e} * {T} + {pert_coeff:.1e} * ({delta_T:.2f})^2")
print(f"E_gamma = {e_gamma_linear:.5e} + {e_gamma_pert:.5e}")
print(f"E_gamma = {e_gamma:.5e} GeV\n")

# --- Step 2: Calculate the threshold energy for the proton (E_p) ---
print("Step 2: Calculate the proton threshold energy E_p.")
print("The threshold energy is given by the formula: E_p = (m_delta^2 - m_proton^2) / (4 * E_gamma)")

# Calculate numerator and denominator separately for clarity.
m_delta_sq = m_delta**2
m_proton_sq = m_proton**2
numerator = m_delta_sq - m_proton_sq
denominator = 4 * e_gamma
e_proton = numerator / denominator

print("Plugging in the numbers:")
print(f"E_p = ({m_delta}^2 - {m_proton}^2) / (4 * {e_gamma:.5e})")
print(f"E_p = ({m_delta_sq:.5f} - {m_proton_sq:.5f}) / (4 * {e_gamma:.5e})")
print(f"E_p = {numerator:.5f} / {denominator:.5e}")
print(f"E_p = {e_proton:.5e} GeV\n")

# --- Step 3: Provide the final answer in the requested format ---
print("Final Answer:")
print("The average threshold energy for the proton, in scientific notation rounded to three decimal places, is:")
print(f"{e_proton:.3e} GeV")

# Final answer in specified format
final_answer = f"{e_proton:.3e}"
# print(f"<<<{final_answer}>>>")