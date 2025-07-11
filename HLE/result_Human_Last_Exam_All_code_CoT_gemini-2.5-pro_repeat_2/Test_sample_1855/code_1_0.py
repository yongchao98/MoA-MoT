import math

# --- 1. Define Constants ---
# Given values from the problem
m_delta = 1.233  # Mass of Delta baryon in GeV/c^2
T = 2.73         # CMB temperature in Kelvin
T0 = 2.7         # Base temperature for perturbation in Kelvin
C_prop = 2.71    # Proportionality coefficient for the main energy term
C_pert = 1.0e-10 # Perturbative coefficient in GeV*K^-2

# Physical constants, rounded to three decimal places as instructed
# Proton mass in GeV/c^2 (PDG value is ~0.938272 GeV/c^2)
m_p = 0.938
# Boltzmann's constant in GeV/K (PDG value is ~8.617333e-14 GeV/K)
k = 8.617e-14

# --- 2. Calculate CMB Photon Energy (E_gamma) ---
delta_T = T - T0
# Calculate the two components of the energy
main_energy_term = C_prop * k * T
pert_energy_term = C_pert * (delta_T**2)
E_gamma = main_energy_term + pert_energy_term

print("Step 1: Calculate the mean CMB photon energy (E_gamma)")
print(f"E_gamma = (C_prop * k * T) + (C_pert * (T - T0)^2)")
print(f"E_gamma = ({C_prop} * {k:.3e} GeV/K * {T} K) + ({C_pert:.1e} GeV/K^2 * ({T} K - {T0} K)^2)")
print(f"E_gamma = {main_energy_term} GeV + {pert_energy_term} GeV")
print(f"E_gamma = {E_gamma} GeV\n")


# --- 3. Calculate Proton Threshold Energy (Ep_th) ---
# Calculate the numerator (difference of squared masses)
m_delta_sq = m_delta**2
m_p_sq = m_p**2
numerator = m_delta_sq - m_p_sq

# Calculate the denominator
denominator = 4 * E_gamma

# Calculate the threshold energy
Ep_th = numerator / denominator

print("Step 2: Calculate the proton threshold energy (E_p)")
print(f"E_p = (m_delta^2 - m_p^2) / (4 * E_gamma)")
print(f"E_p = (({m_delta} GeV)^2 - ({m_p} GeV)^2) / (4 * {E_gamma} GeV)")
print(f"E_p = ({m_delta_sq} GeV^2 - {m_p_sq} GeV^2) / ({denominator} GeV)")
print(f"E_p = {numerator} GeV^2 / {denominator} GeV")
print(f"E_p = {Ep_th} GeV\n")


# --- 4. Format and Print Final Answer ---
# Format the final answer to scientific notation with three decimal places
final_answer_formatted = f"{Ep_th:.3e}"

print("Final Answer:")
print(f"The average threshold energy for the proton is {final_answer_formatted} GeV.")

# Final answer in the required format
print(f"<<<{final_answer_formatted}>>>")