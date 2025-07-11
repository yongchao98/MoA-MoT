import math

# Plan:
# 1. Define physical constants and given parameters with required precision.
# 2. Calculate the average energy of a CMB photon, including the perturbative effect.
# 3. Use the kinematic threshold energy formula for the reaction p + gamma -> Delta.
# 4. Compute the final threshold energy for the proton.
# 5. Print the steps and the final answer in the specified format.

# Step 1: Define constants and parameters
# Masses are in GeV, Temperature in Kelvin.
m_delta = 1.233  # Mass of Delta baryon in GeV
m_proton = 0.938  # Mass of proton in GeV (rounded to 3 decimal places)

# Boltzmann's constant in eV/K, rounded to 3 decimal places, then converted to GeV/K
k_boltzmann_eV_K = 8.617e-5
k_boltzmann_GeV_K = k_boltzmann_eV_K * 1e-9

T = 2.73  # Temperature of the CMB in K
T_base = 2.7 # Base temperature for perturbation in K

# Coefficients for the photon energy calculation
coeff_proportionality = 2.71
coeff_perturbative = 1.0e-10  # GeV * K^-2

print("--- Step 1: Defined Constants and Parameters ---")
print(f"Mass of Delta baryon (m_delta): {m_delta} GeV")
print(f"Mass of proton (m_proton): {m_proton} GeV")
print(f"Boltzmann's constant (k): {k_boltzmann_GeV_K:.3e} GeV/K")
print(f"CMB Temperature (T): {T} K")
print("-" * 50)

# Step 2: Calculate the average CMB photon energy (E_gamma)
delta_T = T - T_base
main_energy_term = coeff_proportionality * k_boltzmann_GeV_K * T
perturbative_term = coeff_perturbative * (delta_T**2)
E_gamma_mean = main_energy_term + perturbative_term

print("--- Step 2: Calculating Average CMB Photon Energy (E_gamma) ---")
print("Formula: E_gamma = (2.71 * k * T) + (1e-10 * (T - 2.7)^2)")
print(f"Plugging in values: E_gamma = (2.71 * {k_boltzmann_GeV_K:.3e} * {T}) + ({coeff_perturbative:.1e} * ({T} - {T_base})^2)")
print(f"Resulting E_gamma: {E_gamma_mean:.6e} GeV")
print("-" * 50)

# Step 3: Calculate the proton threshold energy (E_p)
m_delta_sq = m_delta**2
m_proton_sq = m_proton**2
mass_diff_sq = m_delta_sq - m_proton_sq

E_proton_threshold = mass_diff_sq / (4 * E_gamma_mean)

print("--- Step 3: Calculating Proton Threshold Energy (E_p) ---")
print("Formula: E_p = (m_delta^2 - m_proton^2) / (4 * E_gamma)")
print(f"Plugging in values: E_p = ({m_delta}^2 - {m_proton}^2) / (4 * {E_gamma_mean:.6e})")
print(f"E_p = ({m_delta_sq:.6f} - {m_proton_sq:.6f}) / {4 * E_gamma_mean:.6e}")
print(f"E_p = {mass_diff_sq:.6f} / {4 * E_gamma_mean:.6e}")
print(f"Resulting E_p: {E_proton_threshold:.6e} GeV")
print("-" * 50)


# Step 4: Final Answer Formatting
final_answer_sci = f"{E_proton_threshold:.3e}"

print("--- Final Answer ---")
print("The average threshold energy for the proton is:")
print(f"{final_answer_sci} GeV")

# Final answer in the specified format
print(f"\n<<<{final_answer_sci}>>>")
