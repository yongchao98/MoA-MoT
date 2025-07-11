import math

# --- Given Information & Constants ---

# Ionization chamber current in Amperes (C/s)
I_pA = 2.0  # pA
I_A = I_pA * 1e-12  # C/s

# Photon beam size at focus (horizontal x vertical) in cm
beam_h_mm = 0.3  # mm
beam_v_mm = 6.0  # mm
beam_h_cm = beam_h_mm / 10.0
beam_v_cm = beam_v_mm / 10.0

# Length of the ionization chamber in cm
L_cm = 15.1

# Density of air in the ionization chamber in kg/cm^3
rho_air_mg_cm3 = 1.293 # mg/cm^3
rho_air_kg_cm3 = rho_air_mg_cm3 * 1e-6 # kg/cm^3

# Exposure time per scan point in seconds
t_exp_s = 0.02

# Standard physical constant: Average energy to create an ion pair in air (W/e)
W_air_J_per_C = 33.97  # J/C

# The assumption is D_tissue = D_air, so no specific tissue coefficient is needed.

# --- Step 1: Calculate the Dose Rate in Air ---

print("Step 1: Calculating the Dose Rate in Air (D_rate)...")

# Calculate the beam's cross-sectional area in cm^2
A_beam_cm2 = beam_h_cm * beam_v_cm
print(f"  - Photon beam area at focus (A_beam) = {beam_h_cm:.3f} cm * {beam_v_cm:.1f} cm = {A_beam_cm2:.4f} cm^2")

# Calculate the volume of air being irradiated in the chamber in cm^3
V_air_cm3 = A_beam_cm2 * L_cm
print(f"  - Irradiated air volume (V_air) = {A_beam_cm2:.4f} cm^2 * {L_cm:.1f} cm = {V_air_cm3:.4f} cm^3")

# Calculate the mass of the irradiated air in kg
m_air_kg = V_air_cm3 * rho_air_kg_cm3
print(f"  - Mass of irradiated air (m_air) = {V_air_cm3:.4f} cm^3 * {rho_air_kg_cm3:.4e} kg/cm^3 = {m_air_kg:.4e} kg")

# Calculate the energy absorbed per second (Power) in J/s
power_J_s = I_A * W_air_J_per_C
print(f"  - Energy absorbed per second (Power) = {I_A:.1e} C/s * {W_air_J_per_C:.2f} J/C = {power_J_s:.4e} J/s")

# Calculate the dose rate in Gy/s (J/kg/s)
# Equation: D_rate = Power / m_air = (I * W_air) / (A_beam * L * rho_air)
dose_rate_Gy_s = power_J_s / m_air_kg
print(f"\nEquation for Dose Rate:")
print(f"D_rate = (Current * W_air) / (Beam Area * Chamber Length * Air Density)")
print(f"D_rate = ({I_A:.1e} C/s * {W_air_J_per_C:.2f} J/C) / ({A_beam_cm2:.4f} cm^2 * {L_cm:.1f} cm * {rho_air_kg_cm3:.4e} kg/cm^3)")
print(f"D_rate = {power_J_s:.4e} J/s / {m_air_kg:.4e} kg")
print(f"D_rate = {dose_rate_Gy_s:.4e} Gy/s\n")

# --- Step 2: Calculate the Cumulative Dose ---

print("Step 2: Calculating the Cumulative Surface Dose (D_cumulative)...")
# Equation: D_cumulative = D_rate * t_exp
cumulative_dose_Gy = dose_rate_Gy_s * t_exp_s
cumulative_dose_mGy = cumulative_dose_Gy * 1000

print(f"Equation for Cumulative Dose:")
print(f"D_cumulative = Dose Rate * Exposure Time")
print(f"D_cumulative = {dose_rate_Gy_s:.4e} Gy/s * {t_exp_s:.2f} s")
print(f"D_cumulative = {cumulative_dose_Gy:.4e} Gy\n")


# --- Final Answer ---

print("="*40)
print("Final Answer:")
print(f"The cumulative surface dose to the tissue is {cumulative_dose_mGy:.5f} mGy.")
print("="*40)

<<<0.00387>>>