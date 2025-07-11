import math

# --- Given Information ---
# Ionization chamber current in Amperes (C/s)
I_pA = 2.0  # pA
I = I_pA * 1e-12  # A or C/s

# Photon beam size at focus (horizontal x vertical) in meters
beam_h_mm = 0.3  # mm
beam_v_mm = 6.0  # mm
beam_h = beam_h_mm * 1e-3  # m
beam_v = beam_v_mm * 1e-3  # m

# Density of air in the ionization chamber in kg/m^3
rho_air_mg_cm3 = 1.293  # mg/cm^3
# 1 mg/cm^3 = 1 kg/m^3
rho_air = rho_air_mg_cm3  # kg/m^3

# Length of ionization chamber in meters
L_cm = 15.1  # cm
L = L_cm * 1e-2  # m

# Effective exposure time from the scan parameters in seconds
t_effective = 0.02  # s

# --- Physical Constants ---
# Mean energy required to produce an ion pair in air (J/C)
W_air = 33.97  # J/C

# --- Step 1: Calculate the Dose Rate ---
# Calculate the cross-sectional area of the beam at the focus
A_beam = beam_h * beam_v

# Calculate the volume of air in the chamber irradiated by the beam
V_air = A_beam * L

# Calculate the mass of the irradiated air
m_air = rho_air * V_air

# Calculate the rate of energy deposition in the air
energy_deposition_rate = I * W_air

# Calculate the dose rate (Energy per mass per time) in Gy/s
# 1 Gy = 1 J/kg
dose_rate = energy_deposition_rate / m_air

# --- Step 2: Calculate the Cumulative Dose ---
# Cumulative dose is the dose rate multiplied by the effective exposure time
cumulative_dose_Gy = dose_rate * t_effective

# Convert the dose from Grays (Gy) to milliGrays (mGy)
cumulative_dose_mGy = cumulative_dose_Gy * 1000

# --- Print the results and the final equation ---
print("--- Step 1: Dose Rate Calculation ---")
print(f"The dose rate is calculated as: (I * W_air) / (rho_air * A_beam * L)")
print(f"Dose Rate = ({I:.1e} C/s * {W_air} J/C) / ({rho_air} kg/m^3 * {A_beam:.1e} m^2 * {L} m)")
print(f"Dose Rate = {energy_deposition_rate:.4e} J/s / {m_air:.4e} kg")
print(f"Dose Rate = {dose_rate:.4e} Gy/s\n")

print("--- Step 2: Cumulative Dose Calculation ---")
print(f"The cumulative dose is calculated as: Dose Rate * Effective Exposure Time")
print(f"Cumulative Dose = {dose_rate:.4e} Gy/s * {t_effective} s")
print(f"Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")
print(f"Cumulative Dose = {cumulative_dose_mGy:.5f} mGy\n")

print("Final Answer:")
print(f"The cumulative surface dose is {cumulative_dose_mGy:.3f} mGy.")

# The final answer in the required format
final_answer = round(cumulative_dose_mGy, 3)
print(f'<<<{final_answer}>>>')