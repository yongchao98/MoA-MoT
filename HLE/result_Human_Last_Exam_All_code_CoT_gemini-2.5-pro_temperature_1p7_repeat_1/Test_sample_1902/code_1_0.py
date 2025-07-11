import math

# --- Given Information & Constants ---

# Beam properties
beam_h_mm = 0.3  # mm
beam_v_mm = 6.0  # mm

# Ionization chamber properties
chamber_current_pA = 2.0  # pA
air_density_mg_cm3 = 1.293 # mg/cm^3
chamber_length_cm = 15.1  # cm

# Scan property
effective_exposure_time_s = 0.02 # s

# Physical constant
W_air_J_per_C = 33.97  # Average energy to produce an ion pair in air (J/C)

# --- Unit Conversions for Calculation ---

# Convert beam dimensions from mm to cm
beam_h_cm = beam_h_mm / 10.0
beam_v_cm = beam_v_mm / 10.0

# Convert air density from mg/cm^3 to kg/cm^3
air_density_kg_cm3 = air_density_mg_cm3 / 1000.0 / 1000.0

# Convert chamber current from pA to A (C/s)
chamber_current_A = chamber_current_pA * 1e-12

# --- Step-by-Step Calculation ---

# 1. Calculate the Mass of Irradiated Air
# Volume = cross-sectional area of the beam * length of the chamber
beam_area_cm2 = beam_h_cm * beam_v_cm
irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm

# Mass = density * volume
mass_air_kg = air_density_kg_cm3 * irradiated_volume_cm3

# 2. Calculate the Energy Deposition Rate (Power)
# Power (J/s) = Current (C/s) * W_air (J/C)
energy_deposition_rate_J_s = chamber_current_A * W_air_J_per_C

# 3. Calculate the Dose Rate
# Dose Rate (Gy/s) = Energy Deposition Rate (J/s) / Mass (kg)
# (1 Gy = 1 J/kg)
dose_rate_Gy_s = energy_deposition_rate_J_s / mass_air_kg

# 4. Calculate the Cumulative Dose
# Cumulative Dose (Gy) = Dose Rate (Gy/s) * Effective Exposure Time (s)
cumulative_dose_Gy = dose_rate_Gy_s * effective_exposure_time_s

# --- Output the Results ---
print("--- Calculation Steps ---")
print(f"1. Mass of Irradiated Air:")
print(f"   Beam Area = {beam_h_cm:.3f} cm * {beam_v_cm:.1f} cm = {beam_area_cm2:.4f} cm^2")
print(f"   Irradiated Volume = {beam_area_cm2:.4f} cm^2 * {chamber_length_cm:.1f} cm = {irradiated_volume_cm3:.4f} cm^3")
print(f"   Mass of Air = {air_density_kg_cm3:.4e} kg/cm^3 * {irradiated_volume_cm3:.4f} cm^3 = {mass_air_kg:.4e} kg")
print("\n2. Energy Deposition Rate:")
print(f"   Rate = {chamber_current_A:.2e} C/s * {W_air_J_per_C:.2f} J/C = {energy_deposition_rate_J_s:.4e} J/s")
print("\n3. Dose Rate:")
print(f"   Rate = {energy_deposition_rate_J_s:.4e} J/s / {mass_air_kg:.4e} kg = {dose_rate_Gy_s:.6f} Gy/s")
print("\n4. Final Cumulative Dose Calculation:")
print(f"   Cumulative Dose = Dose Rate * Exposure Time")
print(f"   Cumulative Dose = {dose_rate_Gy_s:.6f} Gy/s * {effective_exposure_time_s:.2f} s")
print(f"   Cumulative Dose = {cumulative_dose_Gy:.8f} Gy")

# For automated checking, output the final answer in scientific notation
# print(f"\nFinal Answer (Scientific Notation): {cumulative_dose_Gy:.3e} Gy")

<<<3.86568e-06>>>