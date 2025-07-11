import math

# Step 1: Define the given constants and parameters.

# Beam dimensions
beam_width_focus_mm = 0.3  # mm
beam_height_focus_mm = 6.0 # mm

# Ionization chamber parameters
chamber_current_pA = 2.0  # pA
air_density_mg_cm3 = 1.293 # mg/cm^3
chamber_length_cm = 15.1 # cm

# Physical constants and conversion factors
# Average energy to produce an ion pair in air (standard value)
W_air_J_per_C = 33.97 # J/C

# Timing information for cumulative dose calculation
# Interpreted as the total time a single point is irradiated.
total_exposure_time_s = 0.02 # s

# Unit conversions to SI standard (meters, kilograms, seconds, Coulombs)
# Convert beam dimensions from mm to meters
beam_width_focus_m = beam_width_focus_mm / 1000.0
beam_height_focus_m = beam_height_focus_mm / 1000.0

# Convert chamber length from cm to meters
chamber_length_m = chamber_length_cm / 100.0

# Convert current from pA to A (C/s)
chamber_current_A = chamber_current_pA * 1e-12

# Convert air density from mg/cm^3 to kg/m^3
# 1 mg/cm^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
air_density_kg_m3 = air_density_mg_cm3 # The numeric value is the same

# Step 2: Calculate the mass of the irradiated air in the chamber.
beam_area_m2 = beam_width_focus_m * beam_height_focus_m
irradiated_volume_m3 = beam_area_m2 * chamber_length_m
irradiated_air_mass_kg = irradiated_volume_m3 * air_density_kg_m3

# Step 3: Calculate the dose rate in air (in Gray/second).
# Energy deposited per second (Watts or J/s) = Current (C/s) * W_air (J/C)
energy_per_second_J_s = chamber_current_A * W_air_J_per_C

# Dose rate (Gy/s = J/kg/s) = Energy per second (J/s) / mass (kg)
dose_rate_Gy_s = energy_per_second_J_s / irradiated_air_mass_kg

# As per the problem statement, dose rate in tissue is equal to dose rate in air.
dose_rate_tissue_Gy_s = dose_rate_Gy_s

# Step 4: Calculate the cumulative surface dose.
# Cumulative Dose (Gy) = Dose Rate (Gy/s) * Total Exposure Time (s)
cumulative_dose_Gy = dose_rate_tissue_Gy_s * total_exposure_time_s

# Convert the final dose to milligrays (mGy) for readability.
cumulative_dose_mGy = cumulative_dose_Gy * 1000.0

# Print the final equation and the result
dose_rate_mGy_s = dose_rate_tissue_Gy_s * 1000
print("Calculation Breakdown:")
print(f"Irradiated Air Mass (kg): {irradiated_air_mass_kg:.3e}")
print(f"Dose Rate (mGy/s): {dose_rate_mGy_s:.4f}")
print("\nFinal Cumulative Dose Equation:")
print(f"Cumulative Dose (mGy) = Dose Rate (mGy/s) * Exposure Time (s)")
print(f"Cumulative Dose (mGy) = {dose_rate_mGy_s:.4f} mGy/s * {total_exposure_time_s} s")
print(f"\nCalculated Cumulative Surface Dose: {cumulative_dose_mGy:.4f} mGy")

# Final answer in the requested format
final_answer = round(cumulative_dose_mGy, 4)