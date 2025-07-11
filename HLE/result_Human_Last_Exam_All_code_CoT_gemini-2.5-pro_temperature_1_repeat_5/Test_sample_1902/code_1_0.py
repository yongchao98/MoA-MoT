import math

# --- Given Information ---
# Photon beam size at focus (horizontal x vertical)
beam_width_focus_mm = 0.3  # in mm
beam_height_focus_mm = 6.0   # in mm

# Ionization chamber parameters
chamber_current_pA = 2.0  # in picoamperes (pA)
air_density_mg_cm3 = 1.293 # in mg/cm^3
chamber_length_cm = 15.1  # in cm

# Exposure time for a point on the surface
exposure_time_s = 0.02  # in seconds

# --- Physical Constants ---
# Energy to create an ion pair in air (W_air)
W_air_eV_per_ion_pair = 33.97  # in eV

# Charge of a single electron/ion
e_coulombs = 1.602e-19  # in Coulombs

# --- Unit Conversions ---
# Convert all measurements to a consistent system (meters, kilograms, seconds)
# Beam dimensions
beam_width_m = beam_width_focus_mm / 1000.0
beam_height_m = beam_height_focus_mm / 1000.0
# Chamber length
chamber_length_m = chamber_length_cm / 100.0
# Air density
air_density_kg_m3 = air_density_mg_cm3 * (1e-6 / 1e-6) # mg/cm^3 -> kg/m^3
# Chamber current
chamber_current_A = chamber_current_pA * 1e-12  # pA to A (C/s)
# Joules per eV
J_per_eV = 1.602e-19

# Step 1: Calculate the mass of the irradiated air
beam_area_m2 = beam_width_m * beam_height_m
irradiated_volume_m3 = beam_area_m2 * chamber_length_m
irradiated_air_mass_kg = irradiated_volume_m3 * air_density_kg_m3

print(f"Step 1: Calculating the mass of irradiated air")
print(f"Beam area at focus = {beam_width_m:.4f} m * {beam_height_m:.4f} m = {beam_area_m2:.2e} m^2")
print(f"Irradiated volume = {beam_area_m2:.2e} m^2 * {chamber_length_m:.3f} m = {irradiated_volume_m3:.3e} m^3")
print(f"Mass of irradiated air = {irradiated_volume_m3:.3e} m^3 * {air_density_kg_m3:.4f} kg/m^3 = {irradiated_air_mass_kg:.3e} kg\n")

# Step 2: Calculate the energy deposition rate (Power) in Joules/second
# Power (J/s) = (I / e) * W_air * J_per_eV
# Note: I is in C/s, e is in C/ion_pair, W_air is in eV/ion_pair, J_per_eV is in J/eV
# The units combine to J/s.
# A simpler way is Power (J/s) = Current (C/s) * W_air (J/C), where W_air in J/C is W_air in eV/ip * J_per_eV / e_coulombs
# This simplifies to Power (J/s) = Current (A) * W_air (eV/ip)
power_J_per_s = chamber_current_A * W_air_eV_per_ion_pair

print(f"Step 2: Calculating energy deposition rate (Power)")
print(f"Power = {chamber_current_A:.1e} C/s * {W_air_eV_per_ion_pair:.2f} J/C = {power_J_per_s:.3e} J/s\n")

# Step 3 & 4: Calculate the dose rate in tissue (in Gray/second)
# Dose Rate = Power (J/s) / Mass (kg). Result is in Gy/s (since 1 Gy = 1 J/kg)
# We assume dose rate in tissue is the same as in air.
dose_rate_Gy_s = power_J_per_s / irradiated_air_mass_kg

print(f"Step 3 & 4: Calculating dose rate")
print(f"Dose Rate = {power_J_per_s:.3e} J/s / {irradiated_air_mass_kg:.3e} kg = {dose_rate_Gy_s:.4f} Gy/s\n")

# Step 5: Calculate the final cumulative dose
cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s
cumulative_dose_mGy = cumulative_dose_Gy * 1000

print(f"Step 5: Calculating cumulative dose")
print(f"Cumulative Dose (Gy) = Dose Rate (Gy/s) * Exposure Time (s)")
print(f"Cumulative Dose (Gy) = {dose_rate_Gy_s:.4f} * {exposure_time_s:.2f} = {cumulative_dose_Gy:.6f} Gy")
print(f"This is equivalent to {cumulative_dose_mGy:.2f} mGy.")

# Final Answer
# The requested format is to show the numbers in the final equation.
# Here we represent the final calculation for the cumulative dose in Gy.
print("\n--- Final Calculation ---")
print(f"Cumulative Surface Dose = {dose_rate_Gy_s:.4f} Gy/s * {exposure_time_s:.2f} s")
final_answer = dose_rate_Gy_s * exposure_time_s
print(f"Result = {final_answer:.6f} Gy")
<<<0.003865>>>