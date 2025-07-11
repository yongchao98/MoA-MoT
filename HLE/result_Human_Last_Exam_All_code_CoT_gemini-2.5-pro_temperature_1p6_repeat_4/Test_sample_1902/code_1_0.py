import math

# --- Given Information ---
# Photon beam size at focus (horizontal x vertical)
beam_width_focus_mm = 0.3  # in mm
beam_height_focus_mm = 6.0  # in mm

# Ionization chamber parameters
chamber_length_cm = 15.1  # in cm
ionization_current_pa = 2.0  # in pA
air_density_mg_cm3 = 1.293  # in mg/cm^3

# Total exposure time for a point on the tissue
total_exposure_time_s = 0.02  # in s

# --- Physical Constants ---
# Mean energy to create an ion pair in air (W_air / e)
W_air_over_e_J_C = 33.97  # in J/C

# --- Unit Conversions ---
# Convert all units to a consistent system (cm, kg, s) for calculation
# Convert beam dimensions from mm to cm
beam_width_focus_cm = beam_width_focus_mm / 10.0
beam_height_focus_cm = beam_height_focus_mm / 10.0

# Convert current from pA to A (C/s)
ionization_current_A = ionization_current_pa * 1e-12

# Convert air density from mg/cm^3 to kg/cm^3
air_density_kg_cm3 = air_density_mg_cm3 / 1e6

# --- Step 1: Calculate the Mass of Irradiated Air ---
# Calculate the beam's cross-sectional area at the focus
beam_area_cm2 = beam_width_focus_cm * beam_height_focus_cm

# Calculate the volume of air in the chamber that is irradiated by the beam
irradiated_volume_cm3 = beam_area_cm2 * chamber_length_cm

# Calculate the mass of the irradiated air
mass_of_air_kg = irradiated_volume_cm3 * air_density_kg_cm3

# --- Step 2: Calculate the Dose Rate in Air (and Tissue) ---
# Dose Rate (Gy/s) = (Current [C/s] / Mass [kg]) * (W_air/e [J/C])
# Note: 1 Gy = 1 J/kg
dose_rate_Gy_s = (ionization_current_A / mass_of_air_kg) * W_air_over_e_J_C

# --- Step 3: Calculate the Cumulative Surface Dose ---
# Cumulative Dose (Gy) = Dose Rate (Gy/s) * Total Exposure Time (s)
cumulative_dose_Gy = dose_rate_Gy_s * total_exposure_time_s

# --- Step 4: Display the Results ---
print("--- Calculation of Cumulative Surface Dose ---")
print("\nStep 1: Calculate the mass of irradiated air (m).")
print(f"m = Air Density * Beam Width * Beam Height * Chamber Length")
print(f"m = {air_density_kg_cm3:.4g} kg/cm^3 * {beam_width_focus_cm:.2f} cm * {beam_height_focus_cm:.2f} cm * {chamber_length_cm:.2f} cm")
print(f"m = {mass_of_air_kg:.4g} kg")

print("\nStep 2: Calculate the dose rate (D_rate).")
print(f"D_rate = (Ionization Current / Mass of Air) * (W_air / e)")
print(f"D_rate = ({ionization_current_A:.2g} A / {mass_of_air_kg:.4g} kg) * {W_air_over_e_J_C} J/C")
print(f"D_rate = {dose_rate_Gy_s:.4g} Gy/s")

print("\nStep 3: Calculate the cumulative dose (D_cumulative).")
print(f"D_cumulative = Dose Rate * Total Exposure Time")
print(f"D_cumulative = {dose_rate_Gy_s:.4g} Gy/s * {total_exposure_time_s} s")
print(f"D_cumulative = {cumulative_dose_Gy:.4g} Gy")

# Convert final answer to microGrays (µGy) for readability
cumulative_dose_uGy = cumulative_dose_Gy * 1e6
print("\n--- Final Answer ---")
print(f"The cumulative surface dose is {cumulative_dose_uGy:.3f} µGy.")

# Final equation with all numbers
print("\nFinal Equation:")
print(f"{cumulative_dose_Gy:.4g} Gy = ( ( {ionization_current_A:.2g} A / ( {air_density_kg_cm3:.4g} kg/cm^3 * {beam_width_focus_cm:.2f} cm * {beam_height_focus_cm:.2f} cm * {chamber_length_cm:.2f} cm ) ) * {W_air_over_e_J_C} J/C ) * {total_exposure_time_s} s")
final_answer = cumulative_dose_uGy
# The final answer is wrapped in <<<>>>
# print(f'<<<{final_answer:.3f}>>>')