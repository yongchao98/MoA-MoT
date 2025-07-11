import math

# --- Given Information ---
# Beam size at focus (horizontal x vertical) in mm
beam_width_focus_mm = 0.3
beam_height_focus_mm = 6.0

# Ionization chamber current in picoamperes (pA)
current_pa = 2.0

# Density of air in mg/cm^3
density_air_mg_cm3 = 1.293

# Length of ionization chamber in cm
chamber_length_cm = 15.1

# Total exposure time for a single point in seconds
exposure_time_s = 0.02

# --- Physical Constants ---
# Mean energy to create an ion pair in air (W/e) in Joules/Coulomb
W_per_e_air_JC = 33.97

# --- Step 1: Calculate the area and volume of irradiated air ---
# Convert beam dimensions from mm to cm
beam_width_focus_cm = beam_width_focus_mm / 10.0
beam_height_focus_cm = beam_height_focus_mm / 10.0

# Calculate beam area in cm^2
beam_area_cm2 = beam_width_focus_cm * beam_height_focus_cm

# Calculate the volume of irradiated air in cm^3
air_volume_cm3 = beam_area_cm2 * chamber_length_cm

# --- Step 2: Calculate the mass of irradiated air ---
# Convert air density from mg/cm^3 to kg/cm^3
# 1 mg = 1e-6 kg
density_air_kg_cm3 = density_air_mg_cm3 * 1e-6

# Calculate the mass of the irradiated air in kg
air_mass_kg = air_volume_cm3 * density_air_kg_cm3

# --- Step 3: Calculate the absorbed dose rate in air ---
# Convert current from pA to A
# 1 pA = 1e-12 A (which is 1e-12 C/s)
current_A = current_pa * 1e-12

# Calculate energy deposited per second (Power) in J/s (Watts)
energy_per_second_Js = current_A * W_per_e_air_JC

# Calculate dose rate in Gy/s (J/kg/s)
dose_rate_gy_s = energy_per_second_Js / air_mass_kg

# --- Step 4 & 5: Calculate the cumulative dose ---
# Dose to tissue is assumed to be the same as dose to air.
# Cumulative dose = Dose Rate * Exposure Time
cumulative_dose_gy = dose_rate_gy_s * exposure_time_s

# --- Print the results and the final equation ---
print("--- Calculation Steps ---")
print(f"1. Mass of irradiated air (m):")
print(f"   Beam area = {beam_width_focus_cm:.3f} cm * {beam_height_focus_cm:.3f} cm = {beam_area_cm2:.4f} cm^2")
print(f"   Air volume = {beam_area_cm2:.4f} cm^2 * {chamber_length_cm:.1f} cm = {air_volume_cm3:.4f} cm^3")
print(f"   Air mass = {air_volume_cm3:.4f} cm^3 * {density_air_kg_cm3:.4e} kg/cm^3 = {air_mass_kg:.4e} kg")
print("\n2. Dose Rate (D_rate):")
print(f"   Energy per second = Ionization Current (I) * (W/e)")
print(f"   Energy per second = {current_A:.1e} C/s * {W_per_e_air_JC:.2f} J/C = {energy_per_second_Js:.4e} J/s")
print(f"   Dose Rate = Energy per second / Air mass")
print(f"   Dose Rate = {energy_per_second_Js:.4e} J/s / {air_mass_kg:.4e} kg = {dose_rate_gy_s:.4e} Gy/s")
print("\n3. Cumulative Dose (D_cumulative):")
print(f"   Cumulative Dose = Dose Rate * Exposure Time (t)")
print(f"   Final Equation: D_cumulative = ({current_A:.1e} C/s * {W_per_e_air_JC:.2f} J/C) / ({air_volume_cm3:.4f} cm^3 * {density_air_kg_cm3:.4e} kg/cm^3) * {exposure_time_s:.2f} s")
print(f"   Cumulative Dose = {dose_rate_gy_s:.4e} Gy/s * {exposure_time_s:.2f} s = {cumulative_dose_gy:.4e} Gy")

print("\n--- Final Answer ---")
print(f"The cumulative surface dose to the tissue is {cumulative_dose_gy:.3e} Gy.")

# The final answer in the required format.
# Converting Gy to microGrays (µGy) for better readability
cumulative_dose_uGy = cumulative_dose_gy * 1e6
# Report to two significant figures based on the input values (2.0 pA, 0.02 s)
final_answer_str = f"{cumulative_dose_uGy:.1f} uGy" # (or {cumulative_dose_gy:.2e} Gy)
final_value_for_format = f"{cumulative_dose_gy:.2e}"

# The problem asks for the answer directly.
# The calculation gives ~3.86e-6 Gy
# Printing it in scientific notation with two significant figures.
final_value = "{:.2e}".format(cumulative_dose_gy)
# Let's provide the numerical value as requested
final_numerical_value = float(final_value)

<<<3.87e-06>>>