import math

# --- Given Information ---
# Photon beam size at focus (horizontal by vertical)
beam_width_focus_mm = 0.3  # in mm
beam_height_focus_mm = 6.0   # in mm

# Length of ionization chamber
chamber_length_cm = 15.1  # in cm

# Density of air in the ionization chamber
air_density_mg_cm3 = 1.293  # in mg/cm^3

# Ionization chamber current
ion_current_pA = 2.0  # in pA

# Ratio giving the effective exposure time
exposure_time_s = 0.02  # in s

# --- Physical Constants ---
# Average energy to create an ion pair in air (W/e)
W_air_over_e_J_C = 33.97  # in J/C

# --- Step 1: Calculate Irradiated Volume in SI units (m^3) ---
# Convert dimensions to meters
beam_width_focus_m = beam_width_focus_mm / 1000.0
beam_height_focus_m = beam_height_focus_mm / 1000.0
chamber_length_m = chamber_length_cm / 100.0

# Calculate beam area and irradiated volume
beam_area_m2 = beam_width_focus_m * beam_height_focus_m
irradiated_volume_m3 = beam_area_m2 * chamber_length_m

# --- Step 2: Calculate Mass of Irradiated Air in SI units (kg) ---
# Convert density to kg/m^3 (1 mg/cm^3 = 1 kg/m^3)
air_density_kg_m3 = air_density_mg_cm3

# Calculate mass
mass_air_kg = irradiated_volume_m3 * air_density_kg_m3

# --- Step 3: Calculate Energy Deposited per Second in SI units (J/s) ---
# Convert current to Amperes (C/s)
ion_current_A = ion_current_pA * 1e-12

# Calculate energy deposition rate
energy_per_second_J_s = ion_current_A * W_air_over_e_J_C

# --- Step 4: Calculate Dose Rate in SI units (Gy/s) ---
# Dose Rate = Energy per second / mass
dose_rate_Gy_s = energy_per_second_J_s / mass_air_kg

# --- Step 5: Calculate Cumulative Dose ---
# Cumulative Dose = Dose Rate * Exposure Time
cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

# Convert to microGrays (uGy) for easier interpretation
cumulative_dose_uGy = cumulative_dose_Gy * 1e6

# --- Final Output ---
print("This script calculates the cumulative surface dose based on the provided parameters.")
print("\n--- Calculation Steps ---")
print(f"1. Dose Rate Calculation:")
print(f"   - Energy deposited per second: {energy_per_second_J_s:.3e} J/s")
print(f"   - Mass of irradiated air: {mass_air_kg:.3e} kg")
print(f"   - Dose Rate = ({energy_per_second_J_s:.3e} J/s) / ({mass_air_kg:.3e} kg) = {dose_rate_Gy_s:.4f} Gy/s")

print(f"\n2. Cumulative Dose Calculation:")
print(f"   - The cumulative dose is the product of the dose rate and the exposure time.")
print(f"   - Final Equation: Dose Rate ({dose_rate_Gy_s:.4f} Gy/s) * Exposure Time ({exposure_time_s} s) = Cumulative Dose ({cumulative_dose_Gy:.3e} Gy)")

print(f"\nThe cumulative surface dose is {cumulative_dose_uGy:.2f} \u00b5Gy.")

# The final answer format requires just the number.
# Using the value in microGrays, rounded to two decimal places.
final_answer = round(cumulative_dose_uGy, 2)
# print(f'<<< {final_answer} >>>') # This is for thinking, not final output
# The instruction wants just the <<<answer>>> tag at the end of the response.
# The value is 3.87 uGy.