import math

# Given parameters from the problem
beam_width_mm = 0.3
beam_height_mm = 6.0
chamber_length_cm = 15.1
chamber_current_pA = 2.0
air_density_mg_cm3 = 1.293
exposure_time_s = 0.02

# Physical constants
W_air_eV_per_ion_pair = 33.97  # Mean energy to create an ion pair in dry air (eV)
e_charge_C = 1.602176634e-19  # Elementary charge in Coulombs
eV_to_J = 1.602176634e-19    # Conversion factor from eV to Joules

# --- Convert all units to SI (meters, kilograms, seconds, Amperes) ---

# Convert current from picoamperes to amperes
current_A = chamber_current_pA * 1e-12

# Convert dimensions from mm/cm to meters
beam_width_m = beam_width_mm / 1000.0
beam_height_m = beam_height_mm / 1000.0
chamber_length_m = chamber_length_cm / 100.0

# Convert air density from mg/cm^3 to kg/m^3
# 1 mg/cm^3 = (1e-6 kg) / (1e-6 m^3) = 1 kg/m^3
air_density_kg_m3 = air_density_mg_cm3

# --- Step-by-step Calculation ---

# 1. Calculate the rate of energy absorption (Power) in J/s
# Power = (Number of ion pairs per second) * (Energy per ion pair in Joules)
# Number of ion pairs per second = current_A / e_charge_C
# Energy per ion pair in Joules = W_air_eV_per_ion_pair * eV_to_J
W_air_J_per_ion_pair = W_air_eV_per_ion_pair * eV_to_J
rate_energy_J_s = (current_A / e_charge_C) * W_air_J_per_ion_pair

# 2. Calculate the mass of the irradiated air in kg
# Volume = width * height * length
irradiated_volume_m3 = beam_width_m * beam_height_m * chamber_length_m
# Mass = Volume * Density
mass_air_kg = irradiated_volume_m3 * air_density_kg_m3

# 3. Calculate the dose rate in Gy/s (1 Gy = 1 J/kg)
dose_rate_Gy_s = rate_energy_J_s / mass_air_kg

# 4. Calculate the cumulative dose in Gy
cumulative_dose_Gy = dose_rate_Gy_s * exposure_time_s

# --- Print the results and the final equation ---
print("Calculation Steps for Cumulative Surface Dose:\n")

print("1. Rate of Energy Absorption (Power):")
print(f"Power (J/s) = (Current ({current_A:.1e} A) / Elementary Charge ({e_charge_C:.6e} C)) * W-value ({W_air_J_per_ion_pair:.4e} J/ion-pair)")
print(f"Power = {rate_energy_J_s:.4e} J/s\n")

print("2. Mass of Irradiated Air:")
print(f"Mass (kg) = (Beam Width ({beam_width_m} m) * Beam Height ({beam_height_m} m) * Chamber Length ({chamber_length_m} m)) * Density ({air_density_kg_m3} kg/m^3)")
print(f"Mass = {mass_air_kg:.4e} kg\n")

print("3. Dose Rate:")
print(f"Dose Rate (Gy/s) = Power ({rate_energy_J_s:.4e} J/s) / Mass ({mass_air_kg:.4e} kg)")
print(f"Dose Rate = {dose_rate_Gy_s:.4e} Gy/s\n")

print("4. Cumulative Surface Dose:")
print(f"Cumulative Dose (Gy) = Dose Rate ({dose_rate_Gy_s:.4e} Gy/s) * Exposure Time ({exposure_time_s} s)")
print(f"Final Cumulative Dose = {cumulative_dose_Gy:.4e} Gy")

# You can also express the answer in microGray (uGy) for easier reading
cumulative_dose_uGy = cumulative_dose_Gy * 1e6
print(f"Final Cumulative Dose = {cumulative_dose_uGy:.3f} microGray (uGy)")
print(f"\nFinal numerical answer in Gray (Gy):")
print(f"{cumulative_dose_Gy}")

<<<3.866465355605051e-06>>>