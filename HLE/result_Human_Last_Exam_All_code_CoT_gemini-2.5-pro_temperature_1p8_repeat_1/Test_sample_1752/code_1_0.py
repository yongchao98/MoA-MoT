import math

# Define the given variables
weight_kg = 22
drug_dose_per_m2 = 25  # mg/m²/day
bsa_m2 = 0.8
drug_admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate the total daily dose of the chemotherapy drug in mg
daily_drug_dose_mg = drug_dose_per_m2 * bsa_m2

# Step 2: Calculate the daily volume of the drug in ml
daily_drug_volume_ml = daily_drug_dose_mg / drug_admin_concentration

# Step 3: Calculate total daily fluid requirement using Holliday-Segar method
fluid_for_first_10kg = 10 * 100
fluid_for_next_10kg = 10 * 50
fluid_for_remaining_kg = (weight_kg - 20) * 20
total_daily_fluid_ml = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg

# Step 4: Determine total fluid intake from other sources
total_other_fluids_ml = daily_drug_volume_ml + enteral_feeding_ml_day

# Step 5: Calculate the remaining fluid volume to be given as maintenance fluid
remaining_maintenance_fluid_ml_day = total_daily_fluid_ml - total_other_fluids_ml

# Step 6: Calculate the hourly rate and round to the nearest whole number
hourly_rate = remaining_maintenance_fluid_ml_day / 24
rounded_hourly_rate = round(hourly_rate)

# Output the step-by-step calculation
print("Calculation of the Maintenance Fluid Rate:\n")

print(f"Daily Drug Dose: {drug_dose_per_m2} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg")
print(f"Daily Drug Volume: {daily_drug_dose_mg} mg / {drug_admin_concentration} mg/ml = {daily_drug_volume_ml} ml")
print(f"Total Daily Fluid Requirement (Holliday-Segar): {fluid_for_first_10kg} ml + {fluid_for_next_10kg} ml + {fluid_for_remaining_kg} ml = {total_daily_fluid_ml} ml")
print(f"Fluid From Other Sources: {daily_drug_volume_ml} ml (Drug) + {enteral_feeding_ml_day} ml (Milk) = {total_other_fluids_ml} ml")
print(f"Remaining Maintenance Fluid: {total_daily_fluid_ml} ml - {total_other_fluids_ml} ml = {remaining_maintenance_fluid_ml_day} ml")
print(f"Hourly Rate (Unrounded): {remaining_maintenance_fluid_ml_day} ml / 24 hours = {hourly_rate:.1f} ml/hr")
print(f"\nFinal Rounded Hourly Rate: {rounded_hourly_rate} ml/hr")

print("\n--- Final Equation ---")
print(f"({total_daily_fluid_ml} - ({daily_drug_volume_ml} + {enteral_feeding_ml_day})) / 24 = {rounded_hourly_rate}")
