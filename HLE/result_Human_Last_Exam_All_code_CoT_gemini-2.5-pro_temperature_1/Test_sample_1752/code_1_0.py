import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_mg_per_m2 = 25
drug_admin_conc_mg_per_ml = 1
milk_volume_ml_day = 500
hours_in_day = 24

# --- Step 1: Calculate total daily maintenance fluid using Holliday-Segar method ---
if weight_kg <= 10:
    daily_fluid_maintenance = weight_kg * 100
    holliday_segar_calc_str = f"({weight_kg} * 100)"
elif weight_kg <= 20:
    daily_fluid_maintenance = (10 * 100) + ((weight_kg - 10) * 50)
    holliday_segar_calc_str = f"(10 * 100 + ({weight_kg} - 10) * 50)"
else:
    daily_fluid_maintenance = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    holliday_segar_calc_str = f"(10 * 100 + 10 * 50 + ({weight_kg} - 20) * 20)"

print(f"Step 1: Calculate total daily maintenance fluid requirement for a {weight_kg} kg child.")
print(f"Holliday-Segar Calculation: {holliday_segar_calc_str} = {daily_fluid_maintenance} ml/day\n")

# --- Step 2: Calculate fluid from other sources (drug + milk) ---
# Daily drug dose in mg
daily_drug_dose_mg = drug_dose_mg_per_m2 * bsa_m2
# Daily drug volume in ml
daily_drug_volume_ml = daily_drug_dose_mg / drug_admin_conc_mg_per_ml

total_other_fluids_ml = daily_drug_volume_ml + milk_volume_ml_day

print("Step 2: Calculate total fluid from other sources.")
print(f"Daily drug volume = (Dose {drug_dose_mg_per_m2} mg/m² * BSA {bsa_m2} m²) / Concentration {drug_admin_conc_mg_per_ml} mg/ml = {daily_drug_volume_ml} ml/day")
print(f"Enteral feeding volume = {milk_volume_ml_day} ml/day")
print(f"Total other fluids = {daily_drug_volume_ml} ml + {milk_volume_ml_day} ml = {total_other_fluids_ml} ml/day\n")

# --- Step 3: Calculate remaining fluid to be given as maintenance ---
remaining_fluid_ml = daily_fluid_maintenance - total_other_fluids_ml

print("Step 3: Calculate the remaining volume for maintenance fluid.")
print(f"Remaining volume = Total maintenance {daily_fluid_maintenance} ml - Other fluids {total_other_fluids_ml} ml = {remaining_fluid_ml} ml/day\n")

# --- Step 4: Calculate the hourly rate and round ---
maintenance_rate_ml_hr = remaining_fluid_ml / hours_in_day
rounded_rate = round(maintenance_rate_ml_hr)

print("Step 4: Calculate the hourly rate for the maintenance fluid.")
print(f"Rate (ml/hr) = Remaining volume {remaining_fluid_ml} ml / {hours_in_day} hours = {maintenance_rate_ml_hr:.2f} ml/hr")
print(f"Rounding to the nearest whole number gives: {rounded_rate} ml/hr\n")


# --- Final Answer Section ---
print("--- Final Equation ---")
# Build the final equation string with all the numbers
final_equation = f"Rate (ml/hr) = ( (10 * 100 + 10 * 50 + ({weight_kg} - 20) * 20) - (({drug_dose_mg_per_m2} * {bsa_m2}) / {drug_admin_conc_mg_per_ml}) - {milk_volume_ml_day} ) / {hours_in_day}"
print(f"The calculation for the maintenance fluid rate is:\n{final_equation}")
print(f"Result: {rounded_rate}")

# The final answer in the required format
print(f"\n<<<{rounded_rate}>>>")