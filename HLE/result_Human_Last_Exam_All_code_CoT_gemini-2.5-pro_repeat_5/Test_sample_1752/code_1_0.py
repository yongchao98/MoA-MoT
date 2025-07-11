import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1.0  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    total_daily_fluid_req = weight_kg * 100
elif weight_kg <= 20:
    total_daily_fluid_req = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_daily_fluid_req = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

print(f"Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar)")
print(f"Based on a weight of {weight_kg} kg, the total daily fluid requirement is {int(total_daily_fluid_req)} ml/day.\n")

# Step 2: Calculate fluid from other sources (drug + feeds)
# Calculate daily drug dose in mg
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
# Calculate daily drug volume in ml
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration
# Calculate total fluid from other sources
total_other_fluids_ml = daily_drug_volume_ml + enteral_feeding_ml_day

print(f"Step 2: Calculate Fluid from Other Sources")
print(f"Daily drug dose = {drug_dose_per_m2_per_day} mg/m^2 * {bsa_m2} m^2 = {int(daily_drug_dose_mg)} mg")
print(f"Daily drug volume = {int(daily_drug_dose_mg)} mg / {admin_concentration} mg/ml = {int(daily_drug_volume_ml)} ml")
print(f"Total fluid from other sources = {int(daily_drug_volume_ml)} ml (drug) + {enteral_feeding_ml_day} ml (feeds) = {int(total_other_fluids_ml)} ml/day.\n")

# Step 3: Calculate the remaining volume for maintenance fluid
maintenance_iv_volume_ml_day = total_daily_fluid_req - total_other_fluids_ml

print(f"Step 3: Calculate Remaining Volume for Maintenance IV Fluid")
print(f"Maintenance fluid volume = {int(total_daily_fluid_req)} ml - {int(total_other_fluids_ml)} ml = {int(maintenance_iv_volume_ml_day)} ml/day.\n")

# Step 4 & 5: Calculate the hourly rate and round to the nearest whole number
maintenance_rate_ml_hr = maintenance_iv_volume_ml_day / 24
rounded_rate = round(maintenance_rate_ml_hr)

print(f"Step 4 & 5: Calculate and Round the Hourly Rate")
print(f"The required rate for the maintenance fluid is {maintenance_iv_volume_ml_day} ml / 24 hours = {maintenance_rate_ml_hr:.1f} ml/hr.")
print(f"Rounded to the nearest whole number, the final rate is {rounded_rate} ml/hr.")

# Final Answer
# This print is for the final answer block and is not part of the step-by-step explanation.
# print(f"<<<{rounded_rate}>>>")
<<<43>>>