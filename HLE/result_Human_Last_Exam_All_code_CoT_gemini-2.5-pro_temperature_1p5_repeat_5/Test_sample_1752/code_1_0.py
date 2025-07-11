import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_day = 25
admin_concentration_mg_ml = 1
milk_volume_day_ml = 500

# Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    total_maintenance_fluid_day_ml = weight_kg * 100
elif weight_kg <= 20:
    total_maintenance_fluid_day_ml = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_maintenance_fluid_day_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

print(f"Step 1: Calculate total daily maintenance fluid for a {weight_kg} kg child.")
print(f"Total Daily Maintenance Fluid Requirement: {total_maintenance_fluid_day_ml} ml/day\n")

# Step 2: Calculate the daily chemotherapy infusion volume
daily_drug_dose_mg = bsa_m2 * drug_dose_per_m2_day
chemo_volume_day_ml = daily_drug_dose_mg / admin_concentration_mg_ml

print("Step 2: Calculate the daily chemotherapy infusion volume.")
print(f"Daily Drug Dose: {bsa_m2} m² * {drug_dose_per_m2_day} mg/m²/day = {daily_drug_dose_mg} mg/day")
print(f"Daily Chemo Infusion Volume: {daily_drug_dose_mg} mg / {admin_concentration_mg_ml} mg/ml = {chemo_volume_day_ml} ml/day\n")

# Step 3: Sum up fluids from other sources
total_other_fluids_ml = milk_volume_day_ml + chemo_volume_day_ml
print("Step 3: Calculate total daily fluid from other sources.")
print(f"Total from other sources: {milk_volume_day_ml} ml (milk) + {chemo_volume_day_ml} ml (chemo) = {total_other_fluids_ml} ml/day\n")

# Step 4: Calculate the volume for IV maintenance fluid
maintenance_volume_to_give_ml = total_maintenance_fluid_day_ml - total_other_fluids_ml
print("Step 4: Calculate remaining volume for maintenance fluid.")
print(f"Remaining maintenance volume: {total_maintenance_fluid_day_ml} ml - {total_other_fluids_ml} ml = {maintenance_volume_to_give_ml} ml/day\n")


# Step 5: Calculate the hourly rate
hourly_rate_ml_hr = maintenance_volume_to_give_ml / 24
final_rate_ml_hr = round(hourly_rate_ml_hr)
print("Step 5 & 6: Calculate and round the hourly rate.")
print(f"The final calculation is: ( {total_maintenance_fluid_day_ml} ml - {milk_volume_day_ml} ml - {chemo_volume_day_ml} ml ) / 24 hours = {hourly_rate_ml_hr:.1f} ml/hr")
print(f"The final rate rounded to the nearest whole number is: {final_rate_ml_hr} ml/hr")

<<<43>>>