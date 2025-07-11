import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25 # mg/m^2/day
admin_concentration_mg_per_ml = 1.0 # mg/ml

# Other fluid sources
milk_volume_per_day = 500 # ml/day

# Step 1: Calculate the total daily dose of the drug
total_daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"Step 1: The total daily drug dose is {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_drug_dose_mg} mg/day.")

# Step 2: Calculate the daily volume of the drug fluid
drug_volume_per_day = total_daily_drug_dose_mg / admin_concentration_mg_per_ml
print(f"Step 2: The daily volume of the drug fluid is {total_daily_drug_dose_mg} mg / {admin_concentration_mg_per_ml} mg/ml = {int(drug_volume_per_day)} ml/day.")

# Step 3: Calculate total daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    total_maintenance_fluid_per_day = weight_kg * 100
elif weight_kg <= 20:
    total_maintenance_fluid_per_day = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_maintenance_fluid_per_day = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
print(f"Step 3: The total daily maintenance fluid requirement for a {weight_kg} kg child is {int(total_maintenance_fluid_per_day)} ml/day.")

# Step 4: Calculate the required volume of IV maintenance fluid
iv_maintenance_fluid_per_day = total_maintenance_fluid_per_day - drug_volume_per_day - milk_volume_per_day
print(f"Step 4: The required IV fluid is {int(total_maintenance_fluid_per_day)} ml - {int(drug_volume_per_day)} ml (drug) - {milk_volume_per_day} ml (milk) = {int(iv_maintenance_fluid_per_day)} ml/day.")

# Step 5 & 6: Calculate the hourly rate and round to the nearest whole number
iv_rate_ml_per_hr = iv_maintenance_fluid_per_day / 24
# Standard rounding (round .5 up)
rounded_iv_rate = int(iv_rate_ml_per_hr + 0.5)

print("\nFinal Calculation:")
# Print the final equation with all the numbers
print(f"The required maintenance fluid rate is calculated as: \n({int(total_maintenance_fluid_per_day)} - {int(drug_volume_per_day)} - {milk_volume_per_day}) / 24 = {rounded_iv_rate} ml/hr")

<<<43>>>