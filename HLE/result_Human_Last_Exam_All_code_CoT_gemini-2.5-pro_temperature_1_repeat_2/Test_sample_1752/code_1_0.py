import math

# Patient and prescription details
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate Total Daily Maintenance Fluid (Holliday-Segar method)
if weight_kg <= 10:
    total_maintenance_fluid_ml_day = weight_kg * 100
elif weight_kg <= 20:
    total_maintenance_fluid_ml_day = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_maintenance_fluid_ml_day = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

print(f"Step 1: Calculate Total Daily Maintenance Fluid")
print(f"Based on the Holliday-Segar method for a {weight_kg} kg child:")
print(f"10 kg * 100 ml/kg + 10 kg * 50 ml/kg + {weight_kg - 20} kg * 20 ml/kg = {int(total_maintenance_fluid_ml_day)} ml/day")
print("-" * 30)

# Step 2: Calculate Daily Drug Volume
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration

print(f"Step 2: Calculate Daily Drug Volume")
print(f"Daily drug dose: {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg/day")
print(f"Daily drug volume: {daily_drug_dose_mg} mg / {admin_concentration} mg/ml = {daily_drug_volume_ml} ml/day")
print("-" * 30)

# Step 3: Calculate Total Fluid from Other Sources
total_other_fluids_ml_day = daily_drug_volume_ml + enteral_feeding_ml_day

print(f"Step 3: Calculate Total Fluid from Other Sources")
print(f"Total from other sources: {daily_drug_volume_ml} ml (drug) + {enteral_feeding_ml_day} ml (milk) = {total_other_fluids_ml_day} ml/day")
print("-" * 30)

# Step 4: Determine Required Maintenance Fluid Volume
remaining_maintenance_fluid_ml_day = total_maintenance_fluid_ml_day - total_other_fluids_ml_day

print(f"Step 4: Determine Required Maintenance Fluid Volume")
print(f"Required maintenance fluid: {int(total_maintenance_fluid_ml_day)} ml (total) - {total_other_fluids_ml_day} ml (other) = {remaining_maintenance_fluid_ml_day} ml/day")
print("-" * 30)

# Step 5: Calculate the Final Hourly Rate
maintenance_fluid_rate_ml_hr = remaining_maintenance_fluid_ml_day / 24
rounded_rate = round(maintenance_fluid_rate_ml_hr)

print(f"Step 5: Calculate the Final Hourly Rate")
print(f"Hourly rate: {remaining_maintenance_fluid_ml_day} ml / 24 hours = {maintenance_fluid_rate_ml_hr:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the final rate is {rounded_rate} ml/hr.")

# Final Answer
final_answer = rounded_rate
print(f"\nFinal Answer: The rate of the maintenance fluid should be {final_answer} ml/hr.")