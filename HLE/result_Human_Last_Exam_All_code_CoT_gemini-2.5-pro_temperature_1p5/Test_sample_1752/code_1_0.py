import math

# Patient and Drug Information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily fluid requirement using Holliday-Segar method
daily_fluid_req = 0
if weight_kg <= 10:
    daily_fluid_req = weight_kg * 100
elif weight_kg <= 20:
    daily_fluid_req = (10 * 100) + ((weight_kg - 10) * 50)
else:
    daily_fluid_req = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

print(f"Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar)")
print(f"Based on a weight of {weight_kg} kg, the total daily fluid requirement is: {int(daily_fluid_req)} ml/day")
print("-" * 30)

# Step 2: Calculate the daily volume of the chemotherapy drug
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration

print(f"Step 2: Calculate Daily Chemotherapy Volume")
print(f"Daily Drug Dose = {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {int(daily_drug_dose_mg)} mg/day")
print(f"Daily Drug Volume = {int(daily_drug_dose_mg)} mg / {admin_concentration} mg/ml = {int(daily_drug_volume_ml)} ml/day")
print("-" * 30)

# Step 3: Sum up all other fluid sources
total_other_fluids_ml_day = enteral_feeding_ml_day + daily_drug_volume_ml

print(f"Step 3: Sum Other Fluid Sources")
print(f"Total Other Fluids = {enteral_feeding_ml_day} ml (milk) + {int(daily_drug_volume_ml)} ml (chemo) = {int(total_other_fluids_ml_day)} ml/day")
print("-" * 30)

# Step 4: Calculate the remaining volume for IV maintenance fluid
iv_maintenance_volume_ml_day = daily_fluid_req - total_other_fluids_ml_day

print(f"Step 4: Calculate IV Maintenance Fluid Volume")
print(f"IV Maintenance Volume = {int(daily_fluid_req)} ml/day - {int(total_other_fluids_ml_day)} ml/day = {int(iv_maintenance_volume_ml_day)} ml/day")
print("-" * 30)

# Step 5 & 6: Convert to hourly rate and round
iv_rate_ml_hr = iv_maintenance_volume_ml_day / 24
rounded_iv_rate = round(iv_rate_ml_hr)

print(f"Step 5 & 6: Calculate Hourly Rate and Round")
print(f"IV Maintenance Rate = {int(iv_maintenance_volume_ml_day)} ml / 24 hr = {iv_rate_ml_hr:.1f} ml/hr")
print(f"The rounded rate for the maintenance fluid is: {rounded_iv_rate} ml/hr")

# Final Answer
print(f"\nFinal Equation: round(({int(daily_fluid_req)} - ({enteral_feeding_ml_day} + {int(daily_drug_volume_ml)})) / 24) = {rounded_iv_rate}")
print("<<<43>>>")