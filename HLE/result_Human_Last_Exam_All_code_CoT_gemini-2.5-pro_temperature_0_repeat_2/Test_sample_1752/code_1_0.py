import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily fluid requirement (Holliday-Segar method)
# For a 22 kg child:
# 10 kg * 100 ml/kg = 1000 ml
# 10 kg * 50 ml/kg  = 500 ml
# 2 kg * 20 ml/kg   = 40 ml
# Total = 1540 ml/day
fluid_for_first_10kg = 10 * 100
fluid_for_next_10kg = 10 * 50
fluid_for_remaining_kg = (weight_kg - 20) * 20
total_fluid_ml_day = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg

print("--- Calculation Steps ---")
print(f"1. Total Daily Fluid Requirement (Holliday-Segar):")
print(f"   ({10} * {100}) + ({10} * {50}) + (({weight_kg} - {20}) * {20}) = {total_fluid_ml_day} ml/day")
print("-" * 25)

# Step 2: Calculate the daily volume of the chemotherapy drug
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
drug_volume_ml_day = daily_drug_dose_mg / admin_concentration

print(f"2. Daily Chemotherapy Fluid Volume:")
print(f"   (({drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m²) / {admin_concentration} mg/ml) = {drug_volume_ml_day} ml/day")
print("-" * 25)

# Step 3: Sum all non-maintenance fluids
total_other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day

print(f"3. Total Fluid from Other Sources (Drug + Milk):")
print(f"   {drug_volume_ml_day} ml + {enteral_feeding_ml_day} ml = {total_other_fluids_ml_day} ml/day")
print("-" * 25)

# Step 4: Calculate the remaining volume for IV maintenance fluid
maintenance_fluid_ml_day = total_fluid_ml_day - total_other_fluids_ml_day

print(f"4. Remaining IV Maintenance Fluid Volume:")
print(f"   {total_fluid_ml_day} ml - {total_other_fluids_ml_day} ml = {maintenance_fluid_ml_day} ml/day")
print("-" * 25)

# Step 5 & 6: Calculate the hourly rate and round
hourly_rate = maintenance_fluid_ml_day / 24
final_rate = round(hourly_rate)

print(f"5. Final Hourly Rate Calculation:")
print(f"   ({maintenance_fluid_ml_day} ml / {24} hours) = {hourly_rate:.1f} ml/hr")
print(f"   Rounding to the nearest whole number gives: {final_rate} ml/hr")
print("-" * 25)

print("\nFinal Answer:")
print(f"The rate of the maintenance fluid should be {final_rate} ml/hr.")
<<<43>>>