import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration_mg_per_ml = 1 # mg/ml
enteral_feeding_ml_per_day = 500 # ml/day

# Step 1: Calculate total daily fluid requirement using Holliday-Segar method
# For the first 10 kg: 100 ml/kg
# For the next 10 kg (11-20 kg): 50 ml/kg
# For weight > 20 kg: 20 ml/kg
fluid_for_first_10kg = 10 * 100
fluid_for_next_10kg = 10 * 50
fluid_for_remaining_kg = (weight_kg - 20) * 20
total_daily_fluid_ml = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg

print(f"Step 1: Calculate the child's total daily fluid requirement for a weight of {weight_kg} kg.")
print(f"Calculation: ({10} * {100}) + ({10} * {50}) + (({weight_kg} - {20}) * {20}) = {fluid_for_first_10kg} + {fluid_for_next_10kg} + {fluid_for_remaining_kg} = {total_daily_fluid_ml} ml/day.")
print("-" * 50)

# Step 2: Calculate the total daily volume of the chemotherapy drug
total_drug_dose_mg_per_day = drug_dose_per_m2_per_day * bsa_m2
drug_volume_ml_per_day = total_drug_dose_mg_per_day / admin_concentration_mg_per_ml

print(f"Step 2: Calculate the total daily volume of the chemotherapy drug.")
print(f"Daily drug dose (mg): {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_drug_dose_mg_per_day} mg/day.")
print(f"Daily drug volume (ml): {total_drug_dose_mg_per_day} mg / {admin_concentration_mg_per_ml} mg/ml = {drug_volume_ml_per_day} ml/day.")
print("-" * 50)

# Step 3: Calculate the total fluid from non-maintenance sources
other_fluids_ml_per_day = drug_volume_ml_per_day + enteral_feeding_ml_per_day

print(f"Step 3: Calculate the total fluid from non-maintenance sources (Drug + Milk).")
print(f"Calculation: {drug_volume_ml_per_day} ml/day + {enteral_feeding_ml_per_day} ml/day = {other_fluids_ml_per_day} ml/day.")
print("-" * 50)

# Step 4: Calculate the remaining volume for maintenance fluid
maintenance_fluid_ml_per_day = total_daily_fluid_ml - other_fluids_ml_per_day

print(f"Step 4: Calculate the remaining volume to be given as maintenance fluid.")
print(f"Calculation: {total_daily_fluid_ml} ml/day - {other_fluids_ml_per_day} ml/day = {maintenance_fluid_ml_per_day} ml/day.")
print("-" * 50)

# Step 5: Calculate the hourly rate for the maintenance fluid
maintenance_rate_ml_per_hr = maintenance_fluid_ml_per_day / 24
final_rate = round(maintenance_rate_ml_per_hr)

print(f"Step 5: Calculate the hourly rate and round to the nearest whole number.")
print(f"Calculation: {maintenance_fluid_ml_per_day} ml / 24 hours = {maintenance_rate_ml_per_hr:.1f} ml/hr.")
print(f"The final rate for the maintenance fluid is {final_rate} ml/hr.")
<<<43>>>