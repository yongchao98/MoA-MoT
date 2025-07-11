import math

# Patient and prescription data
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration_mg_ml = 1  # mg/ml
enteral_feeding_ml_day = 500

print("Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar Method)")
# For a 22 kg child:
# 10 kg * 100 ml/kg = 1000 ml
# 10 kg * 50 ml/kg = 500 ml
# 2 kg * 20 ml/kg = 40 ml
# Total = 1000 + 500 + 40 = 1540 ml/day
first_10_kg_fluid = 10 * 100
next_10_kg_fluid = 10 * 50
remaining_weight_fluid = (weight_kg - 20) * 20
total_daily_maintenance_fluid_ml = first_10_kg_fluid + next_10_kg_fluid + remaining_weight_fluid
print(f"For a {weight_kg} kg child:")
print(f"  Fluid for first 10 kg: 10 * 100 = {first_10_kg_fluid} ml")
print(f"  Fluid for 11-20 kg: 10 * 50 = {next_10_kg_fluid} ml")
print(f"  Fluid for remaining {weight_kg - 20} kg: {weight_kg - 20} * 20 = {remaining_weight_fluid} ml")
print(f"Total Daily Fluid Requirement: {total_daily_maintenance_fluid_ml} ml/day\n")

print("Step 2: Calculate Daily Chemotherapy Volume")
# Calculate total daily dose in mg
total_daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"  Daily Drug Dose (mg): {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_drug_dose_mg} mg")
# Calculate volume of this dose
drug_volume_ml_day = total_daily_drug_dose_mg / admin_concentration_mg_ml
print(f"  Daily Drug Volume (ml): {total_daily_drug_dose_mg} mg / {admin_concentration_mg_ml} mg/ml = {drug_volume_ml_day} ml\n")

print("Step 3: Calculate Total Intake from Other Sources")
total_other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day
print(f"  Chemotherapy Volume + Enteral Feed Volume = {drug_volume_ml_day} ml + {enteral_feeding_ml_day} ml = {total_other_fluids_ml_day} ml/day\n")

print("Step 4: Determine Remaining Fluid for IV Maintenance")
remaining_maintenance_volume_ml_day = total_daily_maintenance_fluid_ml - total_other_fluids_ml_day
print(f"  Total Requirement - Other Sources = {total_daily_maintenance_fluid_ml} ml - {total_other_fluids_ml_day} ml = {remaining_maintenance_volume_ml_day} ml/day\n")

print("Step 5: Calculate the Final Infusion Rate")
maintenance_rate_ml_hr = remaining_maintenance_volume_ml_day / 24
final_rate_rounded = round(maintenance_rate_ml_hr)
print(f"The final calculation is:")
print(f"({total_daily_maintenance_fluid_ml} ml - {drug_volume_ml_day} ml - {enteral_feeding_ml_day} ml) / 24 hours = {final_rate_rounded} ml/hr")
<<<43>>>