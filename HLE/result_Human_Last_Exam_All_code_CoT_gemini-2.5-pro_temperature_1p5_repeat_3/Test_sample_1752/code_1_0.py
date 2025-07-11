import math

# Patient and Drug Information
weight_kg = 22
recommended_dose_mg_per_m2_day = 25
bsa_m2 = 0.8
admin_concentration_mg_per_ml = 1.0
enteral_feeding_ml_day = 500.0

# --- Step 1: Calculate total daily fluid requirement (Holliday-Segar Method) ---
# For the first 10 kg
fluid_for_first_10kg = 10 * 100
# For the next 10 kg (from 11 to 20 kg)
fluid_for_next_10kg = 10 * 50
# For the remaining weight (above 20 kg)
remaining_weight = weight_kg - 20
fluid_for_remaining_weight = remaining_weight * 20

total_daily_fluid_requirement = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_weight
print("1. Calculate Total Daily Fluid Requirement (Holliday-Segar):")
print(f"For a {weight_kg} kg child:")
print(f"({10} kg * {100} ml/kg) + ({10} kg * {50} ml/kg) + ({remaining_weight} kg * {20} ml/kg) = {total_daily_fluid_requirement} ml/day")
print("-" * 30)

# --- Step 2: Calculate the volume of the daily chemotherapy infusion ---
daily_drug_dose_mg = recommended_dose_mg_per_m2_day * bsa_m2
chemo_fluid_volume_day = daily_drug_dose_mg / admin_concentration_mg_per_ml
print("2. Calculate Fluid Volume from Chemotherapy:")
print(f"Daily Dose = {recommended_dose_mg_per_m2_day} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg")
print(f"Chemo Fluid Volume = {daily_drug_dose_mg} mg / {admin_concentration_mg_per_ml} mg/ml = {chemo_fluid_volume_day} ml/day")
print("-" * 30)

# --- Step 3: Sum all other fluid sources ---
total_other_fluids_day = enteral_feeding_ml_day + chemo_fluid_volume_day
print("3. Calculate Total Fluid from Other Sources:")
print(f"Other Fluids = Enteral Feeding + Chemo Volume = {enteral_feeding_ml_day} ml + {chemo_fluid_volume_day} ml = {total_other_fluids_day} ml/day")
print("-" * 30)

# --- Step 4: Calculate remaining maintenance fluid volume ---
remaining_iv_fluid_day = total_daily_fluid_requirement - total_other_fluids_day
print("4. Calculate Remaining Maintenance Fluid Volume:")
print(f"Remaining Volume = Total Requirement - Other Fluids = {total_daily_fluid_requirement} ml - {total_other_fluids_day} ml = {remaining_iv_fluid_day} ml/day")
print("-" * 30)

# --- Step 5: Calculate and round the hourly rate ---
iv_rate_ml_hr = remaining_iv_fluid_day / 24
final_rate_rounded = round(iv_rate_ml_hr)
print("5. Calculate Final Hourly Rate:")
print(f"The equation to find the final hourly rate is:")
print(f"(({total_daily_fluid_requirement}) - ({enteral_feeding_ml_day} + {chemo_fluid_volume_day})) / 24")
print(f"= ({remaining_iv_fluid_day}) / 24 = {iv_rate_ml_hr:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the maintenance fluid rate is: {final_rate_rounded} ml/hr")

# --- Final Answer ---
final_answer = final_rate_rounded
print(f'<<<{final_answer}>>>')