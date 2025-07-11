import math

# --- Given Information ---
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# --- Step 1: Calculate total daily maintenance fluid using Holliday-Segar method ---
if weight_kg <= 10:
    daily_fluid_maintenance_ml = weight_kg * 100
    calc_1 = f"({weight_kg} kg * 100 ml/kg)"
elif weight_kg <= 20:
    daily_fluid_maintenance_ml = (10 * 100) + ((weight_kg - 10) * 50)
    calc_1 = f"(10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg)"
else:
    daily_fluid_maintenance_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    calc_1 = f"(10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({weight_kg - 20} kg * 20 ml/kg)"

print(f"1. Calculate total daily maintenance fluid requirement:")
print(f"   Calculation: {calc_1} = {daily_fluid_maintenance_ml} ml/day")
print("-" * 30)

# --- Step 2: Calculate the daily volume of the chemotherapy drug ---
total_daily_drug_dose_mg = bsa_m2 * drug_dose_per_m2_per_day
drug_volume_ml_day = total_daily_drug_dose_mg / admin_concentration

print(f"2. Calculate the daily volume of the chemotherapy drug:")
print(f"   Total daily drug dose = {bsa_m2} m² * {drug_dose_per_m2_per_day} mg/m² = {total_daily_drug_dose_mg} mg/day")
print(f"   Drug volume = {total_daily_drug_dose_mg} mg / {admin_concentration} mg/ml = {drug_volume_ml_day} ml/day")
print("-" * 30)

# --- Step 3: Determine the volume needed for IV maintenance fluids ---
total_other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day
iv_maintenance_volume_ml_day = daily_fluid_maintenance_ml - total_other_fluids_ml_day

print(f"3. Determine the remaining volume for IV maintenance fluids:")
print(f"   Total intake from other sources = {drug_volume_ml_day} ml (chemo) + {enteral_feeding_ml_day} ml (milk) = {total_other_fluids_ml_day} ml/day")
print(f"   IV maintenance volume = {daily_fluid_maintenance_ml} ml (total) - {total_other_fluids_ml_day} ml (other) = {iv_maintenance_volume_ml_day} ml/day")
print("-" * 30)

# --- Step 4 & 5: Calculate the hourly infusion rate and round it ---
iv_rate_ml_hr = iv_maintenance_volume_ml_day / 24
rounded_iv_rate_ml_hr = int(round(iv_rate_ml_hr))

print(f"4. Calculate the hourly infusion rate for maintenance fluid:")
print(f"   IV rate = {iv_maintenance_volume_ml_day} ml / 24 hours = {iv_rate_ml_hr:.2f} ml/hr")
print(f"5. Round to the nearest whole number:")
print(f"   Rounded IV rate = {rounded_iv_rate_ml_hr} ml/hr")

# --- Final Answer ---
print("\nFinal Answer:")
final_equation = f"Round((({calc_1}) - ({bsa_m2} * {drug_dose_per_m2_per_day} / {admin_concentration}) - {enteral_feeding_ml_day}) / 24) = {rounded_iv_rate_ml_hr}"
print(f"The rate of the maintenance fluid (ml/hr) is calculated as: {final_equation}")
print(f"<<<{rounded_iv_rate_ml_hr}>>>")