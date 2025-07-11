import math

# --- Given Information ---
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2 = 25  # mg/m²/day
drug_admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500 # ml/day

# Step 1: Calculate total daily fluid requirement using Holliday-Segar method
if weight_kg <= 10:
    total_maintenance_fluid_ml_day = weight_kg * 100
    print(f"Calculating maintenance fluid for first 10 kg: {weight_kg} kg * 100 ml/kg = {total_maintenance_fluid_ml_day} ml/day")
elif weight_kg <= 20:
    total_maintenance_fluid_ml_day = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Calculating maintenance fluid: (10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg) = {total_maintenance_fluid_ml_day} ml/day")
else: # weight > 20
    calc_first_10 = 10 * 100
    calc_next_10 = 10 * 50
    calc_remainder = (weight_kg - 20) * 20
    total_maintenance_fluid_ml_day = calc_first_10 + calc_next_10 + calc_remainder
    print("Calculating maintenance fluid for a patient over 20 kg:")
    print(f"  First 10 kg: 10 kg * 100 ml/kg = {calc_first_10} ml")
    print(f"  Next 10 kg: 10 kg * 50 ml/kg = {calc_next_10} ml")
    print(f"  Remaining {weight_kg - 20} kg: {weight_kg - 20} kg * 20 ml/kg = {calc_remainder} ml")
    print(f"Total Daily Fluid Requirement: {calc_first_10} + {calc_next_10} + {calc_remainder} = {total_maintenance_fluid_ml_day} ml/day\n")

# Step 2: Calculate the fluid volume from the chemotherapy drug
total_daily_drug_dose_mg = drug_dose_per_m2 * bsa_m2
drug_volume_ml_day = total_daily_drug_dose_mg / drug_admin_concentration
print("Calculating fluid volume from medication:")
print(f"  Daily Drug Dose (mg): {drug_dose_per_m2} mg/m² * {bsa_m2} m² = {total_daily_drug_dose_mg} mg/day")
print(f"  Daily Drug Volume (ml): {total_daily_drug_dose_mg} mg / {drug_admin_concentration} mg/ml = {drug_volume_ml_day} ml/day\n")

# Step 3: Calculate total fluid intake from other sources
other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day
print("Calculating total fluid from other sources:")
print(f"  Total = Drug Volume + Enteral Feed Volume")
print(f"  Total = {drug_volume_ml_day} ml/day + {enteral_feeding_ml_day} ml/day = {other_fluids_ml_day} ml/day\n")

# Step 4: Determine the remaining fluid volume for IV maintenance
iv_maintenance_ml_day = total_maintenance_fluid_ml_day - other_fluids_ml_day
print("Calculating remaining IV maintenance fluid:")
print(f"  Remaining Volume = Total Requirement - Other Fluids")
print(f"  Remaining Volume = {total_maintenance_fluid_ml_day} ml/day - {other_fluids_ml_day} ml/day = {iv_maintenance_ml_day} ml/day\n")

# Step 5 & 6: Calculate and round the hourly rate
iv_rate_ml_hr = iv_maintenance_ml_day / 24
rounded_iv_rate_ml_hr = round(iv_rate_ml_hr)
print("Calculating the final infusion rate:")
print(f"  Rate (ml/hr) = Remaining Volume / 24 hours")
print(f"  Rate (ml/hr) = {iv_maintenance_ml_day} ml / 24 hr = {iv_rate_ml_hr:.2f} ml/hr")
print(f"\nRounding to the nearest whole number, the final rate for the maintenance fluid is {rounded_iv_rate_ml_hr} ml/hr.")
print("<<<43>>>")