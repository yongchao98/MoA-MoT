import math

# Patient and Drug Information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2 = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate the total daily dose of the chemotherapy drug
print("Step 1: Calculate the total daily dose (mg/day)")
daily_dose_mg = drug_dose_per_m2 * bsa_m2
print(f"Daily Dose = {drug_dose_per_m2} mg/m² * {bsa_m2} m² = {daily_dose_mg} mg/day")
print("-" * 30)

# Step 2: Calculate the volume of the drug solution to be administered daily
print("Step 2: Calculate the daily drug volume (ml/day)")
drug_volume_ml_day = daily_dose_mg / admin_concentration
print(f"Drug Volume = {daily_dose_mg} mg / {admin_concentration} mg/ml = {drug_volume_ml_day} ml/day")
print("-" * 30)

# Step 3: Calculate the total daily maintenance fluid requirement using Holliday-Segar method
print("Step 3: Calculate total daily maintenance fluid (ml/day) using Holliday-Segar")
if weight_kg <= 10:
    maintenance_fluid_day = weight_kg * 100
    print(f"Total Daily Fluid = {weight_kg} kg * 100 ml/kg = {maintenance_fluid_day} ml/day")
elif weight_kg <= 20:
    maintenance_fluid_day = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Total Daily Fluid = (10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg) = 1000 + {(weight_kg - 10) * 50} = {maintenance_fluid_day} ml/day")
else:
    part1 = 10 * 100
    part2 = 10 * 50
    part3 = (weight_kg - 20) * 20
    maintenance_fluid_day = part1 + part2 + part3
    print(f"Total Daily Fluid = (10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({weight_kg - 20} kg * 20 ml/kg)")
    print(f"= {part1} + {part2} + {part3} = {maintenance_fluid_day} ml/day")
print("-" * 30)

# Step 4: Calculate total fluid intake from other sources
print("Step 4: Calculate total fluid from other sources (drug + milk)")
other_fluids_day = drug_volume_ml_day + enteral_feeding_ml_day
print(f"Other Fluids = {drug_volume_ml_day} ml (drug) + {enteral_feeding_ml_day} ml (milk) = {other_fluids_day} ml/day")
print("-" * 30)

# Step 5: Calculate the remaining volume needed for IV maintenance fluids
print("Step 5: Calculate remaining IV fluid volume needed per day")
remaining_iv_volume_day = maintenance_fluid_day - other_fluids_day
print(f"Remaining IV Volume = {maintenance_fluid_day} ml - {other_fluids_day} ml = {remaining_iv_volume_day} ml/day")
print("-" * 30)

# Step 6: Calculate the IV infusion rate in ml/hr
print("Step 6: Calculate the IV infusion rate (ml/hr) and round to the nearest whole number")
iv_rate_ml_hr = remaining_iv_volume_day / 24
rounded_iv_rate = round(iv_rate_ml_hr)
print(f"IV Rate = {remaining_iv_volume_day} ml / 24 hours = {iv_rate_ml_hr:.1f} ml/hr")
print(f"Rounded IV Rate = {rounded_iv_rate} ml/hr")

# Final Answer
print(f"\nThe final rate of the maintenance fluid is {rounded_iv_rate} ml/hr.")
print(f"<<<{rounded_iv_rate}>>>")