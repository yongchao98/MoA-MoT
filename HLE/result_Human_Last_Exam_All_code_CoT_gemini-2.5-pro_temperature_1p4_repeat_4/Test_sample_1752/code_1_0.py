import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1.0  # mg/ml
enteral_feeding_ml_day = 500 # ml/day

# Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    total_daily_fluid_ml = weight_kg * 100
    print(f"Calculating maintenance fluid for the first 10 kg: {weight_kg} kg * 100 ml/kg = {total_daily_fluid_ml} ml")
elif weight_kg <= 20:
    total_daily_fluid_ml = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Calculating maintenance fluid for the first 10 kg: 10 kg * 100 ml/kg = 1000 ml")
    print(f"Calculating maintenance fluid for the next {weight_kg - 10} kg: {weight_kg - 10} kg * 50 ml/kg = {(weight_kg - 10) * 50} ml")
else:
    fluid_for_first_10kg = 10 * 100
    fluid_for_next_10kg = 10 * 50
    fluid_for_remaining_kg = (weight_kg - 20) * 20
    total_daily_fluid_ml = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg
    print(f"Calculating maintenance fluid for the first 10 kg: 10 kg * 100 ml/kg = {fluid_for_first_10kg} ml")
    print(f"Calculating maintenance fluid for the next 10 kg: 10 kg * 50 ml/kg = {fluid_for_next_10kg} ml")
    print(f"Calculating maintenance fluid for the remaining {weight_kg - 20} kg: {weight_kg - 20} kg * 20 ml/kg = {fluid_for_remaining_kg} ml")

print(f"Total Daily Maintenance Fluid Requirement = {total_daily_fluid_ml} ml/day\n")

# Step 2: Calculate the total daily dose of the chemotherapy drug in mg
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"Calculating daily drug dose: {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg/day\n")

# Step 3: Calculate the volume of the daily chemotherapy dose
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration
print(f"Calculating daily drug volume: {daily_drug_dose_mg} mg / {admin_concentration} mg/ml = {daily_drug_volume_ml} ml/day\n")

# Step 4: Sum all non-maintenance fluid volumes
total_other_fluids_ml = daily_drug_volume_ml + enteral_feeding_ml_day
print(f"Calculating total fluid from other sources: {daily_drug_volume_ml} ml (drug) + {enteral_feeding_ml_day} ml (milk) = {total_other_fluids_ml} ml/day\n")

# Step 5: Determine the remaining fluid volume for maintenance
remaining_maintenance_fluid_ml = total_daily_fluid_ml - total_other_fluids_ml
print(f"Calculating remaining fluid for maintenance IV: {total_daily_fluid_ml} ml (total) - {total_other_fluids_ml} ml (other sources) = {remaining_maintenance_fluid_ml} ml/day\n")

# Step 6: Calculate the hourly rate
hourly_rate = remaining_maintenance_fluid_ml / 24
rounded_hourly_rate = round(hourly_rate)
print(f"Calculating the hourly rate: {remaining_maintenance_fluid_ml} ml / 24 hours = {hourly_rate:.1f} ml/hr")
print(f"Rounding to the nearest whole number gives {rounded_hourly_rate} ml/hr.")

print(f"\nThe final equation is: round( ( (10 * 100 + 10 * 50 + ( {weight_kg} - 20 ) * 20) - ( ( {drug_dose_per_m2_per_day} * {bsa_m2} ) / {admin_concentration} + {enteral_feeding_ml_day} ) ) / 24 ) = {rounded_hourly_rate}")

<<<43>>>