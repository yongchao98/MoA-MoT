import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
chemo_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500  # ml/day

# Step 1: Calculate daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    daily_maintenance_fluid = weight_kg * 100
    print(f"Calculating maintenance fluid for first 10 kg: {weight_kg} kg * 100 ml/kg = {daily_maintenance_fluid} ml/day")
elif weight_kg <= 20:
    daily_maintenance_fluid = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Calculating maintenance fluid for 22 kg child:")
    print(f"(10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg) = 1000 + {(weight_kg - 10) * 50} = {daily_maintenance_fluid} ml/day")
else: # weight > 20
    fluid_for_first_10kg = 10 * 100
    fluid_for_next_10kg = 10 * 50
    fluid_for_remainder = (weight_kg - 20) * 20
    daily_maintenance_fluid = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remainder
    print(f"Calculating total daily maintenance fluid for a {weight_kg} kg child (Holliday-Segar):")
    print(f"For the first 10 kg: 10 kg * 100 ml/kg = {fluid_for_first_10kg} ml")
    print(f"For 11-20 kg: 10 kg * 50 ml/kg = {fluid_for_next_10kg} ml")
    print(f"For the remaining {weight_kg - 20} kg: {weight_kg - 20} kg * 20 ml/kg = {fluid_for_remainder} ml")
    print(f"Total maintenance fluid = {fluid_for_first_10kg} + {fluid_for_next_10kg} + {fluid_for_remainder} = {daily_maintenance_fluid} ml/day\n")

# Step 2: Calculate the volume of the daily chemotherapy dose
total_daily_chemo_dose_mg = chemo_dose_per_m2_per_day * bsa_m2
chemo_volume_ml_day = total_daily_chemo_dose_mg / admin_concentration
print(f"Calculating daily chemotherapy volume:")
print(f"Total daily dose = {chemo_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_chemo_dose_mg} mg/day")
print(f"Volume for administration = {total_daily_chemo_dose_mg} mg / {admin_concentration} mg/ml = {chemo_volume_ml_day} ml/day\n")

# Step 3: Calculate total fluid from other sources
total_other_fluids = chemo_volume_ml_day + enteral_feeding_ml_day
print(f"Calculating total fluid from other sources:")
print(f"Total = Chemo volume + Enteral feeding = {chemo_volume_ml_day} ml + {enteral_feeding_ml_day} ml = {total_other_fluids} ml/day\n")

# Step 4: Calculate the remaining maintenance fluid volume
remaining_fluid_volume = daily_maintenance_fluid - total_other_fluids
print(f"Calculating remaining fluid for IV maintenance:")
print(f"Remaining fluid = Total fluid requirement - Other fluids = {daily_maintenance_fluid} ml - {total_other_fluids} ml = {remaining_fluid_volume} ml/day\n")

# Step 5: Calculate the hourly rate and round it
hourly_rate = remaining_fluid_volume / 24
final_rate = round(hourly_rate)
print(f"Calculating the final hourly rate:")
print(f"Rate = {remaining_fluid_volume} ml / 24 hours = {hourly_rate:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the final rate is {final_rate} ml/hr.")
