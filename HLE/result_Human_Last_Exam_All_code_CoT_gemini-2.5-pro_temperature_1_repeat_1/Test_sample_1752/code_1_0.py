import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
chemo_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
if weight_kg <= 10:
    daily_maintenance_fluid = weight_kg * 100
    equation_part1 = f"{weight_kg} kg * 100 ml/kg = {daily_maintenance_fluid} ml"
    print(f"Calculating maintenance fluid for first 10kg: {equation_part1}")
    total_daily_fluid = daily_maintenance_fluid
elif weight_kg <= 20:
    fluid_first_10kg = 10 * 100
    fluid_next_kg = (weight_kg - 10) * 50
    total_daily_fluid = fluid_first_10kg + fluid_next_kg
    print(f"Calculating maintenance fluid for first 10kg: 10 kg * 100 ml/kg = {fluid_first_10kg} ml")
    print(f"Calculating maintenance fluid for next {weight_kg - 10} kg: {weight_kg - 10} kg * 50 ml/kg = {fluid_next_kg} ml")
    print(f"Total Daily Maintenance Fluid = {fluid_first_10kg} ml + {fluid_next_kg} ml = {total_daily_fluid} ml/day")
else:
    fluid_first_10kg = 10 * 100
    fluid_next_10kg = 10 * 50
    fluid_remaining_kg = (weight_kg - 20) * 20
    total_daily_fluid = fluid_first_10kg + fluid_next_10kg + fluid_remaining_kg
    print(f"Calculating maintenance fluid for first 10kg: 10 kg * 100 ml/kg = {fluid_first_10kg} ml")
    print(f"Calculating maintenance fluid for next 10kg: 10 kg * 50 ml/kg = {fluid_next_10kg} ml")
    print(f"Calculating maintenance fluid for remaining {weight_kg - 20} kg: {weight_kg - 20} kg * 20 ml/kg = {fluid_remaining_kg} ml")
    print(f"Total Daily Maintenance Fluid = {fluid_first_10kg} ml + {fluid_next_10kg} ml + {fluid_remaining_kg} ml = {total_daily_fluid} ml/day")

# Step 2: Calculate daily chemotherapy dose and administration volume
daily_chemo_dose_mg = chemo_dose_per_m2_per_day * bsa_m2
print(f"\nCalculating daily chemotherapy dose: {chemo_dose_per_m2_per_day} mg/m^2 * {bsa_m2} m^2 = {daily_chemo_dose_mg} mg/day")

daily_chemo_volume_ml = daily_chemo_dose_mg / admin_concentration
print(f"Calculating daily chemotherapy fluid volume: {daily_chemo_dose_mg} mg / {admin_concentration} mg/ml = {daily_chemo_volume_ml} ml/day")

# Step 3: Calculate total fluid from other sources
total_other_sources_ml = daily_chemo_volume_ml + enteral_feeding_ml_day
print(f"\nCalculating total fluid from other sources: {daily_chemo_volume_ml} ml (chemo) + {enteral_feeding_ml_day} ml (milk) = {total_other_sources_ml} ml/day")

# Step 4: Calculate the remaining maintenance fluid volume
remaining_fluid_ml_day = total_daily_fluid - total_other_sources_ml
print(f"\nCalculating remaining fluid for IV maintenance: {total_daily_fluid} ml (total) - {total_other_sources_ml} ml (other sources) = {remaining_fluid_ml_day} ml/day")

# Step 5: Calculate the hourly rate
hourly_rate = remaining_fluid_ml_day / 24
print(f"\nCalculating the hourly rate: {remaining_fluid_ml_day} ml / 24 hours = {hourly_rate:.1f} ml/hr")

# Step 6: Round to the nearest whole number
final_rate = round(hourly_rate)
print(f"\nRounding to the nearest whole number: round({hourly_rate:.1f}) = {final_rate}")
print(f"\nThe final rate for the maintenance fluid is {final_rate} ml/hr.")

print(f"<<<{final_rate}>>>")