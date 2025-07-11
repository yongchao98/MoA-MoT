import math

# Step 1: Define the given variables
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration = 1  # mg/ml
milk_volume_per_day = 500  # ml/day

# Step 2: Calculate the total daily dose of the drug
total_daily_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"The total daily drug dose is {drug_dose_per_m2_per_day} mg/m^2 * {bsa_m2} m^2 = {total_daily_dose_mg} mg.")

# Step 3: Calculate the daily volume of the drug infusion
drug_volume_per_day = total_daily_dose_mg / admin_concentration
print(f"The daily volume of the drug infusion is {total_daily_dose_mg} mg / {admin_concentration} mg/ml = {drug_volume_per_day} ml.")

# Step 4: Calculate the total daily maintenance fluid requirement using Holliday-Segar method
if weight_kg <= 10:
    maintenance_fluid_per_day = weight_kg * 100
    print(f"Daily maintenance fluid: {weight_kg} kg * 100 ml/kg = {maintenance_fluid_per_day} ml.")
elif weight_kg <= 20:
    maintenance_fluid_per_day = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Daily maintenance fluid: (10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg) = {maintenance_fluid_per_day} ml.")
else:
    part1 = 10 * 100
    part2 = 10 * 50
    part3 = (weight_kg - 20) * 20
    maintenance_fluid_per_day = part1 + part2 + part3
    print(f"Daily maintenance fluid: ({part1}) + ({part2}) + ({part3}) = {maintenance_fluid_per_day} ml.")

# Step 5: Calculate the total fluid from other sources
total_other_fluids = drug_volume_per_day + milk_volume_per_day
print(f"Total fluid from other sources (drug + milk) is {drug_volume_per_day} ml + {milk_volume_per_day} ml = {total_other_fluids} ml.")

# Step 6: Determine the remaining fluid volume for IV maintenance
remaining_iv_fluid = maintenance_fluid_per_day - total_other_fluids
print(f"Remaining IV fluid needed is {maintenance_fluid_per_day} ml - {total_other_fluids} ml = {remaining_iv_fluid} ml.")

# Step 7: Calculate the hourly rate and round to the nearest whole number
hourly_rate = remaining_iv_fluid / 24
rounded_hourly_rate = round(hourly_rate)
print(f"The required IV maintenance fluid rate is {remaining_iv_fluid} ml / 24 hours = {hourly_rate:.1f} ml/hr.")
print(f"Rounding to the nearest whole number, the final rate is {rounded_hourly_rate} ml/hr.")

print("Final Answer:")
print(f"({maintenance_fluid_per_day} - ({drug_volume_per_day} + {milk_volume_per_day})) / 24 = {rounded_hourly_rate}")
<<<43>>>