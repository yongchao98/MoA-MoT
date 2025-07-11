import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_mg_per_m2_day = 25
admin_concentration_mg_per_ml = 1
milk_volume_ml_day = 500
hours_in_day = 24

# Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
print("Step 1: Calculate total daily maintenance fluid requirement (Holliday-Segar method).")
fluid_for_first_10kg = 100 * 10
fluid_for_next_10kg = 50 * 10
fluid_for_remaining_kg = 20 * (weight_kg - 20)
total_maintenance_fluid_ml_day = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg

print(f"Based on a weight of {weight_kg} kg:")
print(f"(10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({weight_kg - 20} kg * 20 ml/kg) = {fluid_for_first_10kg} + {fluid_for_next_10kg} + {fluid_for_remaining_kg} = {total_maintenance_fluid_ml_day} ml/day\n")

# Step 2: Calculate the volume of all non-maintenance fluids
print("Step 2: Calculate the volume of all non-maintenance fluids.")
# Calculate daily drug dose and volume
daily_drug_dose_mg = drug_dose_mg_per_m2_day * bsa_m2
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration_mg_per_ml
print(f"Daily drug volume = ({drug_dose_mg_per_m2_day} mg/m² * {bsa_m2} m²) / {admin_concentration_mg_per_ml} mg/ml = {daily_drug_volume_ml} ml")

# Calculate total other fluids
total_other_fluids_ml_day = daily_drug_volume_ml + milk_volume_ml_day
print(f"Total other fluids = Drug volume + Milk volume = {daily_drug_volume_ml} ml + {milk_volume_ml_day} ml = {total_other_fluids_ml_day} ml/day\n")

# Step 3: Determine the remaining fluid volume for IV maintenance
print("Step 3: Determine the remaining fluid volume needed as IV maintenance.")
remaining_iv_fluid_ml_day = total_maintenance_fluid_ml_day - total_other_fluids_ml_day
print(f"Remaining IV fluid = Total daily fluid - Total other fluids")
print(f"Remaining IV fluid = {total_maintenance_fluid_ml_day} ml - {total_other_fluids_ml_day} ml = {remaining_iv_fluid_ml_day} ml/day\n")

# Step 4: Calculate the hourly rate
print("Step 4: Calculate the hourly infusion rate for the maintenance fluid.")
maintenance_rate_ml_hr = remaining_iv_fluid_ml_day / hours_in_day
print(f"Hourly rate = Remaining IV fluid / Hours in a day")
print(f"Hourly rate = {remaining_iv_fluid_ml_day} ml / {hours_in_day} hr = {maintenance_rate_ml_hr:.2f} ml/hr\n")

# Round to the nearest whole number
final_rate = round(maintenance_rate_ml_hr)
print(f"The final calculated rate is {maintenance_rate_ml_hr:.2f} ml/hr, which rounds to {final_rate} ml/hr.")

print(f"\nFinal Answer: The maintenance fluid should be run at a rate of {final_rate} ml/hr.")
print(f"<<<{final_rate}>>>")