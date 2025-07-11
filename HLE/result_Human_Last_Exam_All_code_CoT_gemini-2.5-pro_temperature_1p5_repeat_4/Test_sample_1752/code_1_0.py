import math

# Patient and clinical information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m²/day
admin_concentration_mg_per_ml = 1
enteral_feeding_ml_per_day = 500

print("This script calculates the maintenance fluid rate for the patient.\n")

# Step 1: Calculate the total daily dose of the chemotherapy drug
print("Step 1: Calculate the total daily drug dose in mg.")
total_daily_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"The calculation is: {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_dose_mg} mg/day\n")

# Step 2: Calculate the volume of the drug to be administered each day
print("Step 2: Calculate the daily volume of the drug fluid.")
drug_volume_ml_per_day = total_daily_dose_mg / admin_concentration_mg_per_ml
print(f"The calculation is: {total_daily_dose_mg} mg / {admin_concentration_mg_per_ml} mg/ml = {drug_volume_ml_per_day} ml/day\n")

# Step 3: Calculate the total daily maintenance fluid requirement (Holliday-Segar)
print("Step 3: Calculate the total daily maintenance fluid requirement using the Holliday-Segar method.")
fluid_for_first_10kg = 10 * 100
fluid_for_next_10kg = 10 * 50
remaining_weight = weight_kg - 20
fluid_for_remaining_kg = remaining_weight * 20
maintenance_fluid_ml_per_day = fluid_for_first_10kg + fluid_for_next_10kg + fluid_for_remaining_kg
print(f"For the first 10 kg: 10 kg * 100 ml/kg = {fluid_for_first_10kg} ml")
print(f"For the next 10 kg (up to 20 kg): 10 kg * 50 ml/kg = {fluid_for_next_10kg} ml")
print(f"For the remaining {remaining_weight} kg: {remaining_weight} kg * 20 ml/kg = {fluid_for_remaining_kg} ml")
print(f"Total maintenance fluid: {fluid_for_first_10kg} ml + {fluid_for_next_10kg} ml + {fluid_for_remaining_kg} ml = {maintenance_fluid_ml_per_day} ml/day\n")

# Step 4: Calculate the total volume of fluid from other sources
print("Step 4: Calculate the total fluid volume from other sources (drug + milk).")
total_other_fluids_ml_per_day = drug_volume_ml_per_day + enteral_feeding_ml_per_day
print(f"The calculation is: {drug_volume_ml_per_day} ml (drug) + {enteral_feeding_ml_per_day} ml (milk) = {total_other_fluids_ml_per_day} ml/day\n")

# Step 5: Calculate the remaining volume needed as IV maintenance fluid
print("Step 5: Calculate the remaining fluid volume to be given as IV maintenance fluid.")
remaining_iv_fluid_ml_per_day = maintenance_fluid_ml_per_day - total_other_fluids_ml_per_day
print(f"The calculation is: {maintenance_fluid_ml_per_day} ml (total) - {total_other_fluids_ml_per_day} ml (other sources) = {remaining_iv_fluid_ml_per_day} ml/day\n")

# Step 6: Calculate the rate of the maintenance fluid in ml/hr
print("Step 6: Calculate the hourly rate for the IV maintenance fluid.")
maintenance_fluid_rate_ml_hr = remaining_iv_fluid_ml_per_day / 24
rounded_rate = round(maintenance_fluid_rate_ml_hr)
print(f"The calculation is: {remaining_iv_fluid_ml_per_day} ml / 24 hours = {maintenance_fluid_rate_ml_hr:.2f} ml/hr")
print(f"The final rate rounded to the nearest whole number is {rounded_rate} ml/hr.\n")

# Final summary equation
print("Final Summary Equation:")
# ( ( (10 * 100) + (10 * 50) + ((22 - 20) * 20) ) - ( ((25 * 0.8) / 1) + 500 ) ) / 24
print(f"( ( (10 * {100}) + (10 * {50}) + (({weight_kg} - 20) * {20}) ) - ( (({drug_dose_per_m2_per_day} * {bsa_m2}) / {admin_concentration_mg_per_ml}) + {enteral_feeding_ml_per_day} ) ) / 24 = {rounded_rate}")
<<<43>>>