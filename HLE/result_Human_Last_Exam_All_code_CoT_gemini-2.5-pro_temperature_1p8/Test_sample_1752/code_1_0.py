import math

# Patient and Drug Information
weight_kg = 22
bsa_m2 = 0.8
enteral_feeding_ml_day = 500
drug_dose_per_m2_day = 25
admin_concentration_mg_ml = 1

# Step 1: Calculate the patient's total daily maintenance fluid requirement
# Using the Holliday-Segar method for a 22 kg child:
# 100 ml/kg for the first 10 kg
# 50 ml/kg for the next 10 kg
# 20 ml/kg for the remaining weight
first_10_kg_fluid = 10 * 100
next_10_kg_fluid = 10 * 50
remaining_weight_fluid = (weight_kg - 20) * 20
daily_fluid_req_ml = first_10_kg_fluid + next_10_kg_fluid + remaining_weight_fluid

print("Step 1: Calculate total daily fluid requirement using the Holliday-Segar method.")
print(f"Calculation: (10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({weight_kg - 20} kg * 20 ml/kg) = {first_10_kg_fluid} + {next_10_kg_fluid} + {remaining_weight_fluid} = {daily_fluid_req_ml} ml/day\n")

# Step 2: Calculate the volume of the chemotherapy infusion per day
daily_drug_dose_mg = drug_dose_per_m2_day * bsa_m2
drug_volume_ml_day = daily_drug_dose_mg / admin_concentration_mg_ml

print("Step 2: Calculate the volume of the daily chemotherapy infusion.")
print(f"Daily Drug Dose (mg) = {drug_dose_per_m2_day} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg")
print(f"Daily Infusion Volume (ml) = {daily_drug_dose_mg} mg / {admin_concentration_mg_ml} mg/ml = {drug_volume_ml_day} ml\n")

# Step 3: Calculate total fluid already being received
total_other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day

print("Step 3: Calculate total fluid received from other sources (chemo + milk).")
print(f"Total Other Fluids = {drug_volume_ml_day} ml (chemo) + {enteral_feeding_ml_day} ml (milk) = {total_other_fluids_ml_day} ml/day\n")

# Step 4: Calculate the remaining fluid volume for IV maintenance
remaining_fluid_ml_day = daily_fluid_req_ml - total_other_fluids_ml_day

print("Step 4: Calculate the remaining fluid volume to be given via IV.")
print(f"Remaining IV Volume = {daily_fluid_req_ml} ml/day - {total_other_fluids_ml_day} ml/day = {remaining_fluid_ml_day} ml/day\n")

# Step 5: Calculate the hourly rate and round it
maintenance_rate_ml_hr = remaining_fluid_ml_day / 24
rounded_rate = round(maintenance_rate_ml_hr)

print("Step 5: Calculate the hourly rate for the maintenance fluid and round.")
print(f"IV Maintenance Rate (ml/hr) = {remaining_fluid_ml_day} ml / 24 hr = {maintenance_rate_ml_hr:.1f} ml/hr")
print(f"\nThe final rate for the maintenance fluid, rounded to the nearest whole number, is {rounded_rate} ml/hr.")

# Final Answer
print(f'<<<{rounded_rate}>>>')