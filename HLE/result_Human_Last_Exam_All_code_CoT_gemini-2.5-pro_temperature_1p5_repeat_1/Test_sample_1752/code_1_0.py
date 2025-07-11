import math

# Given information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2 = 25  # mg/m2/day
admin_concentration = 1 # mg/ml
milk_volume_ml_day = 500

# --- Step 1: Calculate total daily maintenance fluid using Holliday-Segar method ---
if weight_kg <= 10:
    daily_maintenance_fluid = weight_kg * 100
    calc_str_1 = f"{weight_kg} kg * 100 ml/kg"
elif weight_kg <= 20:
    daily_maintenance_fluid = (10 * 100) + ((weight_kg - 10) * 50)
    calc_str_1 = f"(10 kg * 100 ml/kg) + (({weight_kg} - 10) kg * 50 ml/kg)"
else:
    daily_maintenance_fluid = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    calc_str_1 = f"(10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + (({weight_kg} - 20) kg * 20 ml/kg)"
first_10_kg_fluid = 10 * 100
next_10_kg_fluid = 10 * 50
remaining_weight = weight_kg - 20
remaining_fluid = remaining_weight * 20

print("Step 1: Calculate the patient's total daily maintenance fluid requirement.")
print(f"Calculation based on Holliday-Segar method for a {weight_kg} kg patient:")
print(f"For the first 10 kg: 10 kg * 100 ml/kg = {first_10_kg_fluid} ml")
print(f"For the next 10 kg: 10 kg * 50 ml/kg = {next_10_kg_fluid} ml")
print(f"For the remaining {remaining_weight} kg: {remaining_weight} kg * 20 ml/kg = {remaining_fluid} ml")
print(f"Total Daily Maintenance Fluid = {first_10_kg_fluid} + {next_10_kg_fluid} + {remaining_fluid} = {daily_maintenance_fluid} ml/day\n")


# --- Step 2: Calculate the fluid volume from the chemotherapy drug ---
daily_drug_dose_mg = drug_dose_per_m2 * bsa_m2
drug_infusion_volume_ml_day = daily_drug_dose_mg / admin_concentration

print("Step 2: Calculate the fluid volume from the daily chemotherapy infusion.")
print(f"Daily Drug Dose = {drug_dose_per_m2} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg")
print(f"Required infusion volume for administration at {admin_concentration} mg/ml:")
print(f"Volume = {daily_drug_dose_mg} mg / {admin_concentration} mg/ml = {drug_infusion_volume_ml_day} ml/day\n")

# --- Step 3: Calculate total fluid from other sources ---
total_other_fluids_ml_day = drug_infusion_volume_ml_day + milk_volume_ml_day

print("Step 3: Calculate total daily fluid from other sources (drug + milk).")
print(f"Total Other Fluids = {drug_infusion_volume_ml_day} ml (chemo) + {milk_volume_ml_day} ml (milk) = {total_other_fluids_ml_day} ml/day\n")

# --- Step 4: Calculate the remaining volume for IV maintenance fluid ---
remaining_iv_fluid_ml_day = daily_maintenance_fluid - total_other_fluids_ml_day

print("Step 4: Calculate the remaining volume needed as maintenance fluid.")
print(f"Remaining Volume = {daily_maintenance_fluid} ml (total) - {total_other_fluids_ml_day} ml (other sources) = {remaining_iv_fluid_ml_day} ml/day\n")

# --- Step 5 & 6: Calculate the hourly rate and round it ---
iv_rate_ml_hr = remaining_iv_fluid_ml_day / 24
rounded_iv_rate_ml_hr = round(iv_rate_ml_hr)

print("Step 5 & 6: Calculate the hourly rate for the maintenance fluid and round it.")
print(f"Maintenance Fluid Rate = {remaining_iv_fluid_ml_day} ml / 24 hours = {iv_rate_ml_hr:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the final rate is {rounded_iv_rate_ml_hr} ml/hr.")

# --- Final Answer ---
print("\nFinal Answer:")
print(f"({first_10_kg_fluid} + {next_10_kg_fluid} + {remaining_fluid} - {drug_infusion_volume_ml_day} - {milk_volume_ml_day}) / 24 = {rounded_iv_rate_ml_hr}")
print(f"({daily_maintenance_fluid} - {total_other_fluids_ml_day}) / 24 = {rounded_iv_rate_ml_hr}")
print(f"{remaining_iv_fluid_ml_day} / 24 = {rounded_iv_rate_ml_hr}")
print(f"The rate of his maintenance fluid should be {rounded_iv_rate_ml_hr} ml/hr.")

<<<43>>>