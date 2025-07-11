import math

# Patient Information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
drug_admin_concentration = 1  # mg/ml
enteral_feeding_ml_day = 500  # ml/day

# Step 1: Calculate total daily fluid requirement using Holliday-Segar method
total_daily_fluid_requirement = 0
if weight_kg <= 10:
    total_daily_fluid_requirement = weight_kg * 100
elif weight_kg <= 20:
    total_daily_fluid_requirement = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_daily_fluid_requirement = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

print(f"1. Calculating Total Daily Fluid Requirement (Holliday-Segar):")
print(f"   - For the first 10 kg: 10 kg * 100 ml/kg = {10 * 100} ml")
print(f"   - For the next 10 kg: 10 kg * 50 ml/kg = {10 * 50} ml")
print(f"   - For the remaining {weight_kg - 20} kg: {weight_kg - 20} kg * 20 ml/kg = {(weight_kg - 20) * 20} ml")
print(f"   - Total Daily Requirement: {10*100} + {10*50} + {(weight_kg - 20) * 20} = {total_daily_fluid_requirement} ml/day\n")


# Step 2: Calculate the volume of fluid from the chemotherapy drug per day
total_drug_dose_mg_day = drug_dose_per_m2_per_day * bsa_m2
drug_volume_ml_day = total_drug_dose_mg_day / drug_admin_concentration

print(f"2. Calculating Daily Fluid Volume from Chemotherapy:")
print(f"   - Total daily dose (mg): {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_drug_dose_mg_day} mg/day")
print(f"   - Drug volume (ml): {total_drug_dose_mg_day} mg / {drug_admin_concentration} mg/ml = {drug_volume_ml_day} ml/day\n")


# Step 3: Sum all known fluid sources
total_known_intake_ml_day = enteral_feeding_ml_day + drug_volume_ml_day
print(f"3. Calculating Total Known Fluid Intake:")
print(f"   - Enteral Feeding + Drug Volume: {enteral_feeding_ml_day} ml/day + {drug_volume_ml_day} ml/day = {total_known_intake_ml_day} ml/day\n")


# Step 4: Determine the remaining fluid volume
remaining_iv_fluid_ml_day = total_daily_fluid_requirement - total_known_intake_ml_day
print(f"4. Calculating Remaining IV Maintenance Fluid:")
print(f"   - Total Requirement - Known Intake: {total_daily_fluid_requirement} ml/day - {total_known_intake_ml_day} ml/day = {remaining_iv_fluid_ml_day} ml/day\n")


# Step 5: Calculate the hourly rate
hourly_rate = remaining_iv_fluid_ml_day / 24
print(f"5. Calculating the Final IV Rate (ml/hr):")
print(f"   - Equation: ({total_daily_fluid_requirement} - {total_known_intake_ml_day}) / 24 hours")
print(f"   - Calculation: {remaining_iv_fluid_ml_day} ml / 24 hr = {hourly_rate:.1f} ml/hr\n")


# Step 6: Round the final answer to the nearest whole number
rounded_hourly_rate = round(hourly_rate)
print(f"The final rate for the maintenance fluid, rounded to the nearest whole number, is {rounded_hourly_rate} ml/hr.")

print(f"\n<<< {rounded_hourly_rate} >>>")