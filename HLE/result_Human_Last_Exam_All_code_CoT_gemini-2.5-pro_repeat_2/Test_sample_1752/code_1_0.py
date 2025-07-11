import math

# Patient and Drug Information
weight_kg = 22
bsa_m2 = 0.8
recommended_dose_mg_m2_day = 25
admin_concentration_mg_ml = 1.0
enteral_feeding_ml_day = 500

# Step 1: Calculate the total daily dose of the drug
daily_dose_mg = recommended_dose_mg_m2_day * bsa_m2
print(f"Step 1: Daily drug dose = {recommended_dose_mg_m2_day} mg/m² * {bsa_m2} m² = {daily_dose_mg} mg/day")

# Step 2: Calculate the total daily volume for drug administration
drug_admin_volume_ml_day = daily_dose_mg / admin_concentration_mg_ml
print(f"Step 2: Daily drug administration volume = {daily_dose_mg} mg / {admin_concentration_mg_ml} mg/ml = {drug_admin_volume_ml_day} ml/day")

# Step 3: Calculate the child's total daily maintenance fluid requirement (Holliday-Segar method)
if weight_kg <= 10:
    maintenance_fluid_ml_day = weight_kg * 100
    print(f"Step 3: Total daily maintenance fluid = {weight_kg} kg * 100 ml/kg = {maintenance_fluid_ml_day} ml/day")
elif weight_kg <= 20:
    maintenance_fluid_ml_day = (10 * 100) + ((weight_kg - 10) * 50)
    print(f"Step 3: Total daily maintenance fluid = (10 kg * 100 ml/kg) + ({weight_kg - 10} kg * 50 ml/kg) = 1000 + {(weight_kg - 10) * 50} = {maintenance_fluid_ml_day} ml/day")
else:
    maintenance_fluid_ml_day = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    print(f"Step 3: Total daily maintenance fluid = (10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({weight_kg - 20} kg * 20 ml/kg) = 1000 + 500 + {(weight_kg - 20) * 20} = {maintenance_fluid_ml_day} ml/day")

# Step 4: Calculate the total fluid volume from other sources
total_other_fluids_ml_day = drug_admin_volume_ml_day + enteral_feeding_ml_day
print(f"Step 4: Total fluid from other sources = {drug_admin_volume_ml_day} ml (drug) + {enteral_feeding_ml_day} ml (milk) = {total_other_fluids_ml_day} ml/day")

# Step 5: Calculate the remaining fluid volume to be given as maintenance fluid
remaining_fluid_ml_day = maintenance_fluid_ml_day - total_other_fluids_ml_day
print(f"Step 5: Remaining maintenance fluid = {maintenance_fluid_ml_day} ml - {total_other_fluids_ml_day} ml = {remaining_fluid_ml_day} ml/day")

# Step 6: Calculate the hourly rate for the maintenance fluid
maintenance_rate_ml_hr = remaining_fluid_ml_day / 24
rounded_rate = round(maintenance_rate_ml_hr)
print(f"Step 6: Maintenance fluid rate = {remaining_fluid_ml_day} ml / 24 hours = {maintenance_rate_ml_hr:.1f} ml/hr")
print(f"\nFinal Answer: The rate of the maintenance fluid should be rounded to the nearest whole number.")
print(f"The final calculated rate is {rounded_rate} ml/hr.")

<<<43>>>