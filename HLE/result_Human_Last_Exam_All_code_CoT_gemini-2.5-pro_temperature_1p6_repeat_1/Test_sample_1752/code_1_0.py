import math

# Patient Information
weight_kg = 22
bsa_m2 = 0.8
enteral_feeding_ml_day = 500

# Drug Information
drug_dose_per_m2_day = 25 # mg/m^2/day
admin_concentration_mg_ml = 1 # mg/ml

# Step 1: Calculate total daily fluid requirement using Holliday-Segar method
if weight_kg <= 10:
    daily_fluid_req_ml = weight_kg * 100
    calculation_str = f"{weight_kg} kg * 100 ml/kg"
elif weight_kg <= 20:
    daily_fluid_req_ml = (10 * 100) + ((weight_kg - 10) * 50)
    calculation_str = f"(10 kg * 100 ml/kg) + ({(weight_kg - 10)} kg * 50 ml/kg)"
else:
    daily_fluid_req_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    first_10_kg_fluid = 10 * 100
    next_10_kg_fluid = 10 * 50
    remaining_weight = weight_kg - 20
    remaining_fluid = remaining_weight * 20
    calculation_str = f"({first_10_kg_fluid}) + ({next_10_kg_fluid}) + ({remaining_fluid})"

print(f"Step 1: Calculate total daily fluid requirement (Holliday-Segar):")
print(f"Calculation: {calculation_str} = {daily_fluid_req_ml} ml/day\n")

# Step 2: Calculate the daily volume from the chemotherapy drug
total_drug_dose_mg_day = drug_dose_per_m2_day * bsa_m2
drug_volume_ml_day = total_drug_dose_mg_day / admin_concentration_mg_ml

print(f"Step 2: Calculate daily drug volume:")
print(f"Total daily drug dose = {drug_dose_per_m2_day} mg/m² * {bsa_m2} m² = {total_drug_dose_mg_day} mg/day")
print(f"Drug volume = {total_drug_dose_mg_day} mg / {admin_concentration_mg_ml} mg/ml = {drug_volume_ml_day} ml/day\n")

# Step 3: Sum up non-maintenance fluid volumes
other_fluids_ml_day = drug_volume_ml_day + enteral_feeding_ml_day
print(f"Step 3: Sum all non-maintenance fluids:")
print(f"Total from other sources = {drug_volume_ml_day} ml (drug) + {enteral_feeding_ml_day} ml (milk) = {other_fluids_ml_day} ml/day\n")

# Step 4: Determine the remaining fluid volume
remaining_ivf_ml_day = daily_fluid_req_ml - other_fluids_ml_day
print(f"Step 4: Calculate remaining IV fluid volume needed:")
print(f"Remaining IV fluid = {daily_fluid_req_ml} ml - {other_fluids_ml_day} ml = {remaining_ivf_ml_day} ml/day\n")

# Step 5: Calculate the hourly rate
hourly_rate = remaining_ivf_ml_day / 24
rounded_hourly_rate = round(hourly_rate)
print(f"Step 5: Calculate the hourly rate:")
print(f"Hourly rate = {remaining_ivf_ml_day} ml / 24 hours = {hourly_rate:.1f} ml/hr")
print(f"Rounded to the nearest whole number, the final rate is: {rounded_hourly_rate} ml/hr")

<<<43>>>