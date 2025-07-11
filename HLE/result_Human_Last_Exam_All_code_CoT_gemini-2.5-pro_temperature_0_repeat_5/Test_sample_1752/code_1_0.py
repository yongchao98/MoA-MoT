import math

# Step 1: Define initial variables
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
drug_admin_concentration = 1  # mg/ml
milk_volume_per_day = 500  # ml/day

# Step 2: Calculate total daily fluid requirement using Holliday-Segar method
if weight_kg <= 10:
    total_daily_fluid = weight_kg * 100
    calc_str = f"{weight_kg} * 100"
elif weight_kg <= 20:
    total_daily_fluid = (10 * 100) + ((weight_kg - 10) * 50)
    calc_str = f"(10 * 100) + ({weight_kg - 10} * 50)"
else:
    total_daily_fluid = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    calc_str = f"(10 * 100) + (10 * 50) + ({weight_kg - 20} * 20)"

print(f"1. Calculate Total Daily Fluid Requirement (Holliday-Segar):")
print(f"   Calculation: {calc_str} = {total_daily_fluid} ml/day\n")

# Step 3: Calculate the daily volume of the chemotherapy drug
daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
daily_drug_volume_ml = daily_drug_dose_mg / drug_admin_concentration

print(f"2. Calculate Daily Drug Volume:")
print(f"   Daily drug dose (mg): {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg/day")
print(f"   Daily drug volume (ml): {daily_drug_dose_mg} mg / {drug_admin_concentration} mg/ml = {daily_drug_volume_ml} ml/day\n")

# Step 4: Calculate total volume from other sources
total_other_fluids = daily_drug_volume_ml + milk_volume_per_day

print(f"3. Calculate Total Other Fluid Sources:")
print(f"   Calculation: {daily_drug_volume_ml} ml (drug) + {milk_volume_per_day} ml (milk) = {total_other_fluids} ml/day\n")

# Step 5: Calculate the remaining volume for maintenance fluid
maintenance_fluid_volume_day = total_daily_fluid - total_other_fluids

print(f"4. Calculate Required Maintenance Fluid Volume:")
print(f"   Calculation: {total_daily_fluid} ml (total) - {total_other_fluids} ml (other) = {maintenance_fluid_volume_day} ml/day\n")

# Step 6: Calculate the hourly rate and round it
maintenance_fluid_rate_hr = maintenance_fluid_volume_day / 24
rounded_rate = round(maintenance_fluid_rate_hr)

print(f"5. Calculate Final Maintenance Fluid Rate:")
print(f"   Hourly rate: {maintenance_fluid_volume_day} ml / 24 hr = {maintenance_fluid_rate_hr:.1f} ml/hr")
print(f"   Rounded rate: {rounded_rate} ml/hr")

print("\n<<<43>>>")