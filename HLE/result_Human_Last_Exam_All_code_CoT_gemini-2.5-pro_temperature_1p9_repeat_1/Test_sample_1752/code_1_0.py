import math

# Patient and Drug Information
weight_kg = 22
bsa_m2 = 0.8
chemo_dose_per_m2_day = 25  # mg/m²/day
chemo_admin_concentration_mg_ml = 1.0  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily fluid requirement (Holliday-Segar method)
# This represents 100% of the maintenance fluid allowance.
if weight_kg <= 10:
    total_maintenance_fluid_day = weight_kg * 100
elif weight_kg <= 20:
    total_maintenance_fluid_day = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_maintenance_fluid_day = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Step 2: Calculate the volume of the chemotherapy drug administration
daily_chemo_dose_mg = chemo_dose_per_m2_day * bsa_m2
chemo_volume_day = daily_chemo_dose_mg / chemo_admin_concentration_mg_ml

# Step 3: Sum all other fluid sources
total_other_fluids_day = enteral_feeding_ml_day + chemo_volume_day

# Step 4: Calculate the remaining volume for IV maintenance fluid
remaining_iv_fluid_day = total_maintenance_fluid_day - total_other_fluids_day

# Step 5: Calculate the IV fluid rate in ml/hr and round it
iv_rate_ml_hr = remaining_iv_fluid_day / 24
final_rate = round(iv_rate_ml_hr)

# Print the calculation steps and the final equation
print("--- Calculation of Maintenance Fluid Rate ---")
print(f"1. Total daily fluid requirement for a {weight_kg} kg child: {total_maintenance_fluid_day} ml/day")
print(f"2. Daily chemotherapy fluid volume: {chemo_volume_day} ml/day")
print(f"3. Daily enteral feeding volume: {enteral_feeding_ml_day} ml/day")
print("\nFinal Rate Calculation:")
print("The rate (ml/hr) is calculated as (Total Daily Fluid - Other Fluids) / 24 hours.")
print(f"( {total_maintenance_fluid_day} ml - ({enteral_feeding_ml_day} ml + {chemo_volume_day} ml) ) / 24 hr = {iv_rate_ml_hr:.1f} ml/hr")
print(f"\nRounded to the nearest whole number, the maintenance fluid rate should be set to: {final_rate} ml/hr")
print("<<<43>>>")