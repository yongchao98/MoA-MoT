import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
admin_concentration_mg_ml = 1  # mg/ml
enteral_feeding_ml_day = 500

# Step 1: Calculate total daily drug dose in mg
daily_dose_mg = drug_dose_per_m2_per_day * bsa_m2

# Step 2: Calculate the volume of the drug infusion per day
daily_chemo_volume_ml = daily_dose_mg / admin_concentration_mg_ml

# Step 3: Calculate total daily maintenance fluid requirement (Holliday-Segar method)
total_daily_fluid_ml = 0
if weight_kg <= 10:
    total_daily_fluid_ml = weight_kg * 100
elif weight_kg <= 20:
    total_daily_fluid_ml = (10 * 100) + ((weight_kg - 10) * 50)
else:
    total_daily_fluid_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Step 4: Calculate total fluid from non-maintenance sources
fluid_from_other_sources_ml = enteral_feeding_ml_day + daily_chemo_volume_ml

# Step 5: Calculate the remaining volume needed for IV maintenance fluid
remaining_iv_fluid_ml_day = total_daily_fluid_ml - fluid_from_other_sources_ml

# Step 6: Calculate the rate in ml/hr
maintenance_rate_ml_hr = remaining_iv_fluid_ml_day / 24

# Step 7: Round to the nearest whole number
# We use int(x + 0.5) for standard rounding half up
final_rate = int(maintenance_rate_ml_hr + 0.5)

# Print the final equation with all the numbers
print("Equation to calculate the maintenance fluid rate (ml/hr):")
print(f"(Total Daily Fluid - (Enteral Feeding + Chemo Volume)) / 24 hours")
print(f"({int(total_daily_fluid_ml)} - ({int(enteral_feeding_ml_day)} + {int(daily_chemo_volume_ml)})) / 24 = {final_rate}")

<<<43>>>