import math

# Patient and drug information
weight_kg = 22
bsa_m2 = 0.8
chemo_dose_per_m2_per_day = 25  # mg/m^2/day
enteral_feeding_ml_day = 500  # ml/day
chemo_admin_concentration_mg_ml = 1  # mg/ml

# Step 1: Calculate total daily fluid requirement (Holliday-Segar method)
# For a 22 kg child:
# 10 kg * 100 ml/kg = 1000 ml
# 10 kg * 50 ml/kg  = 500 ml
# 2 kg * 20 ml/kg   = 40 ml
# Total = 1540 ml/day
daily_fluid_requirement_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
print(f"Step 1: Calculate total daily fluid requirement (Holliday-Segar).")
print(f"   - For a {weight_kg} kg child, the requirement is (10 * 100) + (10 * 50) + (({weight_kg} - 20) * 20) = {daily_fluid_requirement_ml} ml/day.")
print("-" * 60)

# Step 2: Calculate total daily volume of chemotherapy
total_daily_chemo_dose_mg = chemo_dose_per_m2_per_day * bsa_m2
chemo_volume_ml_day = total_daily_chemo_dose_mg / chemo_admin_concentration_mg_ml
print(f"Step 2: Calculate the daily volume of the chemotherapy drug.")
print(f"   - Daily Dose (mg): {chemo_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_chemo_dose_mg} mg")
print(f"   - Daily Volume (ml): {total_daily_chemo_dose_mg} mg / {chemo_admin_concentration_mg_ml} mg/ml = {chemo_volume_ml_day} ml")
print("-" * 60)

# Step 3: Sum all non-maintenance fluid volumes
total_other_fluids_ml_day = chemo_volume_ml_day + enteral_feeding_ml_day
print(f"Step 3: Sum all other fluid sources.")
print(f"   - Total from other sources: {chemo_volume_ml_day} ml (chemo) + {enteral_feeding_ml_day} ml (feeding) = {total_other_fluids_ml_day} ml/day.")
print("-" * 60)

# Step 4: Calculate the remaining maintenance fluid volume
remaining_maintenance_volume_ml_day = daily_fluid_requirement_ml - total_other_fluids_ml_day
print(f"Step 4: Calculate the remaining volume needed for maintenance fluids.")
print(f"   - Remaining volume: {daily_fluid_requirement_ml} ml (total) - {total_other_fluids_ml_day} ml (other sources) = {remaining_maintenance_volume_ml_day} ml/day.")
print("-" * 60)

# Step 5: Calculate the infusion rate in ml/hr and round
infusion_rate_ml_hr = remaining_maintenance_volume_ml_day / 24
rounded_infusion_rate = round(infusion_rate_ml_hr)
print(f"Step 5: Calculate the final infusion rate and round to the nearest whole number.")
print(f"   - Hourly rate equation: ( ( (10 * 100) + (10 * 50) + ( (22 - 20) * 20) ) - ( ( ({chemo_dose_per_m2_per_day} * {bsa_m2}) / {chemo_admin_concentration_mg_ml} ) + {enteral_feeding_ml_day} ) ) / 24")
print(f"   - Unrounded rate: {remaining_maintenance_volume_ml_day} ml / 24 hours = {infusion_rate_ml_hr:.2f} ml/hr")
print(f"   - Final rounded rate: {rounded_infusion_rate} ml/hr")
print("-" * 60)

print(f"<<<{rounded_infusion_rate}>>>")