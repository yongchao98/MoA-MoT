import math

# Step 1: Define the given values
weight_kg = 22
bsa_m2 = 0.8
drug_dose_mg_m2_day = 25
admin_concentration_mg_ml = 1.0
enteral_feeding_ml_day = 500

# Step 2: Calculate the total daily maintenance fluid requirement using the Holliday-Segar method
# For a child > 20 kg: 100ml/kg for first 10kg, 50ml/kg for next 10kg, 20ml/kg for the rest
fluid_first_10kg = 10 * 100
fluid_next_10kg = 10 * 50
fluid_remaining_kg = (weight_kg - 20) * 20
total_maintenance_fluid_ml_day = fluid_first_10kg + fluid_next_10kg + fluid_remaining_kg

# Step 3: Calculate the volume of the chemotherapy drug administered per day
daily_drug_dose_mg = drug_dose_mg_m2_day * bsa_m2
daily_drug_volume_ml = daily_drug_dose_mg / admin_concentration_mg_ml

# Step 4: Calculate the total volume of fluid from other sources
total_other_fluids_ml_day = daily_drug_volume_ml + enteral_feeding_ml_day

# Step 5: Calculate the remaining volume to be given as maintenance fluid
remaining_iv_fluid_ml_day = total_maintenance_fluid_ml_day - total_other_fluids_ml_day

# Step 6: Calculate the hourly rate for the maintenance fluid and round to the nearest whole number
iv_rate_ml_hr = remaining_iv_fluid_ml_day / 24
# Standard rounding: add 0.5 and take the integer part to round .5 up
final_rate_ml_hr = int(iv_rate_ml_hr + 0.5)

# --- Output the results step-by-step ---
print("--- Calculation of Maintenance Fluid Rate ---")
print(f"\n1. Total Daily Fluid Requirement (Holliday-Segar for {weight_kg} kg):")
print(f"   ({10} * {100}) + ({10} * {50}) + (({weight_kg} - 20) * {20}) = {int(total_maintenance_fluid_ml_day)} ml/day")

print(f"\n2. Fluid from Other Sources:")
print(f"   Chemotherapy Volume: ({drug_dose_mg_m2_day} mg/m² * {bsa_m2} m²) / {int(admin_concentration_mg_ml)} mg/ml = {int(daily_drug_volume_ml)} ml/day")
print(f"   Enteral Feeding Volume: {enteral_feeding_ml_day} ml/day")
print(f"   Total from other sources: {int(daily_drug_volume_ml)} ml + {enteral_feeding_ml_day} ml = {int(total_other_fluids_ml_day)} ml/day")

print(f"\n3. Remaining IV Fluid Calculation:")
print(f"   Remaining Volume: {int(total_maintenance_fluid_ml_day)} ml - {int(total_other_fluids_ml_day)} ml = {int(remaining_iv_fluid_ml_day)} ml/day")

print(f"\n4. Final IV Rate Calculation:")
print(f"   Hourly Rate: {int(remaining_iv_fluid_ml_day)} ml / 24 hr = {iv_rate_ml_hr:.1f} ml/hr")
print(f"   Rounded Rate: {final_rate_ml_hr} ml/hr")

print(f"\n<<< {final_rate_ml_hr} >>>")