import math

# --- Given Information ---
weight_kg = 22
bsa_m2 = 0.8
drug_dose_per_m2_per_day = 25  # mg/m^2/day
drug_admin_concentration = 1  # mg/ml
milk_volume_per_day = 500  # ml/day
hours_in_day = 24

# --- Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar Method) ---
print("Step 1: Calculate Total Daily Fluid Requirement (Holliday-Segar Method)")
daily_fluid_req = 0
weight_remaining = weight_kg

if weight_remaining > 20:
    fluid_for_over_20kg = (weight_remaining - 20) * 20
    daily_fluid_req += fluid_for_over_20kg
    print(f"For weight over 20 kg: ({weight_kg} kg - 20 kg) * 20 ml/kg = {fluid_for_over_20kg} ml")
    weight_remaining = 20

if weight_remaining > 10:
    fluid_for_10_to_20kg = (weight_remaining - 10) * 50
    daily_fluid_req += fluid_for_10_to_20kg
    print(f"For weight between 10-20 kg: (20 kg - 10 kg) * 50 ml/kg = {fluid_for_10_to_20kg} ml")
    weight_remaining = 10

if weight_remaining > 0:
    fluid_for_first_10kg = weight_remaining * 100
    daily_fluid_req += fluid_for_first_10kg
    print(f"For first 10 kg: 10 kg * 100 ml/kg = {fluid_for_first_10kg} ml")

print(f"Total Daily Fluid Requirement = 1000 + 500 + {fluid_for_over_20kg} = {daily_fluid_req} ml/day\n")

# --- Step 2: Calculate Daily Drug Volume ---
print("Step 2: Calculate Daily Drug Volume")
total_daily_drug_dose_mg = drug_dose_per_m2_per_day * bsa_m2
print(f"Total daily drug dose (mg) = {drug_dose_per_m2_per_day} mg/m² * {bsa_m2} m² = {total_daily_drug_dose_mg} mg")

drug_volume_per_day = total_daily_drug_dose_mg / drug_admin_concentration
print(f"Total daily drug volume (ml) = {total_daily_drug_dose_mg} mg / {drug_admin_concentration} mg/ml = {drug_volume_per_day} ml\n")


# --- Step 3: Calculate Total Intake from Other Sources ---
print("Step 3: Calculate Total Intake from Other Sources")
total_other_fluids = drug_volume_per_day + milk_volume_per_day
print(f"Total other fluids = {drug_volume_per_day} ml (drug) + {milk_volume_per_day} ml (milk) = {total_other_fluids} ml\n")

# --- Step 4: Calculate Remaining Maintenance Fluid ---
print("Step 4: Calculate Remaining Maintenance Fluid")
remaining_fluid_volume = daily_fluid_req - total_other_fluids
print(f"Remaining fluid = {daily_fluid_req} ml (total) - {total_other_fluids} ml (other) = {remaining_fluid_volume} ml\n")

# --- Step 5: Calculate Hourly Rate ---
print("Step 5: Calculate and Round the Hourly Rate")
hourly_rate = remaining_fluid_volume / hours_in_day
print(f"Unrounded hourly rate = {remaining_fluid_volume} ml / {hours_in_day} hours = {hourly_rate:.2f} ml/hr")

# Rounding to the nearest whole number (.5 rounds up)
final_rate = int(hourly_rate + 0.5)
print(f"\nFinal rounded maintenance fluid rate is calculated as: ({daily_fluid_req} - ({drug_volume_per_day} + {milk_volume_per_day})) / {hours_in_day} = {final_rate} ml/hr")
print(f"The final rate for the maintenance fluid is {final_rate} ml/hr.")
print("<<<43>>>")