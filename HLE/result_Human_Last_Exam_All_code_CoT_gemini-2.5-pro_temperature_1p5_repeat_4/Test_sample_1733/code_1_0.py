import math

# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Phase 1 (Initial Resuscitation) ---
# An initial resuscitation bolus of 30 mL/kg of isotonic saline.
resuscitation_bolus_per_kg = 30
total_resuscitation_volume = weight_kg * resuscitation_bolus_per_kg

# --- Calculation 2: Phase 2 (Daily Maintenance Fluids) ---
# Maintenance fluids calculated using the Holliday-Segar method for 24 hours.
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg (11-20 kg)
# 20 mL/kg for the weight over 20 kg
daily_maintenance_volume = 0
if weight_kg <= 10:
    daily_maintenance_volume = weight_kg * 100
elif weight_kg <= 20:
    daily_maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    daily_maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Note: The additional fluids (antibiotics, milk) and ventilation status are clinical
# considerations but do not alter the calculation of the core maintenance volume
# based on the provided Holliday-Segar formula.

# --- Calculation 3: Phase 3 (Deficit Replacement Fluids) ---
# Deficit replacement fluids estimated at 10% of body weight.
# 1 kg of weight loss from dehydration equals 1000 mL of fluid deficit.
dehydration_percentage = 0.10
fluid_deficit_in_kg = weight_kg * dehydration_percentage
total_deficit_volume = fluid_deficit_in_kg * 1000

# --- Print the final results in the required format ---
# The result format is three exact numbers separated by ","
print(f"{math.ceil(total_resuscitation_volume)},{math.ceil(daily_maintenance_volume)},{math.ceil(total_deficit_volume)}")