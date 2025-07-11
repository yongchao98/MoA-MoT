import math

# Patient Information
weight_kg = 12
dehydration_percentage = 0.10  # 10%

# --- 1. Calculate the initial resuscitation bolus volume ---
# Formula: 30 mL/kg
bolus_volume = 30 * weight_kg

# --- 2. Calculate the daily maintenance fluid volume for 24 hours ---
# Holliday-Segar Method:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg (11-20 kg)
# 20 mL/kg for weight over 20 kg
if weight_kg <= 10:
    maintenance_base = weight_kg * 100
else: # For a 12 kg child
    maintenance_base = (10 * 100) + ((weight_kg - 10) * 50)

# Adjust for mechanical ventilation (reduced insensible losses, standard 25% reduction)
ventilation_adjustment_factor = 0.75 # 100% - 25%
maintenance_adjusted = maintenance_base * ventilation_adjustment_factor

# --- 3. Calculate the total deficit replacement fluid volume ---
# Deficit is 10% of body weight. 1 kg of weight loss ~ 1000 mL of fluid.
deficit_volume = weight_kg * dehydration_percentage * 1000

# --- Print the final results in the specified format ---
# The results are cast to integers as fluid volumes are typically managed in whole numbers.
print(f"{int(bolus_volume)},{int(maintenance_adjusted)},{int(deficit_volume)}")