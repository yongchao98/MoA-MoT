import math

# Patient parameters
weight_kg = 12.0
dehydration_percentage = 0.10  # 10%
resuscitation_bolus_factor_ml_per_kg = 30.0

# --- Calculation for Part 1: Initial Resuscitation Volume ---
# This is a direct calculation: weight * bolus factor.
resuscitation_volume = weight_kg * resuscitation_bolus_factor_ml_per_kg

# --- Calculation for Part 2: Daily Maintenance Fluid Volume ---
# Using the Holliday-Segar method:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg (from 11 to 20 kg)
# 20 mL/kg for weight above 20 kg
if weight_kg <= 10:
    holliday_segar_volume = weight_kg * 100
elif weight_kg <= 20:
    holliday_segar_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    holliday_segar_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Adjusting for mechanical ventilation.
# Standard practice is to reduce maintenance fluids, often to 2/3 of the calculated value.
ventilation_adjustment_factor = 2 / 3
adjusted_maintenance_volume = holliday_segar_volume * ventilation_adjustment_factor

# --- Calculation for Part 3: Total Deficit Replacement Fluid Volume ---
# Calculated as 10% of body weight. 1 kg of body weight loss equals 1000 mL of fluid.
deficit_volume = weight_kg * dehydration_percentage * 1000

# --- Final Output ---
# Print the three calculated values, rounded to the nearest whole number,
# separated by commas as requested.
print(f"{round(resuscitation_volume)},{round(adjusted_maintenance_volume)},{round(deficit_volume)}")