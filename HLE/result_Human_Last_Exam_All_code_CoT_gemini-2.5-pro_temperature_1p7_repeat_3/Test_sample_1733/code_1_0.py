import math

# Patient data
weight_kg = 12
resuscitation_bolus_per_kg = 30  # mL/kg
deficit_percentage = 0.10  # 10%
ventilation_adjustment_factor = 0.75 # Standard reduction for mechanical ventilation

# --- Calculation Step 1: Initial Resuscitation Fluid ---
# Formula: Weight (kg) * Bolus (mL/kg)
resuscitation_volume = weight_kg * resuscitation_bolus_per_kg

# --- Calculation Step 2: Daily Maintenance Fluid (Holliday-Segar Method) ---
# For the first 10 kg: 100 mL/kg
# For weight > 10 kg up to 20 kg: 50 mL/kg for each kg over 10
if weight_kg <= 10:
    standard_maintenance_volume = weight_kg * 100
else:
    standard_maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)

# Adjusting for mechanical ventilation
adjusted_maintenance_volume = standard_maintenance_volume * ventilation_adjustment_factor

# --- Calculation Step 3: Total Deficit Replacement Fluid ---
# Formula: Weight (kg) * 1000 (g/kg) * Deficit Percentage
# 1g of body weight loss is assumed to be 1 mL of fluid deficit
deficit_volume = weight_kg * 1000 * deficit_percentage

# --- Final Output ---
# Printing the three calculated values, rounded to the nearest whole number
# The format is three exact numbers separated by ","
result_1 = math.ceil(resuscitation_volume)
result_2 = math.ceil(adjusted_maintenance_volume)
result_3 = math.ceil(deficit_volume)

print(f"{result_1},{result_2},{result_3}")

# The final answer format
print(f'<<<{result_1},{result_2},{result_3}>>>')