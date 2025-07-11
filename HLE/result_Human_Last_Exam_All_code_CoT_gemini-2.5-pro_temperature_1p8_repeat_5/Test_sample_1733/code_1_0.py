import math

# Patient Information
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Formula: Weight (kg) * Dose (mL/kg)
resuscitation_dose_ml_per_kg = 30
initial_resuscitation_volume = weight_kg * resuscitation_dose_ml_per_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg (11-20 kg)
# 20 mL/kg for weight > 20 kg
if weight_kg <= 10:
    standard_maintenance_volume = weight_kg * 100
elif weight_kg <= 20:
    standard_maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    standard_maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Adjust for mechanical ventilation (reduces insensible losses, typically by 25%)
# A common multiplier is 0.75
ventilation_adjustment_factor = 0.75
adjusted_maintenance_volume = standard_maintenance_volume * ventilation_adjustment_factor


# --- Calculation 3: Total Deficit Replacement Fluid Volume ---
# Formula: Weight (kg) * Dehydration Percentage * 1000 mL/L
dehydration_percentage = 0.10 # 10%
# 1 kg weight loss from dehydration is ~1 L (1000 mL) of fluid deficit
total_deficit_volume = weight_kg * dehydration_percentage * 1000


# --- Final Output ---
# Print the results separated by commas
# The values are rounded to the nearest whole number for clinical clarity.
print(f"{math.ceil(initial_resuscitation_volume)},{math.ceil(adjusted_maintenance_volume)},{math.ceil(total_deficit_volume)}")

<<<360,825,1200>>>