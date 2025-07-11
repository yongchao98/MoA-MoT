import math

# Patient variables
weight_kg = 12

# Phase 1: Resuscitation Bolus Calculation
resuscitation_rate_ml_per_kg = 30
resuscitation_volume_ml = weight_kg * resuscitation_rate_ml_per_kg

# Phase 2: Maintenance Fluid Calculation (Holliday-Segar Method)
if weight_kg <= 10:
    base_maintenance_volume_ml = weight_kg * 100
else:
    base_maintenance_volume_ml = (10 * 100) + ((weight_kg - 10) * 50)

# Adjust for mechanical ventilation (reduced insensible losses, typically 75-80% of standard)
ventilation_adjustment_factor = 0.80
adjusted_maintenance_volume_ml = base_maintenance_volume_ml * ventilation_adjustment_factor

# Phase 3: Deficit Replacement Calculation
dehydration_percent = 0.10  # 10%
# 1 kg of weight loss due to dehydration is equivalent to 1 L (1000 mL) of fluid loss.
deficit_volume_ml = weight_kg * dehydration_percent * 1000

# Print the results as integers separated by commas
# The print statement explicitly casts each value to an integer for a clean output.
print(f"{int(resuscitation_volume_ml)},{int(adjusted_maintenance_volume_ml)},{int(deficit_volume_ml)}")
