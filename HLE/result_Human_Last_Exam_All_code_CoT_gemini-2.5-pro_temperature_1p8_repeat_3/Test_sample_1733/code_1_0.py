# Patient data
weight_kg = 12
dehydration_percent = 0.10  # 10%

# Phase 1: Initial Resuscitation Bolus
# Formula: 30 mL/kg * weight
resuscitation_rate_ml_per_kg = 30
resuscitation_volume_ml = resuscitation_rate_ml_per_kg * weight_kg

# Phase 2: Daily Maintenance Fluids (Holliday-Segar Method)
# For the first 10 kg: 100 mL/kg
# For 11-20 kg: 50 mL/kg
# For >20 kg: 20 mL/kg
if weight_kg <= 10:
    maintenance_volume_ml_24h = weight_kg * 100
else:
    # Patient weight is 12 kg
    # Fluid for the first 10 kg
    fluid_for_first_10kg = 10 * 100
    # Fluid for the weight above 10 kg
    fluid_for_remaining_weight = (weight_kg - 10) * 50
    maintenance_volume_ml_24h = fluid_for_first_10kg + fluid_for_remaining_weight

# Phase 3: Deficit Replacement Fluids
# Formula: weight (kg) * dehydration (%) * 1000 (to convert L to mL)
# This represents a 10% loss of body weight in fluid. 1kg of weight loss = 1L of fluid.
deficit_volume_ml = weight_kg * dehydration_percent * 1000

# Print the final results in the required format
# Ensure the outputs are integers
print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_ml_24h)},{int(deficit_volume_ml)}")
<<<360,1100,1200>>>