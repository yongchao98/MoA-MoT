import math

# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Phase 1: 30 mL/kg bolus
resuscitation_bolus_per_kg = 30
resuscitation_volume = weight_kg * resuscitation_bolus_per_kg

# --- Calculation 2: Daily Maintenance Fluid Volume ---
# Phase 2: Holliday-Segar method for 24 hours
maintenance_volume = 0
if weight_kg <= 10:
    maintenance_volume = weight_kg * 100
elif weight_kg <= 20:
    maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# --- Calculation 3: Total Deficit Replacement Volume ---
# Phase 3: Deficit estimated at 10% of body weight
dehydration_percentage = 0.10
# 1 kg of weight loss equals 1000 mL of fluid deficit
deficit_volume = weight_kg * dehydration_percentage * 1000

# --- Final Output ---
# Print the three values separated by a comma.
# Using math.ceil to ensure whole numbers as partial milliliters are not practical in this context.
print(f"{math.ceil(resuscitation_volume)},{math.ceil(maintenance_volume)},{math.ceil(deficit_volume)}")