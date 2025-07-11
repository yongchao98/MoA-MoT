import math

# Patient data
weight_kg = 12
dehydration_percent = 10
iv_antibiotics_ml_day = 60
enteral_nutrition_ml_day = 100

# --- Calculation Step 1: Initial Resuscitation Bolus ---
# Formula: 30 mL/kg
resuscitation_bolus_ml = 30 * weight_kg

# --- Calculation Step 2: Daily Maintenance Fluid ---
# Holliday-Segar Method
if weight_kg <= 10:
    maintenance_fluid_ml_day = 100 * weight_kg
elif weight_kg <= 20:
    maintenance_fluid_ml_day = (100 * 10) + (50 * (weight_kg - 10))
else:
    maintenance_fluid_ml_day = (100 * 10) + (50 * 10) + (20 * (weight_kg - 20))

# Adjustment for mechanical ventilation (reduce by 25%)
# This reduces insensible losses.
adjusted_maintenance_fluid_ml_day = maintenance_fluid_ml_day * 0.75

# --- Calculation Step 3: Deficit Replacement Fluid ---
# Formula: Deficit (L) = % dehydration * weight (kg) / 100
# Deficit (mL) = % dehydration * weight (kg) * 10
deficit_volume_ml = dehydration_percent * weight_kg * 10

# --- Final Output ---
# The final answer requires three exact numbers separated by ",".
# We will use integer values as is common practice for these calculations.
resuscitation_vol = int(resuscitation_bolus_ml)
maintenance_vol = int(adjusted_maintenance_fluid_ml_day)
deficit_vol = int(deficit_volume_ml)

print(f"{resuscitation_vol},{maintenance_vol},{deficit_vol}")

# Final Answer Block
print(f"<<<{resuscitation_vol},{maintenance_vol},{deficit_vol}>>>")