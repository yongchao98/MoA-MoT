import math

# --- Patient Information ---
weight_kg = 12
dehydration_percentage = 0.10

# --- Phase 1: Initial Resuscitation Bolus ---
# Calculation: 30 mL/kg of isotonic saline
resuscitation_bolus_rate_ml_per_kg = 30
resuscitation_volume = weight_kg * resuscitation_bolus_rate_ml_per_kg

# --- Phase 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
# Calculation for a 24-hour period:
# - 100 mL/kg for the first 10 kg
# - 50 mL/kg for the next 10 kg (from 11 to 20 kg)
# - 20 mL/kg for weight over 20 kg
maintenance_volume = 0
if weight_kg <= 10:
    maintenance_volume = weight_kg * 100
elif weight_kg <= 20:
    maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# --- Phase 3: Deficit Replacement Fluids ---
# Calculation: Deficit is 10% of body weight.
# 1 kg of weight is assumed to be equivalent to 1000 mL of fluid.
deficit_in_kg = weight_kg * dehydration_percentage
deficit_volume = deficit_in_kg * 1000

# --- Output the results ---
# The final output provides the three calculated volumes:
# 1. Total resuscitation volume (mL)
# 2. Daily maintenance volume (mL)
# 3. Total deficit replacement volume (mL)
# Note: Additional fluids from antibiotics (60 mL/day) and enteral nutrition (100 mL/day)
# are important for total fluid management but are not part of these specific phase calculations.
# Similarly, mechanical ventilation might lead a clinician to adjust the maintenance fluid rate,
# but the standard calculated volume is presented here.
print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")
