# This script calculates fluid requirements for a pediatric patient.

# --- Given Data ---
weight_kg = 12
resuscitation_rate_ml_per_kg = 30
dehydration_percentage = 0.10

# --- Calculation 1: Initial Resuscitation Volume ---
# This is a one-time bolus given over 1 hour.
# Formula: Weight (kg) * Rate (mL/kg)
resuscitation_volume_ml = weight_kg * resuscitation_rate_ml_per_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# This calculates the physiological fluid need for 24 hours.
# Formula:
# - 100 mL/kg for the first 10 kg
# - 50 mL/kg for the next 10 kg (11-20 kg)
# - 20 mL/kg for weight above 20 kg
if weight_kg <= 10:
    maintenance_volume_unadjusted = weight_kg * 100
elif weight_kg <= 20:
    maintenance_volume_unadjusted = (10 * 100) + ((weight_kg - 10) * 50)
else:
    maintenance_volume_unadjusted = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Adjust for mechanical ventilation, which reduces insensible losses.
# A standard reduction of 25% is applied (multiplying by 0.75).
ventilation_reduction_factor = 0.75
maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_reduction_factor

# --- Calculation 3: Total Deficit Replacement Volume ---
# This volume is based on the estimated dehydration level (10% of body weight).
# It will be replaced over 48 hours. The calculation is for the TOTAL volume.
# Formula: Weight (kg) * Dehydration Percentage * 1000 (to convert kg to mL)
deficit_volume_ml = weight_kg * dehydration_percentage * 1000

# --- Final Output ---
# Print the three calculated volumes, rounded to the nearest whole number,
# separated by commas as requested.
print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_adjusted)},{int(deficit_volume_ml)}")