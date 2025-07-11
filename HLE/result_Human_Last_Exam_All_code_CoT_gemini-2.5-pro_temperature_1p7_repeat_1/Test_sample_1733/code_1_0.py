import math

# Patient Information
weight_kg = 12
dehydration_percent = 10
bolus_rate_ml_per_kg = 30

# --- Step 1: Calculate Initial Resuscitation Volume ---
# This is an initial fluid bolus given over 1 hour.
# Formula: Volume = Weight (kg) * Bolus Rate (mL/kg)
resuscitation_volume_ml = weight_kg * bolus_rate_ml_per_kg
print(f"Step 1: The resuscitation volume is calculated as {weight_kg} kg * {bolus_rate_ml_per_kg} mL/kg = {int(resuscitation_volume_ml)} mL.")


# --- Step 2: Calculate Daily Maintenance Fluid Volume ---
# Using the Holliday-Segar Method for a 12 kg child.
# - 100 mL/kg for the first 10 kg
# - 50 mL/kg for the next 10 kg (i.e., for weight from 11-20 kg)
maintenance_for_first_10kg = 10 * 100
maintenance_for_next_2kg = (weight_kg - 10) * 50
maintenance_volume_unadjusted = maintenance_for_first_10kg + maintenance_for_next_2kg
print(f"Step 2: The unadjusted daily maintenance fluid volume is ({10} kg * 100 mL/kg) + ({(weight_kg - 10)} kg * 50 mL/kg) = {int(maintenance_volume_unadjusted)} mL.")

# Adjust for mechanical ventilation, which reduces insensible losses by approx. 25%.
# The volumes for antibiotics and milk are accounted for when setting the final IV fluid *rate*,
# but they do not change the total calculated *requirement*.
ventilation_reduction_factor = 0.25
maintenance_volume_adjusted = maintenance_volume_unadjusted * (1 - ventilation_reduction_factor)
print(f"Step 2 (cont.): After adjusting for mechanical ventilation (25% reduction), the total daily maintenance volume is {int(maintenance_volume_unadjusted)} mL * (1 - {ventilation_reduction_factor}) = {int(maintenance_volume_adjusted)} mL.")


# --- Step 3: Calculate Total Deficit Replacement Volume ---
# The deficit is estimated at 10% of body weight.
# Formula: Deficit (mL) = % Dehydration * Weight (kg) * 10
# This total volume is given over 48 hours.
deficit_volume_ml = dehydration_percent * weight_kg * 10
print(f"Step 3: The total deficit replacement volume is {dehydration_percent}% * {weight_kg} kg * 10 = {int(deficit_volume_ml)} mL.")


# --- Final Answer ---
# The final answer is the three calculated numbers separated by a comma.
print("\nFinal calculated values (Resuscitation, Maintenance, Deficit):")
print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_adjusted)},{int(deficit_volume_ml)}")
<<<360,825,1200>>>