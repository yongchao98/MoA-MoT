# Patient's weight in kg
weight_kg = 12

# --- 1. Initial Resuscitation Volume ---
# Formula: 30 mL/kg
bolus_rate_ml_per_kg = 30
resuscitation_volume = weight_kg * bolus_rate_ml_per_kg

# --- 2. Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# Formula:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10-20 kg
# 20 mL/kg for weight > 20 kg
maintenance_volume_unadjusted = 0
if weight_kg <= 10:
    maintenance_volume_unadjusted = weight_kg * 100
else:
    maintenance_volume_unadjusted = 10 * 100  # For the first 10 kg
    maintenance_volume_unadjusted += (weight_kg - 10) * 50 # For the weight above 10 kg

# Adjusting for mechanical ventilation (reducing fluids to 80%)
ventilation_adjustment_factor = 0.8
maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_adjustment_factor

# --- 3. Total Deficit Replacement Volume ---
# Formula: 10% of body weight (in kg), converted to mL (1 kg fluid loss ~ 1 L)
dehydration_percentage = 0.10
# Deficit in Liters = weight * percentage
deficit_liters = weight_kg * dehydration_percentage
# Convert Liters to milliliters
deficit_volume_ml = deficit_liters * 1000

# --- Final Output ---
# Print the results in the format: resuscitation,maintenance,deficit
# The values are converted to integers as fluid volumes are typically managed in whole numbers.
print(f"{int(resuscitation_volume)},{int(maintenance_volume_adjusted)},{int(deficit_volume_ml)}")