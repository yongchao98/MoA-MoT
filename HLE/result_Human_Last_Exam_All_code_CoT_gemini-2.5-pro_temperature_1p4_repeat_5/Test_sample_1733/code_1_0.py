import math

# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Formula: 30 mL/kg
bolus_rate_ml_per_kg = 30
resuscitation_volume = bolus_rate_ml_per_kg * weight_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# Formula: 100 mL/kg for the first 10 kg, 50 mL/kg for the next 10-20 kg.
maintenance_volume = 0
if weight_kg <= 10:
    maintenance_volume = 100 * weight_kg
else:
    # Fluid for the first 10 kg
    maintenance_volume += 100 * 10
    # Fluid for the remaining weight (12kg - 10kg = 2kg)
    remaining_weight = weight_kg - 10
    maintenance_volume += 50 * remaining_weight

# --- Calculation 3: Total Deficit Replacement Fluid Volume ---
# Formula: 10% of body weight, where 1 kg weight loss = 1000 mL fluid deficit.
deficit_percentage = 0.10
fluid_deficit_in_kg = weight_kg * deficit_percentage
# Convert kg to mL (1 kg -> 1000 mL)
deficit_replacement_volume = fluid_deficit_in_kg * 1000

# Print the final calculated values separated by a comma.
# The results are cast to integers as fluid calculations are typically rounded to whole numbers.
print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_replacement_volume)}")