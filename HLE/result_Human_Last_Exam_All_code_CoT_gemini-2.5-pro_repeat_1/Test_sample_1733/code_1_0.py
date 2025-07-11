import math

# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Formula: 30 mL/kg
resuscitation_bolus_per_kg = 30
resuscitation_volume = resuscitation_bolus_per_kg * weight_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# For the first 10 kg: 100 mL/kg
# For 10.1 to 20 kg: 50 mL/kg
# For > 20 kg: 20 mL/kg
if weight_kg <= 10:
    maintenance_volume = weight_kg * 100
elif weight_kg <= 20:
    maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
else:
    maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# --- Calculation 3: Total Deficit Replacement Volume ---
# Formula: 10% of body weight in kg, where 1 kg deficit = 1000 mL
dehydration_percentage = 0.10
deficit_volume = weight_kg * dehydration_percentage * 1000

# --- Print the Final Results ---
# The results are cast to integers as fluid volumes are typically managed in whole numbers.
print(f"1. Resuscitation Volume: {int(resuscitation_volume)} mL")
print(f"2. Daily Maintenance Volume: {int(maintenance_volume)} mL")
print(f"3. Total Deficit Volume: {int(deficit_volume)} mL")
print("\nFinal Answer:")
print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")