import math

# Patient's weight in kg
weight_kg = 12

# --- 1. Initial Resuscitation Volume ---
# Formula: 30 mL/kg
bolus_dose_ml_per_kg = 30
resuscitation_volume = weight_kg * bolus_dose_ml_per_kg

print(f"1. Initial Resuscitation Volume: {weight_kg} kg * {bolus_dose_ml_per_kg} mL/kg = {int(resuscitation_volume)} mL")

# --- 2. Daily Maintenance Fluid Volume (Holliday-Segar method) ---
# Formula: 100 mL/kg for first 10 kg, 50 mL/kg for next 10-20 kg.
if weight_kg <= 10:
    maintenance_volume = weight_kg * 100
    print(f"2. Daily Maintenance Volume: {weight_kg} kg * 100 mL/kg = {int(maintenance_volume)} mL")
elif weight_kg <= 20:
    maintenance_volume_first_10kg = 10 * 100
    maintenance_volume_next_kg = (weight_kg - 10) * 50
    maintenance_volume = maintenance_volume_first_10kg + maintenance_volume_next_kg
    print(f"2. Daily Maintenance Volume: ({10} kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg) = {int(maintenance_volume)} mL")
else:
    maintenance_volume_first_10kg = 10 * 100
    maintenance_volume_next_10kg = 10 * 50
    maintenance_volume_remaining_kg = (weight_kg - 20) * 20
    maintenance_volume = maintenance_volume_first_10kg + maintenance_volume_next_10kg + maintenance_volume_remaining_kg
    print(f"2. Daily Maintenance Volume: ({10} kg * 100 mL/kg) + ({10} kg * 50 mL/kg) + ({weight_kg - 20} kg * 20 mL/kg) = {int(maintenance_volume)} mL")

# --- 3. Total Deficit Replacement Fluid Volume ---
# Formula: 10% deficit is equivalent to 100 mL/kg
deficit_percentage = 10
deficit_ml_per_kg = deficit_percentage * 10 # 1% deficit = 10 mL/kg
deficit_volume = weight_kg * deficit_ml_per_kg

print(f"3. Total Deficit Replacement Volume: {weight_kg} kg * {deficit_ml_per_kg} mL/kg (for {deficit_percentage}% deficit) = {int(deficit_volume)} mL")

# --- Final Answer in requested format ---
print(f"\nFinal Answer:")
print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")
