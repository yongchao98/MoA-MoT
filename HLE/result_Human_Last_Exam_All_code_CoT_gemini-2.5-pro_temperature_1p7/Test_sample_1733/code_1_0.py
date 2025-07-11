import math

# Patient Information
weight_kg = 12

# --- Calculations ---

# 1. Initial Resuscitation Volume
bolus_dose_ml_per_kg = 30
resuscitation_volume = weight_kg * bolus_dose_ml_per_kg

# 2. Daily Maintenance Fluid Volume (Holliday-Segar Method)
if weight_kg <= 10:
    maintenance_volume = weight_kg * 100
    maintenance_calculation_str = f"({weight_kg} kg * 100 mL/kg)"
    maintenance_result_str = f"{maintenance_volume} mL"
elif weight_kg <= 20:
    maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    maintenance_calculation_str = f"(10 kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg) = 1000 + {int((weight_kg - 10) * 50)}"
    maintenance_result_str = f"{int(maintenance_volume)} mL"
else:
    maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    maintenance_calculation_str = f"(10 kg * 100 mL/kg) + (10 kg * 50 mL/kg) + ({weight_kg - 20} kg * 20 mL/kg)"
    maintenance_result_str = f"{int(maintenance_volume)} mL"


# 3. Total Deficit Replacement Volume
dehydration_percentage = 0.10  # 10%
deficit_kg = weight_kg * dehydration_percentage
deficit_volume = deficit_kg * 1000  # 1 kg loss = 1000 mL fluid

# --- Output Results ---
print("Calculations:")
print(f"1. Initial Resuscitation Volume: {weight_kg} kg * {bolus_dose_ml_per_kg} mL/kg = {int(resuscitation_volume)} mL")
print(f"2. Daily Maintenance Fluid Volume: {maintenance_calculation_str} = {maintenance_result_str}")
print(f"3. Total Deficit Replacement Volume: {weight_kg} kg * {int(dehydration_percentage*100)}% dehydration * 1000 mL/kg = {int(deficit_volume)} mL")

final_answer = f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}"
print(f"\nFinal Answer (Resuscitation, Maintenance, Deficit):")
print(final_answer)

print(f"<<<{final_answer}>>>")