import math

# --- Patient Data ---
weight_kg = 12
bolus_rate_ml_per_kg = 30
dehydration_percent = 10
# Adjustment factor for mechanical ventilation
ventilation_adjustment_factor = 2/3

# --- 1. Resuscitation Bolus Calculation ---
resuscitation_volume_ml = bolus_rate_ml_per_kg * weight_kg

print("--- 1. Resuscitation Bolus Volume ---")
print(f"The calculation is the bolus rate multiplied by the patient's weight.")
print(f"Formula: {bolus_rate_ml_per_kg} mL/kg * {weight_kg} kg")
print(f"Result: {int(resuscitation_volume_ml)} mL")
print("-" * 40)


# --- 2. Daily Maintenance Fluid Calculation (Holliday-Segar) ---
# Calculate standard maintenance first
if weight_kg <= 10:
    maintenance_standard_ml = weight_kg * 100
    maintenance_calc_str = f"({weight_kg} kg * 100 mL/kg)"
else: # This branch applies to a 12 kg child
    maintenance_standard_ml = (10 * 100) + ((weight_kg - 10) * 50)
    maintenance_calc_str = f"(10 kg * 100 mL/kg) + (({weight_kg} - 10) kg * 50 mL/kg)"

# Adjust for mechanical ventilation
maintenance_adjusted_ml = maintenance_standard_ml * ventilation_adjustment_factor

print("--- 2. Daily Maintenance Fluid Volume ---")
print("First, calculate standard maintenance fluid using the Holliday-Segar method:")
print(f"Formula: {maintenance_calc_str}")
print(f"Result: {int(maintenance_standard_ml)} mL/day")
print("\nNext, adjust for mechanical ventilation (apply 2/3 restriction):")
print(f"Formula: {int(maintenance_standard_ml)} mL * (2/3)")
# Round to two decimal places for the final answer
print(f"Result: {round(maintenance_adjusted_ml, 2)} mL/day")
print("-" * 40)


# --- 3. Deficit Replacement Volume Calculation ---
# Deficit volume (mL) = weight (kg) * percent dehydration * 10
deficit_volume_ml = weight_kg * dehydration_percent * 10

print("--- 3. Total Deficit Replacement Volume ---")
print("The calculation is based on 10% dehydration of body weight (1 kg weight loss ≈ 1000 mL fluid loss).")
print(f"Formula: {weight_kg} kg * {dehydration_percent}% * 10")
print(f"Result: {int(deficit_volume_ml)} mL")
print("-" * 40)


# --- Final Answer ---
# Collect the three final numbers
final_resuscitation = int(resuscitation_volume_ml)
final_maintenance = round(maintenance_adjusted_ml, 2)
final_deficit = int(deficit_volume_ml)

print(f"\nThe three calculated values (Resuscitation, Maintenance, Deficit) are:")
print(f"{final_resuscitation},{final_maintenance},{final_deficit}")