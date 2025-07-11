import math

# --- Patient and Prescription Data ---
weight_kg = 12
resuscitation_bolus_rate_ml_per_kg = 30
dehydration_percentage = 0.10 # 10%
# Note: Mechanical ventilation reduces insensible losses. A common adjustment is a 25% reduction in maintenance fluids.
ventilation_reduction_factor = 0.75 # Corresponds to a 25% reduction

# --- Calculation Step 1: Initial Resuscitation Bolus ---
resuscitation_volume = resuscitation_bolus_rate_ml_per_kg * weight_kg
print(f"1. Initial Resuscitation Volume Calculation:")
print(f"   {resuscitation_bolus_rate_ml_per_kg} mL/kg * {weight_kg} kg = {resuscitation_volume} mL")
print("-" * 30)


# --- Calculation Step 2: Daily Maintenance Fluids ---
# Using the Holliday-Segar Method:
# - 100 mL/kg for the first 10 kg
# - 50 mL/kg for the next 10 kg (11-20 kg)
# - 20 mL/kg for the weight over 20 kg
if weight_kg <= 10:
    base_maintenance_volume = weight_kg * 100
    base_calc_str = f"{weight_kg} kg * 100 mL/kg"
elif weight_kg <= 20:
    base_maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    base_calc_str = f"(10 kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg)"
else:
    base_maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    base_calc_str = f"(10 kg * 100 mL/kg) + (10 kg * 50 mL/kg) + ({weight_kg - 20} kg * 20 mL/kg)"

# Adjust for mechanical ventilation
adjusted_maintenance_volume = base_maintenance_volume * ventilation_reduction_factor

print("2. Daily Maintenance Fluid Volume Calculation:")
print(f"   Base volume (Holliday-Segar): {base_calc_str} = {base_maintenance_volume} mL/day")
print(f"   Adjustment for Mechanical Ventilation: {base_maintenance_volume} mL * {ventilation_reduction_factor} = {adjusted_maintenance_volume} mL/day")
print("-" * 30)


# --- Calculation Step 3: Deficit Replacement Volume ---
# Deficit is 10% of body weight. 1 kg of weight loss = 1 L of fluid deficit.
deficit_volume_l = weight_kg * dehydration_percentage
deficit_volume_ml = deficit_volume_l * 1000 # Convert Liters to Milliliters

print("3. Total Deficit Replacement Volume Calculation:")
print(f"   {weight_kg} kg * {dehydration_percentage * 100}% dehydration = {deficit_volume_l} L = {deficit_volume_ml} mL")
print("-" * 30)

# --- Final Answer ---
# Presenting the three required values, rounded to the nearest whole number.
final_resuscitation = int(round(resuscitation_volume))
final_maintenance = int(round(adjusted_maintenance_volume))
final_deficit = int(round(deficit_volume_ml))

print(f"\nFinal Calculated Volumes (Resuscitation, Maintenance, Deficit):")
print(f"{final_resuscitation},{final_maintenance},{final_deficit}")

# The user is instructed not to copy the result, so this final line is for programmatic parsing.
# The answer is presented in the format <<<resuscitation,maintenance,deficit>>>
# <<<360,825,1200>>>