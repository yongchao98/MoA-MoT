import math

# Patient data
weight_kg = 12
resuscitation_rate_ml_kg = 30
dehydration_percent = 0.10  # 10%

# --- Step 1: Calculate Initial Resuscitation Volume ---
resuscitation_volume = weight_kg * resuscitation_rate_ml_kg
print(f"1. Initial Resuscitation Volume Calculation:")
print(f"   {weight_kg} kg * {resuscitation_rate_ml_kg} mL/kg = {resuscitation_volume} mL\n")

# --- Step 2: Calculate Daily Maintenance Fluid Volume ---
# Using Holliday-Segar method
if weight_kg <= 10:
    standard_maintenance_volume = weight_kg * 100
    maintenance_calc_str = f"{weight_kg} kg * 100 mL/kg"
else: # For patients > 10kg
    first_10kg_fluid = 10 * 100
    remaining_weight = weight_kg - 10
    next_fluid_portion = remaining_weight * 50
    standard_maintenance_volume = first_10kg_fluid + next_fluid_portion
    maintenance_calc_str = f"(10 kg * 100 mL/kg) + ({remaining_weight} kg * 50 mL/kg)"

print(f"2. Daily Maintenance Fluid Volume Calculation:")
print(f"   Standard Maintenance = {maintenance_calc_str} = {standard_maintenance_volume} mL/day")

# Adjust for mechanical ventilation (standard practice is ~2/3 of maintenance)
ventilation_factor = 2/3
adjusted_maintenance_volume = standard_maintenance_volume * ventilation_factor
print(f"   Adjustment for Mechanical Ventilation: {standard_maintenance_volume} mL * (2/3) = {round(adjusted_maintenance_volume, 2)} mL/day")
# Round to the nearest whole number for the final answer
final_maintenance_volume = math.ceil(adjusted_maintenance_volume)
print(f"   Final Adjusted Maintenance Volume (rounded up) = {final_maintenance_volume} mL/day\n")


# --- Step 3: Calculate Total Deficit Replacement Volume ---
# 1 kg of body weight deficit is ~1000 mL of fluid deficit
deficit_volume_ml = weight_kg * dehydration_percent * 1000
print(f"3. Total Deficit Replacement Volume Calculation:")
print(f"   {weight_kg} kg * {int(dehydration_percent*100)}% * 1000 mL/L = {int(deficit_volume_ml)} mL\n")


# --- Final Answer ---
# Print the final results in the specified format
print("Final Answer:")
print(f"{resuscitation_volume},{final_maintenance_volume},{int(deficit_volume_ml)}")