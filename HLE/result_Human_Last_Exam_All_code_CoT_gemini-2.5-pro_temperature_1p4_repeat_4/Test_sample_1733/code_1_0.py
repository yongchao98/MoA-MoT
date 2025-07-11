import math

# --- Patient Data ---
weight_kg = 12
dehydration_percent = 0.10  # 10%
iv_antibiotics_ml_day = 60
enteral_milk_ml_day = 100
# Note: Patient is mechanically ventilated, but no reduction factor is specified.
# Standard calculations will be performed.
# Note: IV antibiotics and enteral milk are other sources of intake but do not change the
# calculated resuscitation, maintenance, or deficit volumes themselves.

# --- 1. Phase 1: Initial Resuscitation Volume ---
print("Part 1: Initial Resuscitation Volume Calculation")
resuscitation_rate_ml_per_kg = 30
resuscitation_volume_ml = weight_kg * resuscitation_rate_ml_per_kg
print(f"The calculation is: {weight_kg} kg * {resuscitation_rate_ml_per_kg} mL/kg = {int(resuscitation_volume_ml)} mL")
print("-" * 20)

# --- 2. Phase 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
print("Part 2: Daily Maintenance Fluid Volume Calculation")
maintenance_fluid_24hr_ml = 0
if weight_kg <= 10:
    maintenance_fluid_24hr_ml = weight_kg * 100
    print(f"The calculation is: {weight_kg} kg * 100 mL/kg = {int(maintenance_fluid_24hr_ml)} mL")
elif weight_kg <= 20:
    first_10kg_fluid = 10 * 100
    remaining_weight = weight_kg - 10
    next_kg_fluid = remaining_weight * 50
    maintenance_fluid_24hr_ml = first_10kg_fluid + next_kg_fluid
    print(f"For the first 10 kg: 10 kg * 100 mL/kg = {int(first_10kg_fluid)} mL")
    print(f"For the remaining {remaining_weight} kg: {remaining_weight} kg * 50 mL/kg = {int(next_kg_fluid)} mL")
    print(f"The total calculation is: {int(first_10kg_fluid)} mL + {int(next_kg_fluid)} mL = {int(maintenance_fluid_24hr_ml)} mL")
else: # For weight > 20 kg
    first_10kg_fluid = 10 * 100
    next_10kg_fluid = 10 * 50
    remaining_weight = weight_kg - 20
    final_kg_fluid = remaining_weight * 20
    maintenance_fluid_24hr_ml = first_10kg_fluid + next_10kg_fluid + final_kg_fluid
    print(f"For the first 10 kg: 10 kg * 100 mL/kg = {int(first_10kg_fluid)} mL")
    print(f"For the next 10 kg: 10 kg * 50 mL/kg = {int(next_10kg_fluid)} mL")
    print(f"For the remaining {remaining_weight} kg: {remaining_weight} kg * 20 mL/kg = {int(final_kg_fluid)} mL")
    print(f"The total calculation is: {int(first_10kg_fluid)} mL + {int(next_10kg_fluid)} mL + {int(final_kg_fluid)} mL = {int(maintenance_fluid_24hr_ml)} mL")
print("-" * 20)


# --- 3. Phase 3: Total Deficit Replacement Fluid Volume ---
print("Part 3: Total Deficit Replacement Volume Calculation")
# 1 kg of weight loss equals 1 L (1000 mL) of fluid deficit.
deficit_volume_ml = weight_kg * dehydration_percent * 1000
print(f"The calculation is: {weight_kg} kg * {dehydration_percent:.2f} * 1000 mL/L = {int(deficit_volume_ml)} mL")
print("-" * 20)


# --- Final Answer ---
print("Final calculated volumes (Resuscitation, Daily Maintenance, Total Deficit):")
final_answer = f"{int(resuscitation_volume_ml)},{int(maintenance_fluid_24hr_ml)},{int(deficit_volume_ml)}"
print(final_answer)

# Do not remove the '<<<...>>>' line below, it is used for answer parsing.
print(f'<<<{final_answer}>>>')