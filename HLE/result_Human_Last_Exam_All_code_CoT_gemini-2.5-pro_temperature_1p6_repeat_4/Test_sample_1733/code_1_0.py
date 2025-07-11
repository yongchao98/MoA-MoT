import math

# --- Patient and Prescription Data ---
weight_kg = 12
bolus_rate_ml_per_kg = 30
dehydration_percent = 0.10  # 10%

# --- 1. Calculate Initial Resuscitation Volume ---
print("1. Initial Resuscitation Volume Calculation:")
resuscitation_volume = weight_kg * bolus_rate_ml_per_kg
print(f"   The resuscitation volume is calculated by multiplying the patient's weight by the prescribed bolus rate.")
print(f"   Equation: {weight_kg} kg * {bolus_rate_ml_per_kg} mL/kg = {int(resuscitation_volume)} mL")
print("-" * 50)


# --- 2. Calculate Daily Maintenance Fluid Volume ---
print("2. Daily Maintenance Fluid Volume Calculation (Holliday-Segar method):")
# Holliday-Segar Formula:
# - 100 mL/kg/day for the first 10 kg
# - 50 mL/kg/day for the next 10 kg (from 11-20 kg)
# - 20 mL/kg/day for weight over 20 kg
if weight_kg <= 10:
    maintenance_volume_24hr = weight_kg * 100
    print(f"   For a child <= 10 kg, the rate is 100 mL/kg.")
    print(f"   Equation: {weight_kg} kg * 100 mL/kg = {int(maintenance_volume_24hr)} mL/day")
else: # This branch applies to the 12 kg child
    maintenance_for_first_10kg = 10 * 100
    remaining_weight = weight_kg - 10
    maintenance_for_remaining_weight = remaining_weight * 50
    maintenance_volume_24hr = maintenance_for_first_10kg + maintenance_for_remaining_weight
    print("   For a child > 10 kg, the calculation is split:")
    print(f"   Part 1 (first 10 kg): 10 kg * 100 mL/kg = {int(maintenance_for_first_10kg)} mL")
    print(f"   Part 2 (remaining {remaining_weight} kg): {remaining_weight} kg * 50 mL/kg = {int(maintenance_for_remaining_weight)} mL")
    print(f"   Total Daily Maintenance: {int(maintenance_for_first_10kg)} mL + {int(maintenance_for_remaining_weight)} mL = {int(maintenance_volume_24hr)} mL/day")
print("-" * 50)


# --- 3. Calculate Total Deficit Replacement Volume ---
print("3. Total Deficit Replacement Volume Calculation:")
# Deficit is estimated at 10% of body weight.
# 1 kg of acute weight loss is equivalent to 1000 mL of fluid deficit.
fluid_deficit_in_kg = weight_kg * dehydration_percent
deficit_replacement_volume = fluid_deficit_in_kg * 1000
print(f"   The fluid deficit in kg is the percentage of dehydration multiplied by body weight.")
print(f"   Equation (kg): {weight_kg} kg * {dehydration_percent:.0%} = {fluid_deficit_in_kg:.1f} kg")
print(f"   This is converted to mL (1 kg = 1000 mL).")
print(f"   Equation (mL): {fluid_deficit_in_kg:.1f} kg * 1000 mL/kg = {int(deficit_replacement_volume)} mL")
print("-" * 50)


# --- Final Formatted Answer ---
print("\nFinal Answer (Resuscitation, Maintenance, Deficit):")
print(f"{int(resuscitation_volume)},{int(maintenance_volume_24hr)},{int(deficit_replacement_volume)}")

<<<360,1100,1200>>>