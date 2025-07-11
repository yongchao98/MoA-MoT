import math

# Patient Information
weight_kg = 12

# Phase 1: Initial Resuscitation Bolus
# Formula: 30 mL/kg
bolus_rate_ml_per_kg = 30
resuscitation_volume = weight_kg * bolus_rate_ml_per_kg

print("--- Fluid Calculations ---")
print(f"Patient Weight: {weight_kg} kg\n")

print("1. Initial Resuscitation Volume Calculation:")
print(f"   {bolus_rate_ml_per_kg} mL/kg * {weight_kg} kg = {resuscitation_volume} mL")
print("-" * 25)

# Phase 2: Daily Maintenance Fluids (Holliday-Segar Method)
# Note: A 25% reduction (multiplied by 0.75) is applied for mechanical ventilation.
maintenance_for_first_10kg = 10 * 100
remaining_weight = weight_kg - 10
maintenance_for_remaining_weight = remaining_weight * 50
total_unadjusted_maintenance = maintenance_for_first_10kg + maintenance_for_remaining_weight
ventilation_adjustment_factor = 0.75
adjusted_maintenance_volume = total_unadjusted_maintenance * ventilation_adjustment_factor

print("2. Daily Maintenance Fluid Volume Calculation (for 24 hours):")
print(f"   Unadjusted (Holliday-Segar): ({10} kg * 100 mL/kg) + ({remaining_weight} kg * 50 mL/kg) = {total_unadjusted_maintenance} mL")
print(f"   Adjusted for Ventilation: {total_unadjusted_maintenance} mL * {ventilation_adjustment_factor} = {math.ceil(adjusted_maintenance_volume)} mL")
print("-" * 25)

# Phase 3: Deficit Replacement Fluids
# Formula: 10% of body weight (1 kg weight loss = 1000 mL fluid loss)
dehydration_percentage = 0.10
deficit_volume = weight_kg * dehydration_percentage * 1000

print("3. Total Deficit Replacement Volume Calculation:")
print(f"   {weight_kg} kg * {int(dehydration_percentage * 100)}% dehydration = {int(deficit_volume)} mL")
print("-" * 25)

# Final results formatted as requested
final_resuscitation = int(resuscitation_volume)
final_maintenance = int(math.ceil(adjusted_maintenance_volume))
final_deficit = int(deficit_volume)

print("\n--- Final Answer Summary ---")
print(f"Resuscitation: {final_resuscitation} mL")
print(f"Maintenance: {final_maintenance} mL")
print(f"Deficit: {final_deficit} mL")