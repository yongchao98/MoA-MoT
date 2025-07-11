# Patient's weight in kg
weight_kg = 12

# Phase 1: Initial Resuscitation
# A bolus of 30 mL/kg
bolus_dose_ml_per_kg = 30
resuscitation_volume = weight_kg * bolus_dose_ml_per_kg

# Phase 2: Daily Maintenance Fluids (Holliday-Segar Method)
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10-20 kg
# The patient's weight is 12 kg.
if weight_kg <= 10:
    maintenance_fluid_24h = weight_kg * 100
else:
    maintenance_fluid_24h = (10 * 100) + ((weight_kg - 10) * 50)

# Note: Additional fluid sources like antibiotics (60 ml/day) and enteral nutrition (100 ml/day)
# are part of the total fluid intake but do not alter the calculated maintenance requirement itself.
# Adjustments for mechanical ventilation are also considered clinically but not applied here without specific orders.

# Phase 3: Deficit Replacement
# Deficit is estimated at 10% of body weight (10% = 0.10)
# 1 kg of body weight loss is equivalent to 1000 mL of fluid deficit.
dehydration_percentage = 0.10
deficit_volume = weight_kg * dehydration_percentage * 1000

# Print the final results in the required format
print(f"{int(resuscitation_volume)},{int(maintenance_fluid_24h)},{int(deficit_volume)}")
