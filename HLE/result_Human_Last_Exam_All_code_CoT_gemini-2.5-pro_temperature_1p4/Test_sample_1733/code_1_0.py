import math

# Step 1: Define Patient Data
weight_kg = 12
dehydration_percent = 0.10  # 10%

# Step 2: Calculate Resuscitation Volume (Phase 1)
# Formula: 30 mL/kg
resuscitation_volume = 30 * weight_kg

# Step 3: Calculate Maintenance Volume (Phase 2)
# Holliday-Segar method for a 12 kg child:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the weight above 10 kg
maintenance_base = (10 * 100) + ((weight_kg - 10) * 50)

# Adjust for mechanical ventilation (standard 25% reduction)
# Humidified air from the ventilator reduces insensible losses.
maintenance_adjusted = maintenance_base * 0.75

# Step 4: Calculate Deficit Replacement Volume (Phase 3)
# Formula: Dehydration % * Weight (kg) * 1000 mL
deficit_volume = dehydration_percent * weight_kg * 1000

# Step 5: Format and Print Output
# The results are rounded to the nearest whole number for clinical clarity.
print(f"{math.ceil(resuscitation_volume)},{math.ceil(maintenance_adjusted)},{math.ceil(deficit_volume)}")