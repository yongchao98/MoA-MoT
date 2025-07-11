# Patient's weight in kg
weight_kg = 12

# Percentage of dehydration
dehydration_percent = 0.10

# Adjustment factor for mechanical ventilation (25% reduction in maintenance fluids)
mechanical_ventilation_adjustment_factor = 0.75

# --- 1. Calculate Phase 1: Resuscitation Bolus ---
# The formula is 30 mL per kg of body weight.
resuscitation_bolus_volume = 30 * weight_kg

# --- 2. Calculate Phase 2: Daily Maintenance Fluids ---
# The Holliday-Segar method is used.
# For a 12 kg child: 100 mL/kg for the first 10 kg, and 50 mL/kg for the weight above 10 kg.
# Calculation for the first 10 kg
maintenance_first_10kg = 10 * 100
# Calculation for the remaining weight (12 kg - 10 kg = 2 kg)
maintenance_next_2kg = (weight_kg - 10) * 50
# Total standard maintenance volume
unadjusted_maintenance_volume = maintenance_first_10kg + maintenance_next_2kg
# Adjust for mechanical ventilation
daily_maintenance_volume = unadjusted_maintenance_volume * mechanical_ventilation_adjustment_factor


# --- 3. Calculate Phase 3: Deficit Replacement ---
# The deficit is estimated as 10% of body weight.
# We assume 1 kg of body weight loss equals 1 L (or 1000 mL) of fluid.
deficit_in_kg = weight_kg * dehydration_percent
total_deficit_volume = deficit_in_kg * 1000

# Print the final results in the required format, ensuring they are whole numbers.
print(f"{int(resuscitation_bolus_volume)},{int(daily_maintenance_volume)},{int(total_deficit_volume)}")
