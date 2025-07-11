# Patient and prescription data
patient_weight_kg = 12
bolus_rate_ml_per_kg = 30
dehydration_percent = 10

# --- 1. Calculate the total volume for the initial resuscitation phase ---
# Equation: Volume = Rate (mL/kg) * Weight (kg)
resuscitation_volume_ml = bolus_rate_ml_per_kg * patient_weight_kg

# --- 2. Calculate the daily maintenance fluid volume for 24 hours ---
# Using the Holliday-Segar method:
# - 100 mL/kg for the first 10 kg
# - 50 mL/kg for the next 10 kg (11-20 kg)
# - 20 mL/kg for any weight above 20 kg
maintenance_volume_ml = 0
weight_for_calc = patient_weight_kg

# Calculate for the first 10 kg
first_10kg_fluid = 10 * 100
# Calculate for the portion of weight between 11kg and 12kg
next_kg_fluid = (weight_for_calc - 10) * 50

# Equation: Total Maintenance = (10kg * 100mL/kg) + ((Weight - 10kg) * 50mL/kg)
maintenance_volume_ml = first_10kg_fluid + next_kg_fluid

# --- 3. Calculate the total deficit replacement fluid volume ---
# Equation: Deficit Volume = Weight (kg) * Dehydration % * 1000 mL/kg
# This converts the weight loss (e.g., 10% of 12kg = 1.2kg) to fluid volume in mL (1.2kg -> 1200 mL)
deficit_as_decimal = dehydration_percent / 100
deficit_volume_ml = patient_weight_kg * deficit_as_decimal * 1000

# --- Final Output ---
# Print the three calculated values, separated by commas.
# The numbers are cast to integers to represent whole milliliters.
print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_ml)},{int(deficit_volume_ml)}")