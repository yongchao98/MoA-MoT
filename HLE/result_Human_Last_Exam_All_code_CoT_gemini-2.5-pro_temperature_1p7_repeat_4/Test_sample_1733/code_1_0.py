# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Formula: 30 mL/kg bolus
bolus_rate_ml_per_kg = 30
resuscitation_volume_ml = weight_kg * bolus_rate_ml_per_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# For a 12 kg child:
# 100 mL/kg for the first 10 kg
# 50 mL/kg for the next 10 kg (i.e., from 10.1 kg to 20 kg)
maintenance_for_first_10kg = 10 * 100
maintenance_for_remaining_weight = (weight_kg - 10) * 50
daily_maintenance_volume_ml = maintenance_for_first_10kg + maintenance_for_remaining_weight

# --- Calculation 3: Total Deficit Replacement Volume ---
# Deficit is 10% of body weight. 1 kg of weight loss equals 1000 mL of fluid deficit.
dehydration_percentage = 0.10
deficit_in_kg = weight_kg * dehydration_percentage
deficit_replacement_volume_ml = deficit_in_kg * 1000

# --- Print the final results ---
# The problem asks for the three exact numbers separated by ","
# 1. Total volume for initial resuscitation
# 2. Daily maintenance fluid volume
# 3. Total deficit replacement volume
print(f"{int(resuscitation_volume_ml)},{int(daily_maintenance_volume_ml)},{int(deficit_replacement_volume_ml)}")