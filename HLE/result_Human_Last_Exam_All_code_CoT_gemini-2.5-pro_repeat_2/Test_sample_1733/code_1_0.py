import math

# Patient's weight in kg
weight_kg = 12

# --- Calculation 1: Initial Resuscitation Volume ---
# Formula: 30 mL/kg bolus
resuscitation_rate_ml_per_kg = 30
resuscitation_volume_ml = weight_kg * resuscitation_rate_ml_per_kg

# --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
# For the first 10 kg: 100 mL/kg/day
# For 11-20 kg: 50 mL/kg/day
# For >20 kg: 20 mL/kg/day
maintenance_volume_ml = 0
if weight_kg <= 10:
    maintenance_volume_ml = weight_kg * 100
elif weight_kg <= 20:
    maintenance_volume_ml = (10 * 100) + ((weight_kg - 10) * 50)
else:
    maintenance_volume_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

# Note: In a clinical setting, this calculated maintenance requirement would be adjusted.
# For example, it might be reduced due to mechanical ventilation, and the volume of
# IV antibiotics and enteral milk would be subtracted from the total to determine
# the final IV fluid rate. The question, however, asks for the calculated daily
# maintenance volume itself.

# --- Calculation 3: Total Deficit Replacement Fluid Volume ---
# Formula: 10% of body weight (in kg), where 1 kg weight loss = 1000 mL fluid deficit
dehydration_percentage = 0.10
deficit_volume_ml = weight_kg * dehydration_percentage * 1000 # convert kg deficit to mL

# --- Final Output ---
# The problem asks for three exact numbers separated by a comma.
# 1. Total volume for initial resuscitation.
# 2. Daily maintenance fluid volume for 24 hours.
# 3. Total deficit replacement fluid volume over 48 hours.
print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_ml)},{int(deficit_volume_ml)}")