import math

def calculate_pediatric_fluids():
    """
    Calculates fluid volumes for a pediatric patient based on a three-phase regimen.
    """
    # Patient Information
    weight_kg = 12

    # --- Phase 1: Initial Resuscitation Bolus ---
    # Dosage: 30 mL/kg
    resuscitation_dosage_ml_per_kg = 30
    resuscitation_volume_ml = weight_kg * resuscitation_dosage_ml_per_kg

    # --- Phase 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
    # 100 mL/kg for the first 10 kg
    # 50 mL/kg for 11-20 kg
    # 20 mL/kg for weight > 20 kg
    maintenance_volume_24hr_ml = 0
    if weight_kg <= 10:
        maintenance_volume_24hr_ml = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_24hr_ml = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_24hr_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # --- Phase 3: Deficit Replacement ---
    # Deficit estimated at 10% of body weight
    # 1 kg weight loss is equivalent to 1000 mL fluid deficit
    deficit_percentage = 0.10
    deficit_volume_ml = weight_kg * deficit_percentage * 1000

    # --- Print the final results ---
    # The format required is three exact numbers separated by ",".
    # We use math.ceil to ensure we get whole numbers if calculations result in decimals.
    print(f"{math.ceil(resuscitation_volume_ml)},{math.ceil(maintenance_volume_24hr_ml)},{math.ceil(deficit_volume_ml)}")

# Execute the calculation and print the results
calculate_pediatric_fluids()