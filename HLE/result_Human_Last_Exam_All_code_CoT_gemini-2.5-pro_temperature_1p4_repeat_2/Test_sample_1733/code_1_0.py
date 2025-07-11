def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12
    resuscitation_rate_ml_kg = 30
    dehydration_percent = 0.10  # 10%

    # --- Calculation 1: Initial Resuscitation Volume ---
    # Formula: Rate (mL/kg) * Weight (kg)
    resuscitation_volume_ml = resuscitation_rate_ml_kg * weight_kg

    # --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
    # Formula:
    # 100 mL/kg for the first 10 kg
    # 50 mL/kg for the next 10 kg (10.1-20 kg)
    # 20 mL/kg for weight > 20 kg
    if weight_kg <= 10:
        maintenance_volume_ml_day = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_ml_day = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_ml_day = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # --- Calculation 3: Total Deficit Replacement Volume ---
    # Formula: Weight (kg) * % dehydration * 1000 mL/kg
    deficit_volume_ml = weight_kg * dehydration_percent * 1000

    # Print the results in the required format
    print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_ml_day)},{int(deficit_volume_ml)}")

calculate_pediatric_fluids()