def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12
    resuscitation_rate_ml_per_kg = 30
    dehydration_percent = 0.10  # 10%

    # --- Calculation 1: Initial Resuscitation Volume ---
    resuscitation_volume = weight_kg * resuscitation_rate_ml_per_kg

    # --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
    maintenance_volume_unadjusted = 0
    if weight_kg <= 10:
        maintenance_volume_unadjusted = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_unadjusted = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_unadjusted = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # Adjust for mechanical ventilation (standard 25% reduction)
    ventilation_reduction_factor = 0.75
    maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_reduction_factor

    # --- Calculation 3: Total Deficit Replacement Volume ---
    # 1 kg of body weight loss from dehydration = 1000 mL of fluid deficit
    deficit_volume_ml = weight_kg * dehydration_percent * 1000

    # --- Print the final results ---
    # The results are cast to integers as the calculations result in whole numbers.
    print(f"{int(resuscitation_volume)},{int(maintenance_volume_adjusted)},{int(deficit_volume_ml)}")

calculate_pediatric_fluids()