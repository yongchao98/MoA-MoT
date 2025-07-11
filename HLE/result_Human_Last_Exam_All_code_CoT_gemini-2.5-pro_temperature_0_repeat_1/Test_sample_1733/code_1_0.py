def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    weight_kg = 12
    resuscitation_rate_ml_per_kg = 30
    deficit_percentage = 0.10  # 10%

    # 1. Calculate the total volume for the initial resuscitation phase.
    # Formula: Volume = Rate (mL/kg) * Weight (kg)
    resuscitation_volume = resuscitation_rate_ml_per_kg * weight_kg

    # 2. Calculate the daily maintenance fluid volume using the Holliday-Segar method.
    # For the first 10 kg: 100 mL/kg/day
    # For 11-20 kg: 50 mL/kg/day
    # For >20 kg: 20 mL/kg/day
    maintenance_volume = 0
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # 3. Calculate the total deficit replacement fluid volume.
    # Deficit is 10% of body weight. 1 kg weight loss = 1000 mL fluid deficit.
    deficit_weight_kg = weight_kg * deficit_percentage
    deficit_volume = deficit_weight_kg * 1000

    # Print the final results in the specified format.
    # The results are cast to integers as fluid volumes are typically managed in whole numbers.
    print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")

calculate_pediatric_fluids()
<<<360,1100,1200>>>