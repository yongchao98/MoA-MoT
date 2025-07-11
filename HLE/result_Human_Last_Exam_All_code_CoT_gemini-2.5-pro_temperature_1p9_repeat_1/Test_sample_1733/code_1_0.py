def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase management plan.
    """
    # Patient's weight in kg
    weight_kg = 12

    # --- Task 1: Calculate Initial Resuscitation Volume ---
    # Formula: 30 mL/kg
    resuscitation_volume = weight_kg * 30

    # --- Task 2: Calculate Daily Maintenance Fluid Volume (Holliday-Segar method) ---
    # For the first 10 kg: 100 mL/kg
    # For 11-20 kg: 1000 mL + 50 mL/kg for each kg over 10
    # For >20 kg: 1500 mL + 20 mL/kg for each kg over 20
    if weight_kg <= 10:
        maintenance_volume_24hr = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_24hr = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_24hr = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # --- Task 3: Calculate Total Deficit Replacement Fluid Volume ---
    # Estimated deficit is 10% of body weight.
    # 1 kg of body weight loss corresponds to ~1000 mL of fluid deficit.
    deficit_percentage = 0.10
    deficit_volume = weight_kg * deficit_percentage * 1000

    # Print the results in the specified format
    # The values will be rounded to the nearest whole number as fluid calculations are typically done in whole mL.
    print(f"{int(resuscitation_volume)},{int(maintenance_volume_24hr)},{int(deficit_volume)}")

# Execute the calculation and print the result.
calculate_pediatric_fluids()
<<<360,1100,1200>>>