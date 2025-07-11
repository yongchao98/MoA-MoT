def calculate_pediatric_fluids():
    """
    Calculates fluid volumes for a pediatric patient based on a three-phase regimen.
    """
    weight_kg = 12

    # --- Phase 1: Initial Resuscitation Bolus ---
    # Formula: 30 mL/kg
    resuscitation_volume = 30 * weight_kg

    # --- Phase 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
    # Formula: 100 mL/kg for first 10 kg, 50 mL/kg for next 10 kg.
    if weight_kg <= 10:
        maintenance_volume = 100 * weight_kg
    else:
        # 1000 mL for the first 10 kg
        # plus 50 mL/kg for the remaining weight
        maintenance_volume = (100 * 10) + (50 * (weight_kg - 10))

    # --- Phase 3: Total Deficit Replacement ---
    # Formula: 10% of body weight in kg, converted to mL (1 kg water ~ 1000 mL)
    dehydration_percent = 0.10
    deficit_volume = weight_kg * dehydration_percent * 1000

    # The problem requests the final numbers separated by a comma.
    # We will ensure the output is in integer format for clarity.
    final_output = f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}"
    
    # Print the final result in the required format.
    print(final_output)

calculate_pediatric_fluids()