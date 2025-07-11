def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12

    # --- Phase 1: Initial Resuscitation ---
    # Formula: 30 mL/kg
    resuscitation_bolus_rate_ml_per_kg = 30
    resuscitation_volume_ml = weight_kg * resuscitation_bolus_rate_ml_per_kg

    # --- Phase 2: Daily Maintenance (Holliday-Segar Method) ---
    # Formula for a child > 10 kg and <= 20 kg:
    # 1000 mL (for the first 10 kg) + 50 mL/kg for each kg over 10 kg.
    if weight_kg <= 10:
        maintenance_volume_24h_ml = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_24h_ml = (10 * 100) + ((weight_kg - 10) * 50)
    else: # weight > 20
        maintenance_volume_24h_ml = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # Note: Clinical adjustments for mechanical ventilation or fever are typically made,
    # but the base calculation as per the Holliday-Segar formula is performed here.
    # Other fluid sources (antibiotics, milk) are accounted for in the total fluid intake plan,
    # but do not change the calculated maintenance *requirement*.

    # --- Phase 3: Deficit Replacement ---
    # Formula: 10% of body weight (in kg) deficit. 1 kg loss ~ 1000 mL fluid.
    dehydration_percentage = 0.10  # 10%
    deficit_volume_ml = weight_kg * dehydration_percentage * 1000

    # --- Print the results ---
    # The final output needs to be three exact numbers separated by ","
    print(f"{int(resuscitation_volume_ml)},{int(maintenance_volume_24h_ml)},{int(deficit_volume_ml)}")

calculate_pediatric_fluids()
<<<360,1100,1200>>>