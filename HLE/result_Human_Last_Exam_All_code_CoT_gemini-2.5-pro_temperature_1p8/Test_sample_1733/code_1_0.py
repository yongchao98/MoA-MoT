def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12.0
    
    # --- Phase 1: Initial Resuscitation ---
    # Formula: 30 mL/kg over 1 hour
    bolus_per_kg = 30.0
    resuscitation_volume = weight_kg * bolus_per_kg

    # --- Phase 2: Maintenance Fluids ---
    # Using Holliday-Segar method for a 12 kg child
    # 100 mL/kg for first 10 kg, 50 mL/kg for 11-20 kg
    if weight_kg <= 10:
        maintenance_volume_24hr = weight_kg * 100.0
    else:
        maintenance_volume_24hr = (10 * 100.0) + ((weight_kg - 10) * 50.0)

    # Adjustment for mechanical ventilation (reduces insensible losses)
    # A standard clinical adjustment is a 25% reduction.
    ventilation_reduction_factor = 0.75 
    adjusted_maintenance_volume_24hr = maintenance_volume_24hr * ventilation_reduction_factor

    # --- Phase 3: Deficit Replacement ---
    # Estimated at 10% of body weight, given over 48 hours.
    # 1 kg of weight loss from dehydration ~ 1 L (1000 mL) of fluid loss.
    dehydration_percentage = 0.10
    total_deficit_volume = weight_kg * dehydration_percentage * 1000.0

    # The problem asks for the three exact numbers separated by ","
    # The note about antibiotics and milk would be used to calculate the final IV rate,
    # but not the total physiological requirements asked for here.
    
    # Print the final results as integers, as the calculations result in whole numbers.
    print(f"{int(resuscitation_volume)},{int(adjusted_maintenance_volume_24hr)},{int(total_deficit_volume)}")

calculate_pediatric_fluids()
<<<360,825,1200>>>