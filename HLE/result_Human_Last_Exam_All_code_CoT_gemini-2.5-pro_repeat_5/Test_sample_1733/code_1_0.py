def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient Data
    weight_kg = 12
    dehydration_percentage = 0.10  # 10%

    # --- 1. Initial Resuscitation Bolus Calculation ---
    # Formula: 30 mL/kg * body weight
    resuscitation_volume = 30 * weight_kg

    # --- 2. Daily Maintenance Fluid Calculation (Holliday-Segar Method) ---
    # Formula:
    # - 100 mL/kg for the first 10 kg
    # - 50 mL/kg for the next 10 kg (10.1-20 kg)
    # - 20 mL/kg for weight over 20 kg
    # Note: This calculates the gross daily maintenance requirement. In a clinical setting,
    # this amount might be adjusted for factors like mechanical ventilation or reduced
    # to account for other fluid sources like IV antibiotics and enteral feeds.
    if weight_kg <= 10:
        maintenance_volume_24hr = 100 * weight_kg
    elif weight_kg <= 20:
        maintenance_volume_24hr = (100 * 10) + (50 * (weight_kg - 10))
    else:
        maintenance_volume_24hr = (100 * 10) + (50 * 10) + (20 * (weight_kg - 20))

    # --- 3. Total Deficit Replacement Fluid Calculation ---
    # Formula: Dehydration % * body weight (kg) * 1000 mL/L
    # 1 kg of body weight is equivalent to 1000 mL of fluid deficit.
    total_deficit_volume = dehydration_percentage * weight_kg * 1000

    # --- Print the final formatted answer ---
    # The final output is three exact numbers separated by ","
    print(f"{int(resuscitation_volume)},{int(maintenance_volume_24hr)},{int(total_deficit_volume)}")

# Execute the function
calculate_pediatric_fluids()