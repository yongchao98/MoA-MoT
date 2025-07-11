def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a given scenario.
    """
    weight_kg = 12
    
    # --- 1. Initial Resuscitation Volume ---
    # Formula: 30 mL/kg
    resuscitation_bolus_rate = 30
    resuscitation_volume = weight_kg * resuscitation_bolus_rate
    
    # --- 2. Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
    # Formula:
    # - 100 mL/kg for the first 10 kg
    # - 50 mL/kg for the next 10 kg (11-20 kg)
    # - 20 mL/kg for weight > 20 kg
    maintenance_volume = 0
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
        
    # --- 3. Total Deficit Replacement Fluid Volume ---
    # Formula: 10% of body weight in kg, converted to mL (1 kg water loss = 1000 mL)
    dehydration_percentage = 0.10
    deficit_kg = weight_kg * dehydration_percentage
    deficit_volume = deficit_kg * 1000  # Convert kg to mL
    
    # Print the final calculated values in the required format
    # The values are cast to integers to remove any decimal points.
    print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")

calculate_pediatric_fluids()
<<<360,1100,1200>>>