import math

def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12
    resuscitation_rate_ml_kg = 30
    deficit_percentage = 0.10  # 10% dehydration
    ventilation_reduction_factor = 0.75 # 25% reduction for mechanical ventilation

    # --- Calculation 1: Initial Resuscitation Volume ---
    # Formula: weight * rate
    resuscitation_volume = weight_kg * resuscitation_rate_ml_kg

    # --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
    # 100 mL/kg for the first 10 kg
    # 50 mL/kg for the next 10 kg (11-20 kg)
    # 20 mL/kg for weight over 20 kg
    if weight_kg <= 10:
        maintenance_volume_unadjusted = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_unadjusted = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_unadjusted = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # Adjust for mechanical ventilation (25% reduction)
    maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_reduction_factor

    # --- Calculation 3: Total Deficit Replacement Volume ---
    # Formula: weight * deficit_percentage * 1000 (to convert kg to mL)
    deficit_volume = weight_kg * deficit_percentage * 1000

    # --- Print the final results ---
    # The problem asks for the three calculated values.
    # Note: In a clinical setting, the daily IV fluid rate would be derived from the maintenance and deficit needs,
    # and then other intake like antibiotics and milk would be subtracted from that total to get the final IV rate.
    # However, the question asks for the calculated resuscitation, maintenance, and deficit volumes themselves.
    
    # Using math.ceil to round up to the nearest whole number for practical clinical application
    final_resuscitation = math.ceil(resuscitation_volume)
    final_maintenance = math.ceil(maintenance_volume_adjusted)
    final_deficit = math.ceil(deficit_volume)

    print(f"{final_resuscitation},{final_maintenance},{final_deficit}")

calculate_pediatric_fluids()
<<<360,825,1200>>>