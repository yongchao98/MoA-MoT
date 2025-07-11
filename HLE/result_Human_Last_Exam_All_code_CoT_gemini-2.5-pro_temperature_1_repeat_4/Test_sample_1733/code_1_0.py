import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a clinical scenario.
    """
    # Patient parameters
    weight_kg = 12
    dehydration_percentage = 0.10  # 10%
    
    # --- Calculation 1: Initial Resuscitation Bolus ---
    # Formula: 30 mL/kg
    bolus_per_kg = 30
    resuscitation_volume = bolus_per_kg * weight_kg

    # --- Calculation 2: Daily Maintenance Fluid (Holliday-Segar Method) ---
    maintenance_volume_unadjusted = 0
    if weight_kg <= 10:
        maintenance_volume_unadjusted = weight_kg * 100
    else:
        maintenance_volume_unadjusted = 10 * 100 + (weight_kg - 10) * 50
        
    # Adjust for mechanical ventilation (typically 70-80% of calculated maintenance)
    # We will use a factor of 0.75
    ventilation_adjustment_factor = 0.75
    maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_adjustment_factor

    # --- Calculation 3: Total Deficit Replacement ---
    # Formula: 10% of body weight in kg, converted to mL (1 kg ~ 1000 mL)
    deficit_volume = weight_kg * dehydration_percentage * 1000

    # --- Final Output ---
    # Print the results as three integer numbers separated by commas.
    print(f"{int(resuscitation_volume)},{int(maintenance_volume_adjusted)},{int(deficit_volume)}")

calculate_pediatric_fluids()
<<<360,825,1200>>>