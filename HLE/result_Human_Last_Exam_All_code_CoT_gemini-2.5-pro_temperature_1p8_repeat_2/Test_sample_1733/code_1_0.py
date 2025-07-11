import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen.
    """
    
    # Patient Information
    weight_kg = 12

    # --- Calculation 1: Initial Resuscitation Bolus (Phase 1) ---
    # Formula: 30 mL/kg
    resuscitation_rate_ml_per_kg = 30
    resuscitation_volume = weight_kg * resuscitation_rate_ml_per_kg

    # --- Calculation 2: Daily Maintenance Fluids (Phase 2 - Holliday-Segar) ---
    # Formula for 24 hours:
    # - 100 mL/kg for the first 10 kg
    # - 50 mL/kg for the next 10 kg (10.1-20 kg)
    # - 20 mL/kg for weight over 20 kg
    maintenance_volume = 0
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    # --- Calculation 3: Deficit Replacement Fluids (Phase 3) ---
    # Formula: 10% of body weight in kg, where 1 kg ~ 1000 mL of fluid.
    dehydration_percentage = 0.10 # 10%
    deficit_volume = weight_kg * dehydration_percentage * 1000 # Convert kg deficit to mL

    # --- Final Output ---
    # The final answer requires three exact numbers separated by ","
    # math.trunc is used to ensure we have whole numbers as volumes are typically administered this way.
    print(f"{math.trunc(resuscitation_volume)},{math.trunc(maintenance_volume)},{math.trunc(deficit_volume)}")

# Execute the function to get the result
calculate_pediatric_fluids()