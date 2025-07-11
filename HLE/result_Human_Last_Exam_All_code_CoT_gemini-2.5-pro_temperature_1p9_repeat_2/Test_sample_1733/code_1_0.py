import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a clinical scenario.
    """
    # Patient parameters
    weight_kg = 12
    dehydration_percentage = 0.10  # 10%

    # --- Calculation 1: Initial Resuscitation Volume ---
    bolus_rate_ml_per_kg = 30
    resuscitation_volume = weight_kg * bolus_rate_ml_per_kg
    
    print("1. Initial Resuscitation Volume Calculation:")
    print(f"   {weight_kg} kg * {bolus_rate_ml_per_kg} mL/kg = {int(resuscitation_volume)} mL")
    print("-" * 30)

    # --- Calculation 2: Daily Maintenance Fluid Volume ---
    # Holliday-Segar Method
    if weight_kg <= 10:
        maintenance_base = weight_kg * 100
        maintenance_calc_str = f"({weight_kg} kg * 100 mL/kg)"
    else:
        maintenance_base = (10 * 100) + ((weight_kg - 10) * 50)
        maintenance_calc_str = f"((10 kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg))"
    
    # Adjust for mechanical ventilation (reduced insensible losses)
    # Applying a standard 20% reduction factor (0.8)
    ventilation_adjustment_factor = 0.8
    maintenance_adjusted = maintenance_base * ventilation_adjustment_factor
    
    print("2. Daily Maintenance Fluid Volume Calculation:")
    print(f"   Holliday-Segar base: {maintenance_calc_str} = {int(maintenance_base)} mL")
    print(f"   Adjustment for mechanical ventilation ({int((1-ventilation_adjustment_factor)*100)}% reduction):")
    print(f"   {int(maintenance_base)} mL * {ventilation_adjustment_factor} = {int(maintenance_adjusted)} mL")
    print("-" * 30)
    
    # --- Calculation 3: Total Deficit Replacement Volume ---
    # 1% dehydration = 10 mL/kg fluid deficit
    # 10% dehydration = 100 mL/kg
    deficit_rate_ml_per_kg = 100 
    deficit_volume = weight_kg * deficit_rate_ml_per_kg

    print("3. Total Deficit Replacement Volume Calculation:")
    # We calculate 10% of body weight in mL. 1 kg = 1000 g, and 1 g water ~ 1 mL. So 10% of 12 kg is 1.2 kg = 1200g = 1200 mL
    print(f"   10% of {weight_kg} kg body weight = {weight_kg} kg * {deficit_rate_ml_per_kg} mL/kg = {int(deficit_volume)} mL")
    print("-" * 30)

    # --- Final formatted output ---
    final_answer = f"{int(resuscitation_volume)},{int(maintenance_adjusted)},{int(deficit_volume)}"
    print("Final calculated volumes (Resuscitation, Maintenance, Deficit):")
    print(final_answer)


# Run the calculation and print the results
calculate_pediatric_fluids()