import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen
    for a 12 kg child with severe dehydration.
    """
    # Patient Data
    weight_kg = 12.0
    resuscitation_rate_ml_per_kg = 30.0
    dehydration_percentage = 10.0

    print("--- Calculating Fluid Volumes for a 12 kg Child ---")

    # 1. Phase 1: Initial Resuscitation Bolus
    # Formula: Rate (mL/kg) * Weight (kg)
    resuscitation_volume = resuscitation_rate_ml_per_kg * weight_kg
    print("\n1. Resuscitation Volume Calculation:")
    print(f"   Equation: {resuscitation_rate_ml_per_kg} mL/kg * {weight_kg} kg")
    print(f"   Result: {resuscitation_volume} mL")

    # 2. Phase 2: Daily Maintenance Fluids (Holliday-Segar Method)
    # Formula for 10-20 kg: 1000 mL (for first 10kg) + 50 mL/kg for each kg over 10
    maintenance_first_10kg = 100.0 * 10.0
    maintenance_over_10kg = 50.0 * (weight_kg - 10.0)
    daily_maintenance_volume = maintenance_first_10kg + maintenance_over_10kg
    print("\n2. Daily Maintenance Volume Calculation:")
    print(f"   Equation: (100 mL/kg * 10 kg) + (50 mL/kg * {weight_kg - 10.0} kg)")
    print(f"   Result: {maintenance_first_10kg} mL + {maintenance_over_10kg} mL = {daily_maintenance_volume} mL")

    # 3. Phase 3: Deficit Replacement Fluids
    # Formula for % dehydration: 10 mL/kg for each 1% of dehydration
    deficit_rate_ml_per_kg = 10.0 * dehydration_percentage
    deficit_volume = deficit_rate_ml_per_kg * weight_kg
    print("\n3. Deficit Replacement Volume Calculation:")
    print(f"   Equation: ({dehydration_percentage}% dehydration * 10 mL/kg) * {weight_kg} kg")
    print(f"   Result: {deficit_rate_ml_per_kg} mL/kg * {weight_kg} kg = {deficit_volume} mL")

    # Final formatted output
    print("\n--- Final Answer (Resuscitation, Maintenance, Deficit) ---")
    # Using math.floor to ensure integer output as fractional mL are not practical.
    print(f"{math.floor(resuscitation_volume)},{math.floor(daily_maintenance_volume)},{math.floor(deficit_volume)}")

# Run the calculation
calculate_pediatric_fluids()