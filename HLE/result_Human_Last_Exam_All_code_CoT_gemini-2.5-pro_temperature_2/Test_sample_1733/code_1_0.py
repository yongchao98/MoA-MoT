import math

def calculate_fluid_volumes():
    """
    Calculates fluid volumes for a pediatric patient based on a three-phase regimen.
    """
    weight_kg = 12
    
    # --- Phase 1: Initial Resuscitation ---
    bolus_rate_ml_per_kg = 30
    resuscitation_volume = bolus_rate_ml_per_kg * weight_kg
    print("1. Initial Resuscitation Volume Calculation:")
    print(f"   Formula: Weight (kg) * Bolus Rate (mL/kg)")
    print(f"   Calculation: {weight_kg} kg * {bolus_rate_ml_per_kg} mL/kg = {resuscitation_volume} mL")
    print("-" * 30)

    # --- Phase 2: Daily Maintenance (Holliday-Segar Method) ---
    maintenance_volume_24h = 0
    print("2. Daily Maintenance Fluid Volume Calculation (Holliday-Segar):")
    if weight_kg <= 10:
        maintenance_volume_24h = weight_kg * 100
        print(f"   For first {weight_kg} kg: {weight_kg} kg * 100 mL/kg = {maintenance_volume_24h} mL")
    else:
        # For the first 10 kg
        vol_first_10kg = 10 * 100
        print(f"   For first 10 kg: 10 kg * 100 mL/kg = {vol_first_10kg} mL")
        
        # For the weight above 10 kg
        remaining_weight = weight_kg - 10
        vol_next_kg = remaining_weight * 50
        print(f"   For remaining {remaining_weight} kg: {remaining_weight} kg * 50 mL/kg = {vol_next_kg} mL")
        
        maintenance_volume_24h = vol_first_10kg + vol_next_kg
    
    print(f"   Total Daily Maintenance Volume: {vol_first_10kg} mL + {vol_next_kg} mL = {maintenance_volume_24h} mL")
    print("-" * 30)

    # --- Phase 3: Deficit Replacement ---
    dehydration_percent = 0.10  # 10%
    # 1 kg of weight loss ~ 1 L (1000 mL) of fluid deficit
    deficit_volume = weight_kg * dehydration_percent * 1000
    print("3. Total Deficit Replacement Volume Calculation:")
    print(f"   Formula: Weight (kg) * Dehydration % * 1000 mL/kg")
    print(f"   Calculation: {weight_kg} kg * {dehydration_percent:.2f} * 1000 mL/kg = {deficit_volume} mL")
    print("-" * 30)

    # --- Final formatted output ---
    final_answer = f"{int(resuscitation_volume)},{int(maintenance_volume_24h)},{int(deficit_volume)}"
    print(f"The three calculated volumes are: {final_answer}")
    
    # Adhering to the required final output format
    print(f"<<<{final_answer}>>>")


calculate_fluid_volumes()