import math

def calculate_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on a three-phase regimen.
    """
    # Patient data
    weight_kg = 12
    resuscitation_dose_ml_kg = 30
    dehydration_percent = 0.10  # 10%

    print(f"Patient Weight: {weight_kg} kg\n")

    # --- 1. Resuscitation Volume Calculation ---
    resuscitation_volume = resuscitation_dose_ml_kg * weight_kg
    print("1. Initial Resuscitation Volume (Phase 1):")
    print(f"   Calculation: {resuscitation_dose_ml_kg} mL/kg * {weight_kg} kg = {int(resuscitation_volume)} mL")
    print("-" * 40)

    # --- 2. Daily Maintenance Volume Calculation (Holliday-Segar Method) ---
    maintenance_first_10kg = 10 * 100
    maintenance_next_kg = (weight_kg - 10) * 50
    maintenance_volume = maintenance_first_10kg + maintenance_next_kg
    
    print("2. Daily Maintenance Volume (Phase 2):")
    print(f"   For first 10 kg: 10 kg * 100 mL/kg = {int(maintenance_first_10kg)} mL")
    print(f"   For remaining {weight_kg - 10} kg: {weight_kg - 10} kg * 50 mL/kg = {int(maintenance_next_kg)} mL")
    print(f"   Total Calculation: {int(maintenance_first_10kg)} mL + {int(maintenance_next_kg)} mL = {int(maintenance_volume)} mL/day")
    print("-" * 40)

    # --- 3. Deficit Replacement Volume Calculation ---
    # Deficit in kg = weight (kg) * % dehydration. 1kg deficit ~ 1L fluid.
    deficit_volume = weight_kg * dehydration_percent * 1000
    print("3. Total Deficit Replacement Volume (Phase 3):")
    print(f"   Calculation: {weight_kg} kg * {dehydration_percent:.2f} (10%) = {weight_kg * dehydration_percent:.1f} L")
    print(f"   In milliliters: {weight_kg * dehydration_percent:.1f} L * 1000 mL/L = {int(deficit_volume)} mL")
    print("-" * 40)
    
    # --- Final Answer Format ---
    print("Final result (Resuscitation, Maintenance, Deficit):")
    print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")

calculate_fluids()