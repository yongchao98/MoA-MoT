def calculate_pediatric_fluids():
    """
    Calculates and prints pediatric fluid requirements based on a three-phase regimen.
    """
    # Patient Data
    weight_kg = 12
    resuscitation_dose_ml_kg = 30
    dehydration_percent = 10
    antibiotics_ml_day = 60
    milk_ml_day = 100

    # --- Part 1: Initial Resuscitation Bolus ---
    resuscitation_volume = weight_kg * resuscitation_dose_ml_kg
    print(f"1. Resuscitation Volume Calculation:")
    print(f"   {weight_kg} kg * {resuscitation_dose_ml_kg} mL/kg = {resuscitation_volume} mL")
    print("-" * 20)

    # --- Part 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
    maintenance_first_10kg = 10 * 100
    maintenance_next_kg = (weight_kg - 10) * 50
    total_maintenance_need = maintenance_first_10kg + maintenance_next_kg

    # Adjust for other fluids (antibiotics, milk)
    adjusted_maintenance_volume = total_maintenance_need - antibiotics_ml_day - milk_ml_day

    print(f"2. Daily Maintenance Fluid Volume Calculation:")
    print(f"   Holliday-Segar Total: (10 kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg) = {total_maintenance_need} mL")
    print(f"   Adjusted IV Volume: {total_maintenance_need} mL - {antibiotics_ml_day} mL (antibiotics) - {milk_ml_day} mL (milk) = {adjusted_maintenance_volume} mL")
    print("-" * 20)

    # --- Part 3: Deficit Replacement Fluids ---
    # 1 kg weight deficit is equivalent to 1000 mL fluid deficit.
    # Deficit (kg) = weight * (% dehydration / 100)
    # Deficit (mL) = deficit (kg) * 1000
    deficit_volume = weight_kg * (dehydration_percent / 100) * 1000
    print(f"3. Deficit Replacement Volume Calculation:")
    print(f"   {weight_kg} kg * {dehydration_percent}% dehydration = {weight_kg * dehydration_percent / 100} L = {int(deficit_volume)} mL")
    print("-" * 20)

    # --- Final Formatted Output ---
    print("Final Answer (Resuscitation, Maintenance, Deficit):")
    # Using int() to ensure whole numbers as fluid volumes are typically managed this way.
    print(f"{int(resuscitation_volume)},{int(adjusted_maintenance_volume)},{int(deficit_volume)}")

calculate_pediatric_fluids()