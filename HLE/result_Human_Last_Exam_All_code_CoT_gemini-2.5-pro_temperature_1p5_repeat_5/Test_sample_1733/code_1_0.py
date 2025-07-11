def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen.
    """
    # --- Patient Data ---
    weight_kg = 12
    dehydration_percent = 0.10  # 10%
    bolus_per_kg = 30
    ventilation_adjustment_factor = 0.75 # Standard adjustment for mechanical ventilation

    # --- Calculation 1: Resuscitation Bolus ---
    resuscitation_volume = weight_kg * bolus_per_kg

    print("--- 1. Initial Resuscitation Volume ---")
    print(f"The calculation is: {weight_kg} kg * {bolus_per_kg} mL/kg")
    print(f"Result: {int(resuscitation_volume)} mL\n")

    # --- Calculation 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
    print("--- 2. Daily Maintenance Fluid Volume ---")
    if weight_kg <= 10:
        maintenance_volume_unadjusted = weight_kg * 100
        print(f"Holliday-Segar Calculation: {weight_kg} kg * 100 mL/kg = {int(maintenance_volume_unadjusted)} mL")
    else:
        vol_first_10 = 10 * 100
        vol_next = (weight_kg - 10) * 50
        maintenance_volume_unadjusted = vol_first_10 + vol_next
        print("Holliday-Segar Calculation:")
        print(f"For the first 10 kg: 10 kg * 100 mL/kg = {int(vol_first_10)} mL")
        print(f"For the remaining {weight_kg - 10} kg: {weight_kg - 10} kg * 50 mL/kg = {int(vol_next)} mL")
        print(f"Total standard maintenance: {int(vol_first_10)} mL + {int(vol_next)} mL = {int(maintenance_volume_unadjusted)} mL")

    # Adjust for mechanical ventilation
    maintenance_volume_adjusted = maintenance_volume_unadjusted * ventilation_adjustment_factor
    print("\nApplying adjustment for mechanical ventilation (standard 25% reduction):")
    print(f"Calculation: {int(maintenance_volume_unadjusted)} mL * {ventilation_adjustment_factor}")
    print(f"Result: {int(maintenance_volume_adjusted)} mL\n")

    # --- Calculation 3: Total Deficit Replacement ---
    print("--- 3. Total Deficit Replacement Volume ---")
    # 1 kg weight loss ~ 1 L (1000 mL) fluid loss
    deficit_volume = weight_kg * dehydration_percent * 1000
    print(f"The calculation is: {weight_kg} kg * {dehydration_percent:.2f} (10% dehydration) * 1000 mL/kg")
    print(f"Result: {int(deficit_volume)} mL\n")

    # --- Final Output ---
    final_answer = f"<<<{int(resuscitation_volume)},{int(maintenance_volume_adjusted)},{int(deficit_volume)}>>>"
    print(f"Final calculated volumes (Resuscitation, Maintenance, Deficit):")
    print(final_answer)

if __name__ == '__main__':
    calculate_pediatric_fluids()