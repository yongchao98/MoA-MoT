def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen.
    """
    weight_kg = 12
    resuscitation_rate_ml_per_kg = 30
    dehydration_percentage = 0.10  # 10%

    # --- Phase 1: Resuscitation Volume ---
    resuscitation_volume = resuscitation_rate_ml_per_kg * weight_kg
    print("1. Initial Resuscitation Volume:")
    print(f"{resuscitation_rate_ml_per_kg} mL/kg * {weight_kg} kg = {int(resuscitation_volume)} mL")
    print("-" * 30)

    # --- Phase 2: Daily Maintenance Volume (Holliday-Segar Method) ---
    maintenance_volume = 0
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    else:
        first_10kg_fluid = 10 * 100
        remaining_weight = weight_kg - 10
        next_fluid = remaining_weight * 50
        maintenance_volume = first_10kg_fluid + next_fluid
    
    print("2. Daily Maintenance Fluid Volume:")
    print("Calculation based on Holliday-Segar method:")
    print(f"(10 kg * 100 mL/kg) + ({weight_kg - 10} kg * 50 mL/kg) = {10 * 100} + {int((weight_kg - 10) * 50)} = {int(maintenance_volume)} mL/day")
    print("-" * 30)

    # --- Phase 3: Total Deficit Replacement Volume ---
    # Deficit (in L) = weight (kg) * % dehydration
    # Assuming 1 kg of weight loss = 1 L of fluid deficit
    deficit_volume_ml = weight_kg * dehydration_percentage * 1000
    print("3. Total Deficit Replacement Fluid Volume:")
    print(f"{weight_kg} kg * {int(dehydration_percentage * 100)}% dehydration = {int(deficit_volume_ml)} mL")
    print("-" * 30)

    # --- Final Formatted Answer ---
    final_answer = f"<<<{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume_ml)}>>>"
    return final_answer

# Execute the calculation and print the final result
final_result = calculate_pediatric_fluids()
print(f"Final Answer: {final_result}")
