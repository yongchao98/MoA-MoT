def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen.
    """
    # 1. Patient Data
    weight_kg = 12
    bolus_rate_ml_per_kg = 30
    dehydration_percent = 0.10

    # --- Calculation for Phase 1: Initial Resuscitation ---
    resuscitation_volume = bolus_rate_ml_per_kg * weight_kg

    # --- Calculation for Phase 2: Daily Maintenance ---
    # Holliday-Segar Method for a 12 kg child
    # First 10 kg at 100 mL/kg
    maintenance_part1_weight = 10
    maintenance_part1_rate = 100
    maintenance_fluid_part1 = maintenance_part1_weight * maintenance_part1_rate

    # Next part of weight (11-20 kg) at 50 mL/kg
    maintenance_part2_weight = weight_kg - maintenance_part1_weight
    maintenance_part2_rate = 50
    maintenance_fluid_part2 = maintenance_part2_weight * maintenance_part2_rate
    
    daily_maintenance_volume = maintenance_fluid_part1 + maintenance_fluid_part2

    # --- Calculation for Phase 3: Deficit Replacement ---
    # The total deficit is 10% of body weight
    # 1 kg of body weight loss corresponds to 1000 mL of fluid deficit
    total_deficit_volume = dehydration_percent * weight_kg * 1000

    # --- Print Detailed Steps and Results ---
    print("Here is a breakdown of the fluid calculations:")
    print("-" * 40)

    # 1. Print Resuscitation Calculation
    print("1. Initial Resuscitation Phase Volume:")
    print(f"The calculation is a one-time bolus of {bolus_rate_ml_per_kg} mL/kg.")
    print(f"   Equation: {bolus_rate_ml_per_kg} mL/kg * {weight_kg} kg = {int(resuscitation_volume)} mL")
    print("-" * 40)

    # 2. Print Maintenance Calculation
    print("2. Daily Maintenance Fluid Volume (24 hours):")
    print("The calculation uses the Holliday-Segar method:")
    print(f"   For the first {maintenance_part1_weight} kg: {maintenance_part1_weight} kg * {maintenance_part1_rate} mL/kg = {int(maintenance_fluid_part1)} mL")
    print(f"   For the next {maintenance_part2_weight} kg: {maintenance_part2_weight} kg * {maintenance_part2_rate} mL/kg = {int(maintenance_fluid_part2)} mL")
    print(f"   Total Daily Maintenance = {int(maintenance_fluid_part1)} + {int(maintenance_fluid_part2)} = {int(daily_maintenance_volume)} mL")
    print("-" * 40)

    # 3. Print Deficit Calculation
    print("3. Total Deficit Replacement Fluid Volume:")
    print(f"The calculation is based on an estimated {dehydration_percent:.0%} dehydration.")
    print(f"   Equation: {dehydration_percent} * {weight_kg} kg * 1000 mL/kg = {int(total_deficit_volume)} mL")
    print("-" * 40)

    # --- Final formatted output ---
    final_answer = f"{int(resuscitation_volume)},{int(daily_maintenance_volume)},{int(total_deficit_volume)}"
    print(f"\nThe comma-separated values are: {final_answer}")
    
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    calculate_pediatric_fluids()