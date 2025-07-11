import math

def calculate_fluid_volumes():
    """
    Calculates fluid volumes for a pediatric patient based on a three-phase regimen.
    """
    # Patient information
    weight_kg = 12
    dehydration_percentage = 0.10  # 10%

    # --- Phase 1: Initial Resuscitation ---
    # Formula: 30 mL/kg bolus
    resuscitation_bolus_per_kg = 30
    resuscitation_volume = weight_kg * resuscitation_bolus_per_kg
    print(f"Phase 1: The total resuscitation fluid volume is {int(resuscitation_volume)} mL.")

    # --- Phase 2: Daily Maintenance Fluids (Holliday-Segar) ---
    # Formula:
    # First 10 kg: 100 mL/kg
    # Next 10 kg (10.1-20 kg): 50 mL/kg
    # Above 20 kg: 20 mL/kg
    maintenance_volume = 0
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_for_first_10kg = 10 * 100
        maintenance_for_remainder = (weight_kg - 10) * 50
        maintenance_volume = maintenance_for_first_10kg + maintenance_for_remainder
    else: # weight > 20
        maintenance_for_first_10kg = 10 * 100
        maintenance_for_next_10kg = 10 * 50
        maintenance_for_remainder = (weight_kg - 20) * 20
        maintenance_volume = maintenance_for_first_10kg + maintenance_for_next_10kg + maintenance_for_remainder
    print(f"Phase 2: The daily maintenance fluid volume is {int(maintenance_volume)} mL.")


    # --- Phase 3: Deficit Replacement ---
    # Formula: Deficit (L) = % dehydration * weight (kg)
    # 1 kg fluid deficit is assumed to be 1000 mL
    deficit_volume = weight_kg * dehydration_percentage * 1000
    print(f"Phase 3: The total deficit replacement fluid volume is {int(deficit_volume)} mL.")

    # --- Final Formatted Answer ---
    print("\nFormatted Answer:")
    # Print the three numbers for the final answer format
    print(f"{int(resuscitation_volume)},{int(maintenance_volume)},{int(deficit_volume)}")

calculate_fluid_volumes()