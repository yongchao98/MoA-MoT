import math

def calculate_fluid_volumes():
    """
    Calculates three-phase fluid management volumes for a pediatric patient.
    """
    # Patient Data
    weight_kg = 12
    dehydration_percentage = 0.10  # 10%

    # --- Phase 1: Initial Resuscitation ---
    # Formula: 30 mL/kg
    resuscitation_bolus_volume = 30 * weight_kg

    # --- Phase 2: Daily Maintenance Fluids (Holliday-Segar Method) ---
    # Formula:
    # 100 mL/kg for the first 10 kg
    # 50 mL/kg for the next 10 kg (11-20 kg)
    # 20 mL/kg for weight > 20 kg
    if weight_kg <= 10:
        maintenance_volume = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)
    
    # The notes on antibiotics, milk, and mechanical ventilation are important clinical
    # considerations for adjusting the *infusion rate* or *total intake*, but they do not change
    # the fundamental calculation of the three specific volumes requested by the prompt.

    # --- Phase 3: Deficit Replacement ---
    # Formula: 10% of body weight (1 kg of weight loss = 1000 mL of fluid deficit)
    deficit_kg = weight_kg * dehydration_percentage
    deficit_replacement_volume = deficit_kg * 1000

    # Print the final calculated volumes, formatted as requested.
    # The volumes are converted to integers as fluid calculations are typically rounded.
    print(f"{int(resuscitation_bolus_volume)},{int(maintenance_volume)},{int(deficit_replacement_volume)}")

calculate_fluid_volumes()