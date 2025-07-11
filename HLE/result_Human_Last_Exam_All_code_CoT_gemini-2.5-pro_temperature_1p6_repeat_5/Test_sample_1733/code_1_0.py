import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a clinical scenario.
    """
    # Patient data and prescriptions
    weight_kg = 12
    resuscitation_bolus_per_kg = 30  # mL/kg
    dehydration_percentage = 0.10  # 10%

    # --- Calculation 1: Initial Resuscitation Volume ---
    resuscitation_volume = weight_kg * resuscitation_bolus_per_kg
    print("1. Initial Resuscitation Volume Calculation:")
    print(f"   {weight_kg} kg * {resuscitation_bolus_per_kg} mL/kg = {int(resuscitation_volume)} mL")
    print("-" * 30)

    # --- Calculation 2: Daily Maintenance Fluid Volume (Holliday-Segar Method) ---
    maintenance_fluid_24hr = 0
    weight_remaining = weight_kg
    
    # First 10 kg
    if weight_remaining > 0:
        first_10_kg_calc = min(weight_remaining, 10)
        maintenance_first_10 = first_10_kg_calc * 100
        maintenance_fluid_24hr += maintenance_first_10
        weight_remaining -= first_10_kg_calc

    # Next 10 kg (10.1 to 20 kg)
    if weight_remaining > 0:
        next_10_kg_calc = min(weight_remaining, 10)
        maintenance_next_10 = next_10_kg_calc * 50
        maintenance_fluid_24hr += maintenance_next_10
        weight_remaining -= next_10_kg_calc

    print("2. Daily Maintenance Fluid Volume Calculation (Holliday-Segar):")
    # Using 2kg for the second part as the child is 12kg
    print(f"   (10 kg * 100 mL/kg) + (2 kg * 50 mL/kg) = 1000 + 100 = {int(maintenance_fluid_24hr)} mL")
    print("-" * 30)
    
    # --- Calculation 3: Total Deficit Replacement Volume ---
    # Deficit is 10% of body weight. 1kg of weight loss = 1000 mL fluid deficit.
    deficit_volume = weight_kg * dehydration_percentage * 1000
    print("3. Total Deficit Replacement Volume Calculation:")
    print(f"   {weight_kg} kg * {int(dehydration_percentage * 100)}% dehydration * 1000 mL/kg = {int(deficit_volume)} mL")
    print("-" * 30)

    # --- Final Formatted Answer ---
    final_answer = f"{int(resuscitation_volume)},{int(maintenance_fluid_24hr)},{int(deficit_volume)}"
    print(f"Final calculated values (Resuscitation, Maintenance, Deficit):")
    print(final_answer)

    # Required output format
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    calculate_pediatric_fluids()