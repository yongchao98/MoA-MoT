import math

def calculate_pediatric_fluids():
    """
    Calculates fluid requirements for a pediatric patient based on provided clinical data.
    """
    # --- Patient Data ---
    weight_kg = 12
    bolus_rate_ml_per_kg = 30
    dehydration_percent = 0.10  # 10%

    # --- 1. Initial Resuscitation Bolus ---
    # This is an initial fluid bolus to restore circulation quickly.
    bolus_volume_ml = weight_kg * bolus_rate_ml_per_kg
    print("1. Initial Resuscitation Bolus Volume Calculation:")
    print(f"   {weight_kg} kg * {bolus_rate_ml_per_kg} mL/kg = {int(bolus_volume_ml)} mL")
    print("-" * 50)

    # --- 2. Daily Maintenance Fluid (Holliday-Segar Method) ---
    # This formula calculates the fluid needed for a 24-hour period to cover basal metabolic needs.
    maintenance_volume_ml = 0
    # For a 12kg child:
    # First 10kg: 10kg * 100 mL/kg = 1000 mL
    # Remaining 2kg: 2kg * 50 mL/kg = 100 mL
    # Total = 1000 + 100 = 1100 mL
    first_10kg_fluid = 10 * 100
    remaining_weight = weight_kg - 10
    next_10kg_fluid = remaining_weight * 50
    maintenance_volume_ml = first_10kg_fluid + next_10kg_fluid
    
    print("2. Daily Maintenance Fluid Volume Calculation (Holliday-Segar):")
    print(f"   For the first 10 kg: 10 kg * 100 mL/kg = {int(first_10kg_fluid)} mL")
    print(f"   For the remaining {int(remaining_weight)} kg: {int(remaining_weight)} kg * 50 mL/kg = {int(next_10kg_fluid)} mL")
    print(f"   Total Maintenance Fluid = {int(first_10kg_fluid)} mL + {int(next_10kg_fluid)} mL = {int(maintenance_volume_ml)} mL")
    print("-" * 50)
    
    # --- 3. Total Deficit Replacement ---
    # This is the total volume of fluid lost due to dehydration.
    # It's calculated as a percentage of the child's body weight. 1 kg of lost weight equals 1000 mL of fluid deficit.
    deficit_volume_ml = weight_kg * dehydration_percent * 1000
    print("3. Total Deficit Replacement Fluid Volume Calculation:")
    print(f"   {weight_kg} kg * {int(dehydration_percent * 100)}% dehydration = {weight_kg * dehydration_percent:.1f} L = {int(deficit_volume_ml)} mL")
    print("-" * 50)

    # --- Final Answer ---
    final_answer = f"{int(bolus_volume_ml)},{int(maintenance_volume_ml)},{int(deficit_volume_ml)}"
    print("\nFinal calculated volumes (Resuscitation, Maintenance, Deficit):")
    print(final_answer)

calculate_pediatric_fluids()