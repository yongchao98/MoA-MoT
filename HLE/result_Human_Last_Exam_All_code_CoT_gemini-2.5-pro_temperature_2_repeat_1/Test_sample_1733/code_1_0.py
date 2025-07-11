import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements based on a three-phase regimen.
    """
    # Patient Data
    weight_kg = 12
    dehydration_percent = 10
    bolus_dose_ml_per_kg = 30

    # Phase 1: Initial Resuscitation Bolus
    # Calculation: weight (kg) * bolus dose (mL/kg)
    # The equation is: 12 kg * 30 mL/kg
    resuscitation_volume_ml = weight_kg * bolus_dose_ml_per_kg

    # Phase 2: Daily Maintenance Fluids (Holliday-Segar Method)
    # For a 12 kg child:
    # 100 mL/kg for the first 10 kg = 100 * 10 = 1000 mL
    # 50 mL/kg for the next 2 kg (12 - 10) = 50 * 2 = 100 mL
    # Total = 1000 + 100 = 1100 mL
    maintenance_fluid_24h_ml = (10 * 100) + ((weight_kg - 10) * 50)

    # Phase 3: Total Deficit Replacement Volume
    # Calculation: % dehydration * weight (kg) * 10
    # The equation is: 10 * 12 kg * 10
    deficit_volume_ml = dehydration_percent * weight_kg * 10

    # Output the final results in the required format
    # The numbers in the final equation are: 360, 1100, 1200
    print(f"{math.trunc(resuscitation_volume_ml)},{math.trunc(maintenance_fluid_24h_ml)},{math.trunc(deficit_volume_ml)}")

calculate_pediatric_fluids()