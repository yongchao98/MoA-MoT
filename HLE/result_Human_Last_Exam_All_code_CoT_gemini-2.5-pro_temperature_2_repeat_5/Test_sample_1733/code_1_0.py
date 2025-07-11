import math

def calculate_pediatric_fluids():
    """
    Calculates pediatric fluid requirements for a three-phase regimen.
    """
    # --- Inputs ---
    weight_kg = 12
    resuscitation_dose_ml_kg = 30
    dehydration_percentage = 0.10  # 10%

    print("Patient's weight: {} kg".format(weight_kg))
    print("-" * 30)

    # --- 1. Phase 1: Initial Resuscitation ---
    print("1. Calculating Phase 1: Resuscitation Volume...")
    resuscitation_volume_ml = weight_kg * resuscitation_dose_ml_kg
    print("   Formula: Weight (kg) * Dose (mL/kg)")
    print("   Calculation: {} kg * {} mL/kg = {} mL".format(weight_kg, resuscitation_dose_ml_kg, resuscitation_volume_ml))
    print("-" * 30)

    # --- 2. Phase 2: Daily Maintenance Fluids ---
    print("2. Calculating Phase 2: Daily Maintenance Volume...")
    # Holliday-Segar Method
    if weight_kg <= 10:
        maintenance_volume_unadjusted = weight_kg * 100
    elif weight_kg <= 20:
        maintenance_volume_unadjusted = (10 * 100) + ((weight_kg - 10) * 50)
    else:
        maintenance_volume_unadjusted = (10 * 100) + (10 * 50) + ((weight_kg - 20) * 20)

    print("   Holliday-Segar Calculation (unadjusted):")
    print("   For first 10 kg: 10 kg * 100 mL/kg = 1000 mL")
    print("   For next 2 kg: (12 - 10) kg * 50 mL/kg = 100 mL")
    print("   Total unadjusted maintenance: 1000 + 100 = {} mL/day".format(maintenance_volume_unadjusted))

    # Adjustment for mechanical ventilation (typically 25% reduction)
    ventilation_reduction_factor = 0.25
    maintenance_volume_adjusted = maintenance_volume_unadjusted * (1 - ventilation_reduction_factor)
    print("\n   Note: Patient is mechanically ventilated, reducing insensible losses.")
    print("   Applying a 25% reduction to maintenance fluids.")
    print("   Adjusted Calculation: {} mL * (1 - {}) = {} mL/day".format(maintenance_volume_unadjusted, ventilation_reduction_factor, maintenance_volume_adjusted))
    print("-" * 30)

    # --- 3. Phase 3: Deficit Replacement ---
    print("3. Calculating Phase 3: Total Deficit Replacement Volume...")
    # Deficit is 10% of body weight. 1 kg weight loss ~ 1000 mL fluid deficit.
    deficit_volume_ml = weight_kg * dehydration_percentage * 1000
    print("   Formula: Weight (kg) * Dehydration % * 1000 mL/kg")
    print("   Calculation: {} kg * {} * 1000 mL/kg = {} mL".format(weight_kg, dehydration_percentage, deficit_volume_ml))
    print("-" * 30)

    # --- Final Output ---
    # The final results are rounded to the nearest whole number as is common practice.
    final_resuscitation = int(round(resuscitation_volume_ml))
    final_maintenance = int(round(maintenance_volume_adjusted))
    final_deficit = int(round(deficit_volume_ml))

    # The other fluids (antibiotics, milk) are noted but do not change the calculation for these three specific questions.

    final_answer = "{},{},{}".format(final_resuscitation, final_maintenance, final_deficit)
    print("\nFinal Answer (Resuscitation Volume, Daily Maintenance, Total Deficit):")
    print(final_answer)
    print("<<<{}>>>".format(final_answer))

calculate_pediatric_fluids()