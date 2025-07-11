import math

def solve_fluid_rate():
    """
    Calculates the maintenance fluid rate for a pediatric patient.
    """
    # Patient and drug information
    weight_kg = 22
    enteral_feeding_ml_day = 500
    drug_dose_mg_per_m2 = 25
    bsa_m2 = 0.8
    admin_concentration_mg_ml = 1.0

    # Step 1: Calculate total daily maintenance fluid using Holliday-Segar method
    print("Step 1: Calculate Total Daily Maintenance Fluid (Holliday-Segar)")
    first_10kg_fluid = 10 * 100
    next_10kg_fluid = 10 * 50
    remaining_weight = weight_kg - 20
    remaining_fluid = remaining_weight * 20
    total_maintenance_fluid_ml_day = first_10kg_fluid + next_10kg_fluid + remaining_fluid
    print(f"(10 kg * 100 ml/kg) + (10 kg * 50 ml/kg) + ({remaining_weight} kg * 20 ml/kg) = {total_maintenance_fluid_ml_day} ml/day")
    print("-" * 50)

    # Step 2: Calculate daily chemotherapy drug dose in mg
    print("Step 2: Calculate Daily Chemotherapy Dose")
    daily_drug_dose_mg = drug_dose_mg_per_m2 * bsa_m2
    print(f"{drug_dose_mg_per_m2} mg/m² * {bsa_m2} m² = {daily_drug_dose_mg} mg/day")
    print("-" * 50)

    # Step 3: Calculate daily volume of chemotherapy infusion
    print("Step 3: Calculate Daily Chemotherapy Infusion Volume")
    chemo_infusion_volume_ml_day = daily_drug_dose_mg / admin_concentration_mg_ml
    print(f"{daily_drug_dose_mg} mg / {admin_concentration_mg_ml} mg/ml = {chemo_infusion_volume_ml_day} ml/day")
    print("-" * 50)

    # Step 4: Calculate total fluid from other sources (enteral feeding + chemo)
    print("Step 4: Calculate Total Fluid from Other Sources")
    total_other_fluids_ml_day = enteral_feeding_ml_day + chemo_infusion_volume_ml_day
    print(f"{enteral_feeding_ml_day} ml (milk) + {chemo_infusion_volume_ml_day} ml (chemo) = {total_other_fluids_ml_day} ml/day")
    print("-" * 50)

    # Step 5: Calculate the required daily volume of maintenance fluid
    print("Step 5: Calculate Required Daily Maintenance Fluid Volume")
    required_maintenance_fluid_ml_day = total_maintenance_fluid_ml_day - total_other_fluids_ml_day
    print(f"{total_maintenance_fluid_ml_day} ml (total) - {total_other_fluids_ml_day} ml (other sources) = {required_maintenance_fluid_ml_day} ml/day")
    print("-" * 50)

    # Step 6: Convert daily volume to an hourly rate and round
    print("Step 6: Convert to Hourly Rate and Round")
    maintenance_rate_ml_hr = required_maintenance_fluid_ml_day / 24
    final_rate_rounded = round(maintenance_rate_ml_hr)
    print(f"{required_maintenance_fluid_ml_day} ml / 24 hours = {final_rate_rounded} ml/hr (rounded from {maintenance_rate_ml_hr:.1f})")
    
    return final_rate_rounded

final_answer = solve_fluid_rate()
print(f"\n<<<>>>\n")
print(f"<<<{final_answer}>>>")
